// Copyright (c) 2009-2018 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause
// License.

// Maintainer: atravitz

#include "PopBD.h"
#include "hoomd/GPUArray.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/Index1D.h"
#include "hoomd/VectorMath.h"

#include "hoomd/Saru.h"
using namespace hoomd;
namespace py = pybind11;

/*! \file PopBD.cc
    \brief Contains code for the PopBD class
*/

using namespace std;
/*! \param sysdef SystemDefinition containing the ParticleData to compute forces on
    \param group Group of particles on which to apply this constraint
    \param nlist Neighborlist to use
    \param seed Random number generator seed
    \param period Period at which to call the updater
*/
PopBD::PopBD(std::shared_ptr<SystemDefinition> sysdef,
             std::shared_ptr<ParticleGroup> group,
             std::shared_ptr<NeighborList> nlist,
             int seed,
             Scalar delta_t,
             int period,
             unsigned int table_width)
    : Updater(sysdef),
        m_group(group),
        m_nlist(nlist),
        m_seed(seed),
        m_r_cut(0.0),
        m_delta_t(delta_t),
        m_table_width(table_width)
    {
    m_exec_conf->msg->notice(5) << "Constructing PopBD" << endl;

    assert(m_pdata);
    // access the bond data for later use
    m_bond_data = m_sysdef->getBondData();

    // create arrays for counting loops on particles
    int n_particles = m_pdata->getN();
    m_nloops.resize(n_particles);

    // allocate storage for the tables and parameters
    GPUArray<Scalar2> tables(m_table_width, m_bond_data->getNTypes(), m_exec_conf);
    m_tables.swap(tables);
    GPUArray<Scalar4> params(m_bond_data->getNTypes(), m_exec_conf);
    m_params.swap(params);
    assert(!m_tables.isNull());

    // record all initial bonds

    // helper to compute indices
    Index2D table_value(m_tables.getPitch(),m_bond_data->getNTypes());
    m_table_value = table_value;
    }

PopBD::~PopBD()
    {
    m_exec_conf->msg->notice(5) << "Destroying PopBD" << endl;
    }

/*! \param r_cut cut off distance for computing bonds'
    \param bond_type m_type of bond to be formed or broken

    Sets parameters for the popBD updater
*/
void PopBD::setParams(Scalar r_cut, std::string bond_type, int n_polymer)
    {
    m_r_cut = r_cut;
    std::fill(m_nloops.begin(), m_nloops.end(), n_polymer);
    m_type = 0;
    }

void PopBD::setTable(const std::vector<Scalar> &XB,
                     const std::vector<Scalar> &M,
                     const std::vector<Scalar> &L,
                     Scalar rmin,
                     Scalar rmax)
    {
    // access the arrays
    ArrayHandle<Scalar2> h_tables(m_tables, access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_params(m_params, access_location::host, access_mode::readwrite);

    // range check on the parameters
    if (rmin < 0 || rmax < 0 || rmax <= rmin)
        {
        m_exec_conf->msg->error() << "popbd.table:  rmin, rmax (" << rmin << "," << rmax << ") is invalid." << endl;
        throw runtime_error("Error initializing PopBD");
        }

    if (XB.size() != m_table_width || M.size() != m_table_width || L.size() != m_table_width)
        {
        m_exec_conf->msg->error() << "popbd.table: table provided to setTable is not of the correct size" << endl;
        m_exec_conf->msg->error() << XB.size() << " " << M.size() << " " << L.size() << " " << m_table_width << endl;
        throw runtime_error("Error initializing PopBD");
        }

    // fill out the parameters
    h_params.data[m_type].x = rmin;
    h_params.data[m_type].y = rmax;
    h_params.data[m_type].z = (rmax - rmin) / Scalar(m_table_width - 1);

    // fill out the table
    for (unsigned int i = 0; i < m_table_width; i++)
        {
        h_tables.data[m_table_value(i, m_type)].x = M[i];
        h_tables.data[m_table_value(i, m_type)].y = L[i];
        }
    }

void PopBD::update(unsigned int timestep)
    {
    // start the profile for this compute
    if (m_prof) m_prof->push("PopBD");
    assert(m_pdata);
    assert(m_nlist);

    // start by updating the neighborlist
    m_nlist->compute(timestep);

    // get box dimensions
    const BoxDim &box = m_pdata->getGlobalBox();

    // access the neighbor list
    ArrayHandle<unsigned int> h_n_neigh(m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(m_nlist->getHeadList(), access_location::host, access_mode::read);

    // Access the particle data
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle<Scalar> h_diameter(m_pdata->getDiameters(), access_location::host, access_mode::read);

    // Access bond data
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);

    // access the table data
    ArrayHandle<Scalar2> h_tables(m_tables, access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_params(m_params, access_location::host, access_mode::read);

    Scalar r_cut_sq = m_r_cut * m_r_cut;


    // for each particle
    for (int i = 0; i < (int)m_pdata->getN(); i++)
        {

        // Access the GPU bond table for reading
        // const Index2D &gpu_table_indexer = this->m_bond_data->getGPUTableIndexer();
        ArrayHandle<BondData::members_t> h_gpu_bondlist(this->m_bond_data->getGPUTable(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_gpu_n_bonds(this->m_bond_data->getNGroupsArray(), access_location::host, access_mode::read);

        // Access the CPU bond table for reading
        // ArrayHandle<typename BondData::members_t> h_bonds(m_bond_data->getMembersArray(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> h_bond_tags(m_bond_data->getTags(), access_location::host, access_mode::read);

        // initialize the RNG
        detail::Saru saru(i, timestep, m_seed);

        // access the particle's position and type (MEM TRANSFER: 4 scalars)
        Scalar3 pi = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);

        // loop over all of the neighbors of this particle
        const unsigned int myHead = h_head_list.data[i];
        const unsigned int size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int k = 0; k < size; k++)
            {
            // access the index (j) of neighbor particle (MEM TRANSFER: 1 scalar)
            unsigned int j = h_nlist.data[myHead + k];
            assert(j < m_pdata->getN() + m_pdata->getNGhosts());

            // calculate dr_ji (MEM TRANSFER: 3 scalars / FLOPS: 3)
            Scalar3 pj = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
            Scalar3 dx = pi - pj;

            // apply periodic boundary conditions
            dx = box.minImage(dx);

            // calculate r_ij squared (FLOPS: 5) (center to center dist)
            Scalar rsq = dot(dx, dx);

            if (rsq < r_cut_sq)
                {
                Scalar r = sqrt(rsq);

                // compute index into the table and read in values

                // unsigned int m_type = m_bond_data->getTypeByIndex(i);
                Scalar4 params = h_params.data[m_type];
                Scalar rmin = params.x;
                Scalar rmax = params.y;
                Scalar delta_r = params.z;

                // precomputed term
                Scalar value_f = (r - rmin) / delta_r;

                if (r < rmin)
                    {
                    throw runtime_error("gap is too small to compute a probability!");
                    }

                /// Here we use the table!!
                unsigned int value_i = (unsigned int)floor(value_f);
                Scalar2 ML0 = h_tables.data[m_table_value(value_i, m_type)];
                Scalar2 ML1 = h_tables.data[m_table_value(value_i+1, m_type)];

                // unpack the data
                Scalar M0 = ML0.x;
                Scalar M1 = ML1.x;
                Scalar L0 = ML0.y;
                Scalar L1 = ML1.y;

                // compute the linear interpolation coefficient
                Scalar f = value_f - Scalar(value_i);

                // interpolate to get M and L;
                Scalar M = M0 + f * (M1 - M0);
                Scalar L = L0 + f * (L1 - L0);

                // (1) Compute P_ij, P_ji, and Q_ij

                int nbridges_ij = m_nbonds[std::pair<int,int>(i,j)];

                Scalar p0 = m_delta_t * L;
                Scalar q0 = m_delta_t * M;

                Scalar p_ij = m_nloops[i] * p0 * pow((1 - p0), m_nloops[i] - 1.0);
                Scalar p_ji = m_nloops[j] * p0 * pow((1 - p0), m_nloops[j] - 1.0);
                Scalar q_ij = nbridges_ij * q0 * pow((1 - q0), nbridges_ij - 1.0);



                // check that P and Q are reasonable
                if (p_ij < 0 ||p_ji < 0 || q_ij < 0 || p_ij > 1 ||p_ji > 1 || q_ij > 1)
                    {
                    // cout << "p_ij = " << p_ij << "\n";
                    cout << "m_nloops[i] = " << m_nloops[i] << "\n";
                    // cout << "p_ji = " << p_ji << "\n";
                    cout << "m_nloops[j] = " << m_nloops[j] << "\n";
                    // cout << "q_ij = " << q_ij << "\n";
                    cout << "nbridges[i,j] = " << nbridges_ij << "\n";
                    throw runtime_error("p and q must be between 0 and 1! \n");
                    }

                // (2) generate random numbers
                Scalar rnd1 = saru.s<Scalar>(0, 1);
                Scalar rnd2 = saru.s<Scalar>(0, 1);
                Scalar rnd3 = saru.s<Scalar>(0, 1);
                Scalar rnd4 = saru.s<Scalar>(0, 1);

                // (3) check to see if a loop on i should form a bridge btwn particles i and j
                if (rnd1 < p_ij && m_nloops[i] >= 1)
                    {
                    m_bond_data->addBondedGroup(Bond(m_type, h_tag.data[i], h_tag.data[j]));
                    m_nbonds[std::pair<int,int>(i,j)] += 1;

                    m_nloops[i] -= 1;
                    }

                // (4) check to see if a loop on j should form a bridge btwn particlesi and j
                if (rnd2 < p_ji && m_nloops[j] >= 1)
                    {
                    m_bond_data->addBondedGroup(Bond(m_type, h_tag.data[i], h_tag.data[j]));
                    m_nbonds[std::pair<int,int>(i,j)] += 1;

                    m_nloops[j] -= 1;
                    }

                // (5) check to see if a bond should be broken between particles i and j
                if (rnd3 < q_ij && nbridges_ij >= 1)
                    {
                    // remove one bond between i and j
                    // iterate over each of the bonds in the *system*
                    const unsigned int size = (unsigned int)m_bond_data->getN();
                    for (unsigned int bond_number = 0; bond_number < size; bond_number++)
                        {
                        // look up the tag of both of the particles participating in the bond
                        const BondData::members_t bond = m_bond_data->getMembersByIndex(bond_number);
                        assert(bond.tag[0] < m_pdata->getN());
                        assert(bond.tag[1] < m_pdata->getN());

                        // transform a and b into indices into the particle data arrays
                        // (MEM TRANSFER: 4 integers)
                        unsigned int idx_a = h_rtag.data[bond.tag[0]];
                        unsigned int idx_b = h_rtag.data[bond.tag[1]];
                        assert(idx_a <= m_pdata->getMaximumTag());
                        assert(idx_b <= m_pdata->getMaximumTag());

                        if ((idx_a == i && idx_b == j) || (idx_a == j & idx_b == i))
                            {
                            // remove bond with tag "bond_number" between particles i and j, then leave the loop
                            m_bond_data->removeBondedGroup(h_bond_tags.data[bond_number]);
                            m_nbonds[std::pair<int,int>(i,j)] -= 1;
                            break;
                            }
                        }
                    if (rnd4 <= 0.5)
                        {
                        m_nloops[i] += 1;
                        }
                    else if (rnd4 > 0.5)
                        {
                        m_nloops[j] += 1;
                        }
                    }
                }
            }
        }

    if (m_prof) m_prof->pop();
    }

void export_PopBD(py::module &m)
{
  py::class_<PopBD, std::shared_ptr<PopBD>>(m, "PopBD", py::base<Updater>())
      .def(py::init<std::shared_ptr<SystemDefinition>,
                    std::shared_ptr<ParticleGroup>,
                    std::shared_ptr<NeighborList>,
                    int, Scalar, int, int>())
      .def("setParams", &PopBD::setParams)
      .def("setTable", &PopBD::setTable)
      ;
}
