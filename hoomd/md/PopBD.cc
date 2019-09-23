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

// TODO alyssa: make these tables or pass in?
Scalar feneEnergy(Scalar x, int nKuhn)
{
  Scalar UBi;
  UBi = nKuhn * -1.5 * log(1.0 - x * x); // original Elnaz code
  return UBi;
}

PopBD::PopBD(std::shared_ptr<SystemDefinition> sysdef,
                         std::shared_ptr<ParticleGroup> group,
                         std::shared_ptr<NeighborList> nlist,
                         int seed,
                         Scalar delta_t,
                         int period,
                         unsigned int table_width
                         )
    : Updater(sysdef),
      m_group(group),
      m_nlist(nlist),
      m_seed(seed),
      m_r_cut(0.0),
      m_r_true(0.0),
      m_delta_t(delta_t),
      m_table_width(table_width)
{
  m_exec_conf->msg->notice(5) << "Constructing PopBD" << endl;
  int n_particles = m_pdata->getN();
  m_nloops.resize(n_particles);
}

/*! \param r_cut cut off distance for computing bonds'
    \param bond_type type of bond to be formed or broken

    Sets parameters for the popBD updater
*/
void PopBD::setParams(Scalar r_cut, Scalar r_true, std::string bond_type,
int n_polymer)
{
  m_r_cut = r_cut;
  m_r_true = r_true;
  std::fill(m_nloops.begin(), m_nloops.end(), n_polymer);
}

void PopBD::setTable(const std::vector<Scalar> &XB,
                           const std::vector<Scalar> &M,
                           const std::vector<Scalar> &L,
                           Scalar rmin,
                           Scalar rmax)
  {
    // range check on the parameters
    if (rmin < 0 || rmax < 0 || rmax <= rmin)
        {
        m_exec_conf->msg->error()<< "popbd.table:  rmin, rmax (" << rmin << "," << rmax
             << ") is invalid."  << endl;
        throw runtime_error("Error initializing PopBD");
        }

    if (XB.size() != m_table_width || M.size() != m_table_width || L.size() != m_table_width)
        {
        m_exec_conf->msg->error() << "popbd.table: table provided to setTable is not of the correct size" << endl;
        throw runtime_error("Error initializing PopBD");
        }

  }
PopBD::~PopBD()
{
  m_exec_conf->msg->notice(5) << "Destroying PopBD" << endl;
}

void PopBD::update(unsigned int timestep)
{
  assert(m_pdata);
  assert(m_nlist);
  int m_nK = 200; // TODO: make this r_cut


  // start by updating the neighborlist
  m_nlist->compute(timestep);

  // get box dimensions
  const BoxDim &box = m_pdata->getGlobalBox();

  // start the profile for this compute
  if (m_prof)
    m_prof->push("PopBD");

  // access the neighbor list
  ArrayHandle<unsigned int> h_n_neigh(m_nlist->getNNeighArray(), access_location::host, access_mode::read);
  ArrayHandle<unsigned int> h_nlist(m_nlist->getNListArray(), access_location::host, access_mode::read);
  ArrayHandle<unsigned int> h_head_list(m_nlist->getHeadList(), access_location::host, access_mode::read);

  // Access the particle data
  ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
  ArrayHandle<Scalar> h_diameter(m_pdata->getDiameters(), access_location::host, access_mode::read);

  // Access bond data
  m_bond_data = m_sysdef->getBondData();

  ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
  ArrayHandle<unsigned int> h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);

  Scalar r_cut_sq = m_r_cut * m_r_cut;

  // for each particle
  for (int i = 0; i < (int)m_pdata->getN(); i++)
  {
    // Access the GPU bond table for reading
    const Index2D &gpu_table_indexer = this->m_bond_data->getGPUTableIndexer();
    ArrayHandle<BondData::members_t> h_gpu_bondlist(
        this->m_bond_data->getGPUTable(), access_location::host,
        access_mode::read);
    ArrayHandle<unsigned int> h_gpu_n_bonds(
        this->m_bond_data->getNGroupsArray(), access_location::host,
        access_mode::read);

    // Access the CPU bond table for reading
    ArrayHandle<unsigned int> h_bond_tags(
        m_bond_data->getTags(), access_location::host, access_mode::read);

    // initialize the RNG
    detail::Saru saru(i, timestep, m_seed);

    // access the particle's position and type (MEM TRANSFER: 4 scalars)
    Scalar3 pi =
        make_scalar3(h_pos.data[i].x, h_pos.data[i].y, h_pos.data[i].z);

    // loop over all of the neighbors of this particle
    // TODO: make sure all eligible bonding particles are within the search
    // radius
    const unsigned int myHead = h_head_list.data[i];
    const unsigned int size = (unsigned int)h_n_neigh.data[i];
    for (unsigned int k = 0; k < size; k++)
    {
      // access the index (j) of neighbor particle (MEM TRANSFER: 1 scalar)
      unsigned int j = h_nlist.data[myHead + k];
      assert(j < m_pdata->getN() + m_pdata->getNGhosts());

      // calculate dr_ji (MEM TRANSFER: 3 scalars / FLOPS: 3)
      Scalar3 pj =
          make_scalar3(h_pos.data[j].x, h_pos.data[j].y, h_pos.data[j].z);
      Scalar3 dx = pi - pj;

      // apply periodic boundary conditions
      dx = box.minImage(dx);

      // calculate r_ij squared (FLOPS: 5) (center to center dist)
      Scalar rsq = dot(dx, dx);

      if (rsq < r_cut_sq)
      {
        // count the number of bonds between i and j
        // TODO alyssa: make this a multimap?
        int n_bonds = h_gpu_n_bonds.data[i];
        int nbridges_ij = 0;
        for (int bond_idx = 0; bond_idx < n_bonds; bond_idx++)
        {
          group_storage<2> cur_bond =
              h_gpu_bondlist.data[gpu_table_indexer(i, bond_idx)];
          int bonded_idx = cur_bond.idx[0]; // bonded-particle's index
          if (bonded_idx == j)
          {
            nbridges_ij += 1;
          }
        }

        Scalar r = sqrt(rsq);
        Scalar surf_dist = r - 2 * m_r_true;
        assert(surf_dist > 0.0);

        // (1) Compute P_ij, P_ji, and Q_ij
        // TODO: make this an input?
        Scalar chain_extension = surf_dist / m_nK;
        if (chain_extension < 1.0 &&
            chain_extension >
                0.0) // bond(extension)-bond(extension_rC)<deltaG?
        {
          Scalar p0 = m_delta_t; //* L(d);
          Scalar q0 = m_delta_t; //* M(d);

          Scalar p_ij = m_nloops[i] * p0 * pow((1 - p0), m_nloops[i]);
          Scalar p_ji = m_nloops[j] * p0 * pow((1 - p0), m_nloops[j]);
          Scalar q_ij = nbridges_ij * q0 * pow((1 - q0), (nbridges_ij - 1.0));

          // (2) generate random numbers
          Scalar rnd1 = saru.s<Scalar>(0, 1);
          Scalar rnd2 = saru.s<Scalar>(0, 1);
          Scalar rnd3 = saru.s<Scalar>(0, 1);
          Scalar rnd4 = saru.s<Scalar>(0, 1);

          // (3) check to see if a loop on i should form
          // a bridge btwn particles i and j
          if (rnd1 < p_ij && m_nloops[i] >= 1.0)
          {
            m_bond_data->addBondedGroup(Bond(0, h_tag.data[i], h_tag.data[j]));
            // TODO alyssa: insert tags into multimap
            m_nloops[i] -= 1;
          }

          // (4) check to see if a loop on j should form
          // a bridge btwn particlesi and j
          if (rnd2 < p_ji && m_nloops[j] >= 1.0)
          {
            m_bond_data->addBondedGroup(Bond(0, h_tag.data[i], h_tag.data[j]));
            // TODO alyssa: insert tags into multimap
            m_nloops[j] -= 1;
          }

          // (5) check to see if a bond should be broken
          // between particles i and j
          if (rnd3 < q_ij && nbridges_ij >= 1.0)
          {
            // remove one bond between i and j
            // for each of the bonds in the *system*
            const unsigned int size = (unsigned int)m_bond_data->getN();
            for (unsigned int bond_number = 0; bond_number < size;
                 bond_number++) // turn into hashtable look-up
            {
              // look up the tag of each of the particles
              // participating in the bond
              const BondData::members_t bond =
                  m_bond_data->getMembersByIndex(bond_number);
              assert(bond.tag[0] < m_pdata->getN());
              assert(bond.tag[1] < m_pdata->getN());

              // transform a and b into indices into the particle data arrays
              // (MEM TRANSFER: 4 integers)
              unsigned int idx_a = h_rtag.data[bond.tag[0]];
              unsigned int idx_b = h_rtag.data[bond.tag[1]];
              assert(idx_a <= m_pdata->getMaximumTag());
              assert(idx_b <= m_pdata->getMaximumTag());

              // remove bond with index "bond_number" between particles i and j,
              // then exit the loop
              if ((idx_a == i && idx_b == j) || (idx_a == j & idx_b == i))
              {
                m_bond_data->removeBondedGroup(h_bond_tags.data[bond_number]);
                // TODO alyssa: remove tags from multimap
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
        // remove all bonds if bond goes beyond extension
        else if (chain_extension > 1.0)
        {
          // break all bridges between i and j
          // for each of the bonds in the *system*
          const unsigned int size = (unsigned int)m_bond_data->getN();
          for (unsigned int bond_number = 0; bond_number < size;
               bond_number++) // turn into hashtable look-up
          {
            // look up the tag of each of the particles participating in the
            // bond
            const BondData::members_t bond =
                m_bond_data->getMembersByIndex(bond_number);
            assert(bond.tag[0] < m_pdata->getN());
            assert(bond.tag[1] < m_pdata->getN());

            // transform a and b into indices into the particle data arrays
            // (MEM TRANSFER: 4 integers)
            unsigned int idx_a = h_rtag.data[bond.tag[0]];
            unsigned int idx_b = h_rtag.data[bond.tag[1]];
            assert(idx_a <= m_pdata->getMaximumTag());
            assert(idx_b <= m_pdata->getMaximumTag());

            // remove bond with index "bond_number" between particles i and j,
            // the exit the loop
            if ((idx_a == i && idx_b == j) || (idx_a == j & idx_b == i))
            {
              m_bond_data->removeBondedGroup(h_bond_tags.data[bond_number]);
              // TODO alyssa: remove tags from multimap
              Scalar rnd5 = saru.s<Scalar>(0, 1);
              if (rnd5 <= 0.5)
              {
                m_nloops[i] += 1;
              }
              else if (rnd5 > 0.5)
              {
                m_nloops[i] += 1;
              }
            }
          }
        }
      }
    }
  }

  if (m_prof)
    m_prof->pop();
}

void export_PopBD(py::module &m)
{
  py::class_<PopBD, std::shared_ptr<PopBD>>(m, "PopBD",
                                                        py::base<Updater>())
      .def(py::init<std::shared_ptr<SystemDefinition>,
                    std::shared_ptr<ParticleGroup>,
                    std::shared_ptr<NeighborList>, int, Scalar, int, int>())
      .def("setParams", &PopBD::setParams);
}
