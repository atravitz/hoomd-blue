// Copyright (c) 2009-2018 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause
// License.

// Maintainer: atravitz

#include <math.h>
#include <memory>
#include "hoomd/Index1D.h"
#include "hoomd/ParticleGroup.h"
#include "hoomd/Updater.h"
#include "hoomd/md/NeighborList.h"

#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
/*! \file DynamicBond.h
    \brief Declares a class for computing bond breakage/formation
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#ifndef __DYNAMICBOND_H__
#define __DYNAMICBOND_H__

//! Creates or breaks bonds with a given probability
/*!
 */
class PYBIND11_EXPORT DynamicBond : public Updater
{
public:
    //! Constructs the compute
    DynamicBond(std::shared_ptr<SystemDefinition> sysdef,
                std::shared_ptr<ParticleGroup> group,
                std::shared_ptr<NeighborList> nlist,
                int seed,
                Scalar delta_t,
                int period,
                unsigned int table_width);

    //! Destructor
    virtual ~DynamicBond();

    virtual void setParams(Scalar r_cut, Scalar r_true, std::string bond_type,
                           Scalar delta_G, int n_polymer, int nK);

    virtual void setTable(const std::vector<Scalar> &XB,
                           const std::vector<Scalar> &M,
                           const std::vector<Scalar> &L,
                           Scalar rmin,
                           Scalar rmax);
    //! Take one timestep forward
    virtual void update(unsigned int timestep);



protected:
    std::shared_ptr<ParticleGroup> m_group; //!< Group of particles to operate on
    std::shared_ptr<NeighborList> m_nlist;  //!< neighborlist
    std::shared_ptr<BondData> m_bond_data;  //!< Bond data to use in computing bonds
    int m_seed;                             //!< seed for random number generator
    int period;                             //!< period to create/destroy bonds
    int bond_type;                          //!< bond type to create and break
    Scalar m_r_cut;                         //!< cut off distance for computing bonds
    Scalar m_r_true;                        //!< the "true" radius, i.e. the sticker-colloid energy well
    unsigned int m_table_width;
                                            //!< minimum
    Scalar m_delta_G;                       //!< sticker strength
    Scalar m_delta_t;                       //!< time step from integrator
    std::vector<int> m_nloops;              //!< structure of size N to store number of loops
                                            //!< for each colloid
    int n_polymer;                          //!< number of polymers per colloid
    int m_nK;                               //!< kuhn steps per polymer
};

//! Exports the DynamicBond class to python
void export_DynamicBond(pybind11::module &m);

#endif
