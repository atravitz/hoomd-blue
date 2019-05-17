// Copyright (c) 2009-2018 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: ?

#include "hoomd/ParticleGroup.h"
#include "hoomd/Updater.h"
#include "hoomd/Index1D.h"
#include "hoomd/md/NeighborList.h"
#include <memory>
#include <math.h>

#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
/*! \file DynamicBond.h
    \brief Declares a class for computing bond breakage/formation
*/

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

//
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
                 int period);

        //! Destructor
        virtual ~DynamicBond();

        virtual void setParams(Scalar r_cut,
                            std::string bond_type,
                            Scalar delta_G,
                            int n_polymer,
                            int nK);

        //! Take one timestep forward
        virtual void update(unsigned int timestep);

    protected:
        std::shared_ptr<ParticleGroup> m_group;   //!< Group of particles to which the dynamic bonding is applied
        std::shared_ptr<NeighborList> m_nlist;    //!< neighborlist
        std::shared_ptr<BondData> m_bond_data;    //!< Bond data to use in computing bonds
        int m_seed;                               //!< seed for random number generator
        int period;                               //!< period to create/destroy bonds
        int seed;                                 //!< a seed for the random number generator
        int bond_type;
        Scalar m_r_cut;                           //!< cut off distance for computing bonds
        Scalar m_delta_G;                         //!< sticker strength
        std::vector<int> m_nloops;
        int n_polymer;                            //!< number of polymers per colloid
        int nK;                                   //!< kuhn steps per polymer

    // private:
    };

//! Exports the DynamicBond class to python
void export_DynamicBond(pybind11::module& m);

#endif
