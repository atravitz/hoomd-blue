// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: joaander

#ifndef __BOND_EVALUATOR_DYNAMIC_COHENPADE_H__
#define __BOND_EVALUATOR_DYNAMIC_COHENPADE_H__

#ifndef NVCC
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorBondDynamicCohenPade.h
    \brief Defines the bond evaluator class for dynamic Cohen-Pade bonds
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

//! Class for evaluating the Cohen-Pade potential
/*! Evaluates the Cohen-Pade potential in an identical manner to EvaluatorPairLJ for pair potentials. See that
    class for a full motivation and design specifics.

   [TODO] params.x is the K stiffness parameter, and params.y is the nK equilibrium rest length.
*/
class EvaluatorBondDynamicCohenPade
    {
    public:
        //! Define the parameter type used by this pair potential evaluator
        typedef Scalar2 param_type;

        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _params Per type pair parameters of this potential
        */
        DEVICE EvaluatorBondDynamicCohenPade(Scalar _rsq, const param_type& _params)
            : rsq(_rsq),nK(_params.x), r_colloid(_params.y)
            {
            }

        //! Requires a diameter because polymers are bound to the surface
        DEVICE static bool needsDiameter() { return true; }

        //! Accept the diameter values
        /*! \param da Diameter of particle a
            \param db Diameter of particle b
        */
        DEVICE void setDiameter(Scalar da, Scalar db)
            {
            // diameter_a = da;
            // diameter_b = db;
            }

        //! Cohen-Pade doesn't use charge
        DEVICE static bool needsCharge() { return false; }

        //! Accept the optional charge values
        /*! \param qa Charge of particle a
            \param qb Charge of particle b
        */
        DEVICE void setCharge(Scalar qa, Scalar qb) { }

        //! Evaluate the force and energy
        /*! \param force_divr Output parameter to write the computed force divided by r.
            \param bond_eng Output parameter to write the computed bond energy

            \return True if they are evaluated or false if the bond
                    energy is not defined
        */
        DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& bond_eng)
            {
            Scalar r = sqrt(rsq);
            force_divr = r_colloid * (nK / r - Scalar(1.0));

            // if the result is not finite, it is likely because of a division by 0, setting force_divr to 0 will
            // correctly result in a 0 force in this case
            #ifdef NVCC
            if (!isfinite(force_divr))
            #else
            if (!std::isfinite(force_divr))
            #endif
                {
                force_divr = Scalar(0);
                }
            bond_eng = Scalar(0.5) * r_colloid * (nK - r) * (nK - r);

            return true;
            }

        #ifndef NVCC
        //! Get the name of this potential
        /*! \returns The potential name. Must be short and all lowercase, as this is the name energies will be logged as
            via analyze.log.
        */
        static std::string getName()
            {
            return std::string("dynamic_cohenpade");
            }
        #endif

    protected:
        Scalar rsq;        //!< Stored rsq from the constructor
        Scalar nK;          //!< nK parameter
        Scalar r_colloid;        //!< r_colloid parameter
    };


#endif // __BOND_EVALUATOR_DYNAMIC_COHENPADE__
