// Copyright (c) 2009-2017 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: mphoward

/*!
 * \file mpcd/StreamingGeometry.h
 * \brief Definition of valid MPCD streaming geometries.
 */

#ifndef MPCD_STREAMING_GEOMETRY_H_
#define MPCD_STREAMING_GEOMETRY_H_

#ifndef NVCC
#include <string>
#endif // NVCC

namespace mpcd
{

namespace detail
{

//! Boundary conditions at the surface
enum struct boundary : unsigned char
    {
    no_slip=0,
    slip
    };

class BulkGeometry
    {
    public:
        //! Detect collision between the particle and the boundary
        /*!
         * \param pos Proposed particle position
         * \param vel Proposed particle velocity
         * \param dt Integration time remaining
         *
         * \returns True if a collision occurred, and false otherwise
         *
         * \post The particle position \a pos is moved to the point of reflection, the velocity \a vel is updated
         *       according to the appropriate bounce back rule, and the integration time \a dt is decreased to the
         *       amount of time remaining.
         */
        bool detectCollision(Scalar3& pos, Scalar3& vel, Scalar& dt) const
            {
            dt = Scalar(0);
            return false;
            }

        #ifndef NVCC
        //! Get the unique name of this geometry
        static std::string getName()
            {
            return std::string("Bulk");
            }
        #endif // NVCC
    };

//! Parallel plate (slit) geometry
/*!
 * This class defines the geometry consistent with two infinite parallel plates. When the plates are
 * in relative motion, Couette flow can be generated in the channel. If a uniform body force is applied to
 * the fluid, the parabolic Poiseuille flow profile is created. Both flow profiles require the enforcement
 * of no-slip boundary conditions.
 *
 * The channel geometry is defined by two parameters: the channel half-width \a H, and the velocity of the
 * plates \a V. The total distance between the plates is \f$2H\f$. The plates are stacked in the
 * \a z direction, and are centered about the origin \f$z=0\f$. The upper plate moves in the \f$+x\f$ direction
 * with velocity \a V, and the lower plate moves in the \f$-x\f$ direction with velocity \f$-V\f$. Hence, for
 * no-slip boundary conditions there is a velocity profile:
 *
 * \f[
 *      v_x(z) = \frac{Vz}{H}
 * \f]
 *
 * This gives an effective shear rate \f$\dot\gamma = V/H\f$, and the shear stress is $\sigma_{xz}\f$.
 *
 * The geometry enforces boundary conditions \b only on the MPCD solvent particles. Additional interactions
 * are required with any embedded particles using appropriate wall potentials.
 *
 * The wall boundary conditions can optionally be changed to slip conditions. For these BCs, the previous
 * discussion of the various flow profiles no longer applies.
 */
class SlitGeometry
    {
    public:
        //! Constructor
        /*!
         * \param H Channel half-width
         * \param V Velocity of the wall
         * \param bc Boundary condition at the wall (slip or no-slip)
         */
        SlitGeometry(Scalar H, Scalar V, boundary bc)
            : m_H(H), m_V(V), m_bc(bc)
            { }

        //! Get maximum extent of the geometry
        Scalar3 getMax() const
            {
            return make_scalar3(0,0,m_H);
            }

        //! Get minimum extent of the geometry
        Scalar3 getMin() const
            {
            return make_scalar3(0,0,-m_H);
            }

        //! Detect collision between the particle and the boundary
        /*!
         * \param pos Proposed particle position
         * \param vel Proposed particle velocity
         * \param dt Integration time remaining
         *
         * \returns True if a collision occurred, and false otherwise
         *
         * \post The particle position \a pos is moved to the point of reflection, the velocity \a vel is updated
         *       according to the appropriate bounce back rule, and the integration time \a dt is decreased to the
         *       amount of time remaining.
         */
        bool detectCollision(Scalar3& pos, Scalar3& vel, Scalar& dt) const
            {
            /*
             * Detect if particle has left the box, and try to avoid branching or absolute value calls. The sign used
             * in calculations is +1 if the particle is out-of-bounds in the +z direction, -1 if the particle is
             * out-of-bounds in the -z direction, and 0 otherwise.
             *
             * We intentionally use > / < rather than >= / <= to make sure that spurious collisions do not get detected
             * when a particle is reset to the boundary location. A particle landing exactly on the boundary from the bulk
             * can be immediately reflected on the next streaming step, and so the motion is essentially equivalent up to
             * an epsilon of difference in the channel width.
             */
            const signed char sign = (pos.z > m_H) - (pos.z < -m_H);
            // exit immediately if no collision is found or particle is not moving normal to the wall
            // (since no new collision could have occurred if there is no normal motion)
            if (sign == 0 || vel.z == Scalar(0))
                {
                dt = Scalar(0);
                return false;
                }

            /*
             * Remaining integration time dt is amount of time spent traveling distance out of bounds.
             * If sign = +1, then pos.z > H. If sign = -1, then pos.z < -H, and we need difference in
             * the opposite direction.
             *
             * TODO: if dt < 0, it is a spurious collision. How should it be treated?
             */
            dt = (pos.z - sign*m_H) / vel.z;

            // backtrack the particle for dt to get to point of contact
            pos.x -= vel.x*dt;
            pos.y -= vel.y*dt;
            pos.z = sign*m_H;

            // update velocity according to boundary conditions
            // no-slip requires reflection of the tangential components
            if (m_bc == boundary::no_slip)
                {
                vel.x = -vel.x + Scalar(sign * 2) * m_V;
                vel.y = -vel.y;
                }
            // both slip and no-slip have no penetration of the surface
            vel.z = -vel.z;

            return true;
            }

        #ifndef NVCC
        //! Get the unique name of this geometry
        static std::string getName()
            {
            return std::string("Slit");
            }
        #endif // NVCC

    private:
        Scalar m_H; //!< Half of the channel width
        Scalar m_V; //!< Velocity of the wall
        boundary m_bc;  //!< Boundary condition
    };

} // end namespace detail
} // end namespace mpcd

#endif // MPCD_STREAMING_GEOMETRY_H_