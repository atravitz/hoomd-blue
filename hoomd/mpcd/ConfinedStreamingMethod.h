// Copyright (c) 2009-2017 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: mphoward

/*!
 * \file mpcd/ConfinedStreamingMethod.h
 * \brief Declaration of mpcd::ConfinedStreamingMethod
 */

#ifndef MPCD_CONFINED_STREAMING_METHOD_H_
#define MPCD_CONFINED_STREAMING_METHOD_H_

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

#include "SystemData.h"
#include "hoomd/extern/pybind/include/pybind11/pybind11.h"

namespace mpcd
{

//! MPCD streaming method
/*!
 * This method implements the base version of ballistic propagation of MPCD
 * particles.
 */
template<class Geometry>
class ConfinedStreamingMethod : public mpcd::StreamingMethod
    {
    public:
        //! Constructor
        /*!
         * \param sysdata MPCD system data
         * \param cur_timestep Current system timestep
         * \param period Number of timesteps between collisions
         * \param phase Phase shift for periodic updates
         * \param geom Streaming geometry
         */
        ConfinedStreamingMethod(std::shared_ptr<mpcd::SystemData> sysdata,
                                unsigned int cur_timestep,
                                unsigned int period,
                                int phase,
                                std::shared_ptr<Geometry> geom)
        : mpcd::StreamingMethod(sysdata, cur_timestep, period, phase),
          m_geom(geom)
          {
          }

        //! Implementation of the streaming rule
        virtual void stream(unsigned int timestep);

    protected:
        std::shared_ptr<Geometry> m_geom;   //!< Streaming geometry
    };

/*!
 * \param timestep Current time to stream
 */
template<class Geometry>
void ConfinedStreamingMethod<Geometry>::stream(unsigned int timestep)
    {
    if (!shouldStream(timestep)) return;

    if (m_prof) m_prof->push("MPCD stream");

    const BoxDim& box = m_mpcd_sys->getCellList()->getCoverageBox();

    ArrayHandle<Scalar4> h_pos(m_mpcd_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_mpcd_pdata->getVelocities(), access_location::host, access_mode::readwrite);

    for (unsigned int cur_p = 0; cur_p < m_mpcd_pdata->getN(); ++cur_p)
        {
        const Scalar4 postype = h_pos.data[cur_p];
        Scalar3 pos = make_scalar3(postype.x, postype.y, postype.z);
        const unsigned int type = __scalar_as_int(postype.w);

        const Scalar4 vel_cell = h_vel.data[cur_p];
        Scalar3 vel = make_scalar3(vel_cell.x, vel_cell.y, vel_cell.z);

        // propagate the particle to its new position ballistically
        Scalar dt = m_mpcd_dt;
        bool collide = true;
        bool any_collide = false;
        do
            {
            pos += dt * vel;
            collide = m_geom->detectCollision(pos, vel, dt);
            // track if a collision ever occurs to redirect particle velocity
            any_collide |= collide;
            }
        while (dt > 0 && collide);

        // wrap and update the position
        int3 image = make_int3(0,0,0);
        box.wrap(pos, image);

        h_pos.data[cur_p] = make_scalar4(pos.x, pos.y, pos.z, __int_as_scalar(type));
        // write velocities when they change
        if (any_collide)
            h_vel.data[cur_p] = make_scalar4(vel.x, vel.y, vel.z, mpcd::detail::NO_CELL);
        }

    // particles have moved, so the cell cache is no longer valid
    m_mpcd_pdata->invalidateCellCache();
    if (m_prof) m_prof->pop();
    }

namespace detail
{
//! Export mpcd::StreamingMethod to python
/*!
 * \param m Python module to export to
 */
template<class Geometry>
void export_ConfinedStreamingMethod(pybind11::module& m)
    {
    {
    namespace py = pybind11;

    const std::string name = "ConfinedStreamingMethod" + Geometry::getName();

    py::class_<mpcd::ConfinedStreamingMethod<Geometry>, std::shared_ptr<mpcd::ConfinedStreamingMethod<Geometry>>>
        (m, name.c_str(), py::base<mpcd::StreamingMethod>())
        .def(py::init<std::shared_ptr<mpcd::SystemData>, unsigned int, unsigned int, int, std::shared_ptr<Geometry>>());
    }
    }
} // end namespace detail
} // end namespace mpcd
#endif // MPCD_CONFINED_STREAMING_METHOD_H_