// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: mphoward

/*!
 * \file mpcd/ConfinedStreamingMethodGPU.cu
 * \brief Defines GPU functions and kernels used by mpcd::ConfinedStreamingMethodGPU
 */

#include "ConfinedStreamingMethodGPU.cuh"
#include "StreamingGeometry.h"

#include "ExternalField.h"
#include "hoomd/GPUPolymorph.cuh"

// explicitly instantiating this here seems to fix the illegal instruction error
// this error also doesn't show up when compiled in debug mode
// is this a compiler bug, or are we doing something illegal?
TEMPLATE_DEVICE_NEW_FIELDS;

namespace mpcd
{
namespace gpu
{

//! Template instantiation of bulk geometry streaming
template cudaError_t confined_stream<mpcd::detail::BulkGeometry>
    (const stream_args_t& args, const mpcd::detail::BulkGeometry& geom);

//! Template instantiation of slit geometry streaming
template cudaError_t confined_stream<mpcd::detail::SlitGeometry>
    (const stream_args_t& args, const mpcd::detail::SlitGeometry& geom);

} // end namespace gpu
} // end namespace mpcd
