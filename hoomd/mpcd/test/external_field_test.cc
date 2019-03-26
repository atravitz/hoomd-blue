// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: mphoward

/*!
 * \file test_gpu_polymorph.cc
 * \brief Tests for GPUPolymorph wrapper.
 */

#include "hoomd/ExecutionConfiguration.h"
#include "hoomd/GPUArray.h"

#include "hoomd/GPUPolymorph.h"
#include "external_field_test.cuh"

#include "hoomd/test/upp11_config.h"
HOOMD_UP_MAIN();

void test_external_field(const ExecutionConfiguration::executionMode mode)
    {
    auto exec_conf = std::make_shared<ExecutionConfiguration>(mode);

    // test points
    GPUArray<Scalar3> pos(2, exec_conf);
        {
        ArrayHandle<Scalar3> h_pos(pos, access_location::host, access_mode::overwrite);
        h_pos.data[0] = make_scalar3(1,2,3);
        h_pos.data[1] = make_scalar3(-3,-2,-1);
        }

    // default initialization is empty
    hoomd::GPUPolymorph<mpcd::ExternalField> field(exec_conf);
    field.reset<mpcd::ConstantForce>(make_scalar3(2,2,2));
    UP_ASSERT(field.get(access_location::host) != nullptr);

    // check host evaluation
        {
        const Scalar3 r = make_scalar3(0,0,0);
        const Scalar3 out = field.get(access_location::host)->evaluate(r);
        CHECK_CLOSE(out.x, 2, tol_small);
        CHECK_CLOSE(out.y, 2, tol_small);
        CHECK_CLOSE(out.z, 2, tol_small);
        }


    #ifdef ENABLE_CUDA
    if (exec_conf->isCUDAEnabled())
        {
        GPUArray<Scalar3> out(2, exec_conf);
            {
            ArrayHandle<Scalar3> d_out(out, access_location::device, access_mode::overwrite);
            ArrayHandle<Scalar3> d_pos(pos, access_location::device, access_mode::read);
            test_field(d_out.data, field.get(access_location::device), d_pos.data, 2);
            }
            {
            ArrayHandle<Scalar3> h_out(out, access_location::host, access_mode::read);
            CHECK_CLOSE(h_out.data[0].x, 2, tol_small);
            CHECK_CLOSE(h_out.data[0].y, 2, tol_small);
            CHECK_CLOSE(h_out.data[0].z, 2, tol_small);

            CHECK_CLOSE(h_out.data[1].x, 2, tol_small);
            CHECK_CLOSE(h_out.data[1].y, 2, tol_small);
            CHECK_CLOSE(h_out.data[1].z, 2, tol_small);
            }
        }
    else
        {
        UP_ASSERT(field.get(access_location::device) == nullptr);
        }
    #endif // ENABLE_CUDA
    }

//! Test external field on CPU
UP_TEST( external_field_cpu )
    {
    test_external_field(ExecutionConfiguration::CPU);
    }

#ifdef ENABLE_CUDA
//! Test external field on GPU
UP_TEST( external_field_gpu )
    {
    test_external_field(ExecutionConfiguration::GPU);
    }
#endif // ENABLE_CUDA
