// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: mphoward

#ifndef HOOMD_MPCD_TEST_EXTERNAL_FIELD_TEST_CUH_
#define HOOMD_MPCD_TEST_EXTERNAL_FIELD_TEST_CUH_

#include "hoomd/mpcd/ExternalField.h"
#include "hoomd/HOOMDMath.h"

cudaError_t test_field(Scalar3* out, const mpcd::ExternalField* field, const Scalar3* pos, const unsigned int N);

#endif // HOOMD_MPCD_TEST_EXTERNAL_FIELD_TEST_CUH_
