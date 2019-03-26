#include "external_field_test.cuh"
#include "hoomd/GPUPolymorph.cuh"
#include <stdio.h>

namespace kernel
{
__global__ void test_field(Scalar3* out, const mpcd::ExternalField* field, const Scalar3* pos, const unsigned int N)
    {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    out[idx] = field->evaluate(pos[idx]);
    }
}

cudaError_t test_field(Scalar3* out, const mpcd::ExternalField* field, const Scalar3* pos, const unsigned int N)
    {
    const unsigned int block_size = 32;
    const unsigned int num_blocks = (N + block_size - 1)/block_size;
    kernel::test_field<<<num_blocks,block_size>>>(out, field, pos, N);
    return cudaSuccess;
    }

// if the instantiation is also done here, then everything works fine?
//template mpcd::ConstantForce* hoomd::gpu::device_new(Scalar3);
