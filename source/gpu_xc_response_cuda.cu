#include <cuda_runtime.h>

namespace {

__global__ void xc_response_kernel(
    int nbasis,
    int nstate,
    const double* density,
    const double* kernel,
    double* response) {
  const int total = nbasis * nstate;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < total) {
    // Conservative scaffold: one fused multiply-add-like contribution per
    // packed response slot.  Full XC quadrature/cache integration will replace
    // this smoke path once the Fortran response arrays are wired.
    response[idx] += density[idx] * kernel[idx];
  }
}

}  // namespace

extern "C" int oqp_gpu_xc_response_contract(
    int nbasis,
    int nstate,
    const double* density,
    const double* kernel,
    double* response) {
  if (nbasis <= 0 || nstate <= 0 || density == nullptr || kernel == nullptr || response == nullptr) {
    return 2;
  }

  const int total = nbasis * nstate;
  const int block = 128;
  const int grid = (total + block - 1) / block;
  xc_response_kernel<<<grid, block>>>(nbasis, nstate, density, kernel, response);
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    return static_cast<int>(err);
  }
  err = cudaDeviceSynchronize();
  if (err != cudaSuccess) {
    return static_cast<int>(err);
  }
  return 0;
}
