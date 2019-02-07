### Eigensolver_gpu
GPU Eigensolver for Quantum ESPRESSO package

###
This library implements a generalized eigensolver for symmetric/hermetian-definite eigenproblems with functionality similar to
the DSYGVD/X or ZHEGVD/X functions available within LAPACK/MAGMA. This solver has less dependencies on CPU computation 
than comparable implementations within MAGMA, which may be of benefit to systems with limited CPU resources or to 
users without access to high-performing CPU LAPACK libraries. 

###
This implementation can be considered as a "proof of concept" and has been written to target the Quantum ESPRESSO
code. As such, this implementation is built only to handle one problem configuration of DSYGVD/X or ZHEGVD/X. Specifically, this
solver computes eigenvalues and associated eigenvectors over a specified integer range for a 
symmetric/hermetian-definite eigenproblem in the following form: 

	A * x = lambda * B * x

where `A` and `B` are symmetric/hermetian-matrices and `B` is positive definite. The solver expects the upper-triangular parts of the 
input `A` and `B` arguments to be populated. This configuration corresponds to calling DSYGVX/ZHEGVX within LAPACK with the configuration 
arguments `ITYPE = 1`, `JOBZ = 'V'`, `RANGE = 'I'`, and `UPLO = 'U'`. 

See comments within `dsygvdx_gpu.F90` or `zhegvdx_gpu.F90` for specific details on usage.

For additional information about the solver with some performance results, see presentation at the following link: (will be added
once available publically on the GTC On-Demand website)

### Building
* Compilation of this library requires the PGI compiler version 17.4 or higher.
* Using the provided `Makefile` will generate a static library object `lib_eigsolve.a` which can included in your
target application. 
* Library requires linking to cuBLAS and cuSOLVER. Use `-Mcuda=cublas,cusolver` flag when linking your application to do this.
* This library also requires linking to a CPU LAPACK library with an implementation of the `zstedc` function.
* If NVTX is enabled with `-DUSE_NVTX` flag, also must link to NVTX. Use `-L${CUDAROOT}/lib64 -lnvToolsExt` flag when linking your application to do this
  where `${CUDAROOT}` is the root directory of your CUDA installation.

An example of using this solver in a program can be found in the `test_driver` subdirectory. This program does a little performance testing
and validation against existing functionality in a linked CPU LAPACK library, cuSOLVER, and MAGMA (if available). 

### License
This code is released under an MIT license which can be found in `LICENSE`. 

### Setup:

ldd ./test_zhegvdx
	linux-vdso.so.1 =>  (0x00007ffe41d99000)
	libnvToolsExt.so.1 => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/2018/cuda/10.0/lib64/libnvToolsExt.so.1 (0x00007fcf75816000)
	libiomp5.so => /mnt/2ndHDD/users/srodriguez/perforce/sw/tools/intel/linux/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64/libiomp5.so (0x00007fcf75472000)
	libmkl_intel_lp64.so => /mnt/2ndHDD/users/srodriguez/perforce/sw/tools/intel/linux/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_intel_lp64.so (0x00007fcf74982000)
	libmkl_core.so => /mnt/2ndHDD/users/srodriguez/perforce/sw/tools/intel/linux/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_core.so (0x00007fcf72c3e000)
	libmkl_intel_thread.so => /mnt/2ndHDD/users/srodriguez/perforce/sw/tools/intel/linux/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_intel_thread.so (0x00007fcf70f45000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007fcf70d28000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fcf70a1f000)
	libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007fcf7081b000)
	libcudafor91.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libcudafor91.so (0x00007fcf7060f000)
	libcudafor.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libcudafor.so (0x00007fcf6e102000)
	libcublas.so.10.0 => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/2018/cuda/10.0/lib64/libcublas.so.10.0 (0x00007fcf69b6c000)
	libcusparse.so.10.0 => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/2018/cuda/10.0/lib64/libcusparse.so.10.0 (0x00007fcf66104000)
	libcurand.so.10.0 => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/2018/cuda/10.0/lib64/libcurand.so.10.0 (0x00007fcf61f9d000)
	libcudaforwrapblas.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libcudaforwrapblas.so (0x00007fcf61d64000)
	libcusolver.so.10.0 => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/2018/cuda/10.0/lib64/libcusolver.so.10.0 (0x00007fcf5967d000)
	libcudart.so.10.0 => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/2018/cuda/10.0/lib64/libcudart.so.10.0 (0x00007fcf59403000)
	libcudafor2.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libcudafor2.so (0x00007fcf59202000)
	libpgf90rtl.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libpgf90rtl.so (0x00007fcf58fdd000)
	libpgf90.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libpgf90.so (0x00007fcf58a10000)
	libpgf90_rpm1.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libpgf90_rpm1.so (0x00007fcf5880e000)
	libpgf902.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libpgf902.so (0x00007fcf585fa000)
	libpgftnrtl.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libpgftnrtl.so (0x00007fcf583c7000)
	libpgmp.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libpgmp.so (0x00007fcf58146000)
	libnuma.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libnuma.so (0x00007fcf57f45000)
	libpgmath.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libpgmath.so (0x00007fcf57b67000)
	libpgc.so => /mnt/2ndHDD/users/srodriguez/PGI/linux86-64/18.10/lib/libpgc.so (0x00007fcf57918000)
	librt.so.1 => /lib/x86_64-linux-gnu/librt.so.1 (0x00007fcf57710000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fcf57346000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fcf57130000)
	libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007fcf56dae000)
	/lib64/ld-linux-x86-64.so.2 (0x00007fcf75a1f000)
