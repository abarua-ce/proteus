#!/bin/bash
PETSC_DIR=${PWD} PETSC_ARCH=arch-darwin-noconda-c-opt ./configure --prefix=${HOME}/petsc-venv --CC=/usr/local/bin/mpicc --CXX=/usr/local/bin/mpicxx --FC=/usr/local/bin/mpif90 --with-numpy --with-petsc4py --with-pthread --download-kokkos --download-kokkos-kernels --download-openblas --download-openblas-use-pthreads --download-superlu --download-superlu_dist --download-mumps --download-eigen --download-tetgen --download-gmsh --download-hdf5 --download-hdf5-configure-arguments="--enable-parallel" --download-libceed --download-metis --download-parmetis --download-opencascade --download-suitesparse --download-scalapack --download-szlib --download-zlib --download-hypre --with-python-exec=${HOME}/petsc-venv/bin/python --download-triangle --download-yaml --download-zoltan --with-debugging=0
