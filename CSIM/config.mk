# Compiler and compilation flags
cc=gcc
cxx=g++
#cflags=-Wall -ansi -pedantic -O3 -fopenmp
cflags=-Wall -ansi -pedantic -O3
# Compiler for ImageMagick-dependent executables. These must be defined separately
# due to linking issues on the Mac.
im_cxx=$(cxx)
im_cflags=$(cflags)
im_lflags=

# MPI compiler
mpicxx=mpicxx -Wno-long-long

# LAPACK flags for dense linear algebra
lp_lflags=-llapack -lblas
