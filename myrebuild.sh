#!/bin/bash
# autoreconf -fi
make clean
export CC=mpicc
./configure CFLAGS="-O3" --prefix=$HOME
make
make install
