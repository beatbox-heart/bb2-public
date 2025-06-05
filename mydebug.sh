#!/bin/bash
# autoreconf -fi
make clean
export CC=mpicc
./configure CFLAGS="-g -O0" --prefix=$HOME
make
make install
