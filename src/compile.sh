#!/bin/bash
# Compile the RELI binary
# (try the 'Makefile' at the top level of the repo first)

g++ RELI.cpp RELI_impl.h RELI_impl.cpp  \
    -o RELI -std=c++11 -O3 -w \
    -lgsl -lgslcblas 
