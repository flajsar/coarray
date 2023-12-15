#!/bin/sh

caf -c constants_parameters.f90 -O3 -ffree-line-length-none

caf -c basic_functions.f90 -O3 -ffree-line-length-none

caf -c motion.f90 -O3 -ffree-line-length-none

caf -c start.f90 -O3 -ffree-line-length-none

caf -c results.f90 -O3 -ffree-line-length-none

caf -c termo_barostats.f90 -O3 -ffree-line-length-none

caf -c interactions.f90 -O3 -ffree-line-length-none

caf -c simulation_step.f90 -O3 -ffree-line-length-none

caf -c Coarray.f90 -O3 -ffree-line-length-none

caf *.o -o coarray_release
