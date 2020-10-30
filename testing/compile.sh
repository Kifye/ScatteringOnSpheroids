#!/bin/bash

cd ../parameters
ifort -c constants.f90
ifort -c regime.f90

cd ../spheroidal_functions
#ifort -I../parameters -c vb_prolate.for
#ifort -I../parameters -c vb_oblate.for
ifort -I../parameters -c spheroidal.f90

cd ../pair_integrals
ifort -I../parameters -I../spheroidal_functions -c integrals.f90

cd ../symmetric
ifort -I../parameters -c matrix.f90
ifort -I../parameters -I../spheroidal_functions -I../pair_integrals -c scattering.f90

cd ../testing
ifort -I../parameters -I../symmetric -I../spheroidal_functions -I../pair_integrals -c test_ext_sca_real.f90
ifort test_ext_sca_real.o \
../parameters/regime.o \
../parameters/constants.o \
../spheroidal_functions/spheroidal.o \
../spheroidal_functions/vb_prolate.o \
../spheroidal_functions/vb_oblate.o \
../pair_integrals/integrals.o \
../symmetric/matrix.o \
../symmetric/scattering.o \
-o test_ext_sca_real
