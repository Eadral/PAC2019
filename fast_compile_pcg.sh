#!/bin/bash

MY_HOME=`pwd`

NETCDF="${MY_HOME}/netcdf"
NETCDFINC=${NETCDF}/include
NETCDFLIB=${NETCDF}/lib
cd ./build
#rm -f *
ln -sf ../source/* .

OPTIONS="-O3 -xHost -ipo -parallel -par-affinity=granularity=fine,proclist=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],explicit -no-prec-div -funroll-all-loops -fp-model fast=2 -qopenmp -falign-loops -assume byterecl -ftz  -free -g"

#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c netcdf_mod.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c kinds_mod.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c domain_size.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c constants.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c communicate.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c exit_mod.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c blocks.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c broadcast.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c distribution.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c io_types.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c boundary.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c domain.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c gather_scatter.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c io_netcdf.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c io_binary.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c io.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c global_reductions.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c grid.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c prognostic.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c time_management.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c solver_pcg_mod.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c timers.f90
#mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c initial.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. ${OPTIONS} -c pcg.f90

# pop_9-off

mpiifort ${OPTIONS} -o pcg.x -mcmodel=medium -qopenmp  blocks.o netcdf_mod.o constants.o distribution.o domain.o exit_mod.o grid.o initial.o io_binary.o io.o io_netcdf.o io_types.o kinds_mod.o pcg.o prognostic.o solver_pcg_mod.o domain_size.o boundary.o  broadcast.o communicate.o gather_scatter.o global_reductions.o time_management.o timers.o -L${NETCDFLIB} -lnetcdf

mv pcg.x ../benchmark
