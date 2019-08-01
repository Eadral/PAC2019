#!/bin/bash

MY_HOME=`pwd`

NETCDF="${MY_HOME}/netcdf"
NETCDFINC=${NETCDF}/include
NETCDFLIB=${NETCDF}/lib
cd ./build
rm -f *
ln -sf ../source/* .

mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c netcdf_mod.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c kinds_mod.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c domain_size.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c constants.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c communicate.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c exit_mod.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c blocks.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c broadcast.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c distribution.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c io_types.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c boundary.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c domain.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c gather_scatter.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c io_netcdf.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c io_binary.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c io.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c global_reductions.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c grid.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c prognostic.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c time_management.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c solver_pcg_mod.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c timers.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c initial.f90
mpiifort -mcmodel=medium -qopenmp -I${NETCDFINC} -I. -O3 -xCORE-AVX2 -fp-model source -assume byterecl -ftz  -free -g -c pcg.f90

# pop_9-off

mpiifort -o pcg.x -mcmodel=medium -qopenmp  blocks.o netcdf_mod.o constants.o distribution.o domain.o exit_mod.o grid.o initial.o io_binary.o io.o io_netcdf.o io_types.o kinds_mod.o pcg.o prognostic.o solver_pcg_mod.o domain_size.o boundary.o  broadcast.o communicate.o gather_scatter.o global_reductions.o time_management.o timers.o -L${NETCDFLIB} -lnetcdf

mv pcg.x ../benchmark
