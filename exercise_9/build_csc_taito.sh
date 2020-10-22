#!/bin/bash

################################################################
## * This script builds available configurations of QMCPACK   ##
##   on Taito at CSC                   .                      ##
##                                                            ##
## * Execute this script in trunk/                            ##
##   ./config/build_csc_taito.sh                              ##
##                                                            ##
## Last modified: March, 2019                                 ##
################################################################


module purge
module load gcc/5.4.0
module load mkl/11.3.2
module load intelmpi/5.1.3
module load hdf5-par/1.8.18
module load fftw/3.3.6
module load boost/1.63
module load cmake/3.9.0

CMAKE_FLAGS="-DCMAKE_C_COMPILER=mpicc \ 
             -DCMAKE_CXX_COMPILER=mpicxx"

# Configure and build cpu real AoS
echo ""
echo ""
echo "building qmcpack for cpu AoS real for taito"
mkdir -p build_taito_cpu_real
cd build_taito_cpu_real
cmake $CMAKE_FLAGS ..
make -j 12
cd ..
ln -sf ./build_taito_cpu_real/bin/qmcpack ./qmcpack_taito_cpu_real

# Configure and build cpu complex AoS
echo ""
echo ""
echo "building qmcpack for cpu AoS complex for taito"
mkdir -p build_taito_cpu_comp
cd build_taito_cpu_comp
cmake -DQMC_COMPLEX=1 $CMAKE_FLAGS ..
make -j 12 
cd ..
ln -sf ./build_taito_cpu_comp/bin/qmcpack ./qmcpack_taito_cpu_comp


# Configure and build cpu real SoA
echo ""
echo ""
echo "building qmcpack for cpu SoA real for taito"
mkdir -p build_taito_cpu_real_SoA
cd build_taito_cpu_real_SoA
cmake -DENABLE_SOA=1 $CMAKE_FLAGS ..
make -j 12
cd ..
ln -sf ./build_taito_cpu_real_SoA/bin/qmcpack ./qmcpack_taito_cpu_real_SoA

# Configure and build cpu complex SoA
echo ""
echo ""
echo "building qmcpack for cpu SoA complex for taito"
mkdir -p build_taito_cpu_comp_SoA
cd build_taito_cpu_comp_SoA
cmake -DQMC_COMPLEX=1 -DENABLE_SOA=1 $CMAKE_FLAGS ..
make -j 12 
cd ..
ln -sf ./build_taito_cpu_comp_SoA/bin/qmcpack ./qmcpack_taito_cpu_comp_SoA




