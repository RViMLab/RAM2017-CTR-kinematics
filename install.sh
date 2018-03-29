#!/bin/bash

# This file is part of CTR, a kinematics library for concentric tube robots
#
# Copyright (C) 2017 Konrad Leibrandt <konrad.lei@gmx.de>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


CURRENT_DIR="$(pwd)"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INSTALL_DIR="$SCRIPT_DIR/install"

SCR='\033[0;31m'
SCG='\033[0;32m'
SCN='\033[0m'

cd $SCRIPT_DIR
mkdir -p ctr_external

echo -e "${SCG}Downloading Dependencies: ${SCN}"

cd ctr_external

############ ERL ############
if [ -d "erl" ]; then
   echo -e "${SCG}Erl exists.${SCN}"
else
   echo -e "${SCR}Erl downloading ...${SCN}"
   git clone https://gitlab.com/ggras/Erl.git
   ln -s Erl erl
fi
############ EIGEN 3.3.3 ############
if [ -d "eigen" ]; then
   echo -e "${SCG}Eigen exists.${SCN}"
else
   echo -e "${SCR}Eigen downloading ...${SCN}"
	wget http://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz
	tar -xzf 3.3.3.tar.gz
	mv eigen-eigen* Eigen-3.3.3
	ln -s Eigen-3.3.3 eigen
fi
############ BOOST 1.65.1 ############
if [ -d "boost" ]; then
   echo -e "${SCG}Boost exists.${SCN}"
else
   echo -e "${SCR}Boost downloading ...${SCN}"
   wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz
   tar -xzf boost_1_65_1.tar.gz
   ln -s boost_1_65_1 boost
fi
#######################################################
cd $SCRIPT_DIR
if [ ! -d "build" ]; then
   echo -e "${SCR}Create build directory.${SCN}"
   mkdir -p build
fi
#######################################################
cd build
#######################################################
echo -e "${SCR}Running cmake ...${SCN}"
cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR $SCRIPT_DIR
#cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DCMAKE_C_COMPILER=/usr/local/gcc-6.1.0/bin/gcc -DCMAKE_CXX_COMPILER=/usr/local/gcc-6.1.0/bin/g++ $SCRIPT_DIR
#######################################################
if [ $? -neq 0 ]; then
   exit 1
else
   echo -e "${SCG}Making concentric tube kinematics library ...${SCN}"
fi
make clean
make -j$(nproc)
if [ $? -ne 0 ]; then
   echo -e "${SCR}Make failed.${SCN}"
   exit 1
else
   echo -e "${SCG}Make successful.${SCN}"
fi
#######################################################
make install
#######################################################
export LD_LIBRARY_PATH=$INSTALL_DIR/lib:$LD_LIBRARY_PATH
cd "$INSTALL_DIR/bin"
#######################################################
echo -e "${SCG}Testing simple binary.${SCN}"
./ctr_kinematics_g_test.bin 256 "$SCRIPT_DIR/ctr_resources/ctr_sec3_tub4.xml"
#######################################################
if [ $? -ne 0 ]; then
   echo -e "${SCR}CTR test failed.${SCN}"
   exit 1
else
   echo -e "${SCG}CTR test  successful.${SCN}"
fi
exit 0
