#!/bin/bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup cmake v3_14_3
setup genie v3_00_06i -q e19:prof

cd /Package
mkdir build
cd build
cmake -GNinja -DCMAKE_CXX_FLAGS="-fdiagnostics-color=always" .. && \
ninja && \
ninja install
