#!/bin/bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup cmake v3_14_3
genie_xsec     	v3_00_04a -q e19:prof
genie_phyopt    v3_00_04 -q e19:prof
