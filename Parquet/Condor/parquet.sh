#!/bin/bash

ls

export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/sft.cern.ch/lcg/views/LCG_97/x86_64-centos7-gcc8-opt/setup.sh

python2 ./convert_AA.py

ls

xrdcp test_Feb26_2023_AA_0p1GeV.HToAA_MA-0p1GeV.parquet.0 root://cmseos.fnal.gov//store/user/kyungmip/CONDOR_TEST/

rm *.parquet.0

ls
