#!/bin/sh
run="root -l -b -q efficiencyloop.cpp++"
clean="rm -r efficiencyloop_*"
eval $run
eval $clean
