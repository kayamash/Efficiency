#!/bin/sh
run="root -l -b -q efficiencyloop.cpp++ > log.txt"
clean="rm -r efficiencyloop_*"
eval $run
eval $clean
