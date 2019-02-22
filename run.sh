#!/bin/sh
run="nohup root -l -b -q efficiencyloop.cpp++ > /data/data3/zp/kayamash/log/log.txt &"
clean="rm -r efficiencyloop_*"
eval $run
eval $clean
