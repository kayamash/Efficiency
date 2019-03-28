#!/bin/sh
run="nohup root -l -b -q efficiencyloop.cpp++ > /home/kayamash/log/logdata18_349327.txt &"
clean="rm -r efficiencyloop_*"
eval $run
eval $clean
