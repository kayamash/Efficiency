#!/bin/bash
run="root -l -b -q efficiencyloop.cpp++"
bsub="bsub -q 4h -o ~/log/outnewClacnoselection.log -e ~/log/errnewCalnocselection.log "
clean="rm efficiencyloop_*"
eval $bsub$run
eval $clean

