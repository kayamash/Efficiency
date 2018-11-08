#!/bin/bash
run="root -l -b -q efficiencyloop.cpp++"
bsub="bsub -q 4h -o ~/log/out.log -e ~/log/err.log "
clean="rm efficiencyloop_*"
eval $bsub$run
eval $clean

