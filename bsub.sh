#!/bin/bash
run="root -l -b -q efficiencyloop.cpp++"
bsub="bsub -q 1d -o ./out.log -e ./err.log "
clean="rm efficiencyloop_*"
eval $bsub$run
eval $clean

