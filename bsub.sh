#!/bin/bash
run="root -l -b -q efficiencyloop.cpp++"
bsub="bsub -q 12h -o ./out.log -e ./err.log "
clean="rm efficiencyloop_*"
eval $bsub$run
eval $clean

