#!/bin/bash
run="root -l -b -q efficiencyloop.cpp++"
bsub="bsub -q 4h -o ./out.log -e ./err.log "
clean="rm efficiencyloop_*"
reset="rm *.log"
eval $reset
eval $bsub$run
eval $clean

