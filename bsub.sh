#!/bin/bash
run="root -l -b -q efficiencyloop.cpp++"
bsub="bsub -q 4h -o ./outnoMdt.log -e ./errnoMt.log "
clean="rm efficiencyloop_*"
eval $bsub$run
eval $clean

