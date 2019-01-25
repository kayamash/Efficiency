#!/bin/sh
scp CalcEff.* Efficiency.* efficiencyloop.* kayamash@lxatut01.cern.ch:/home/kayamash/code/Efficiency/
#scp -r -o "ProxyCommand ssh kayamash@lxplus.cern.ch -W %h:%p" CalcEff.* Efficiency.* efficiencyloop.* kayamash@lxatut01.cern.ch:/home/kayamash/code/Efficiency/ 
