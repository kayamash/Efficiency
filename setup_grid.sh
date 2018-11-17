#!/bin/sh
setupATLAS
# 2015b
#asetup AtlasProduction 20.1.8.1 64,gcc48,here
# 2015c
#asetup AtlasProduction 20.7.3.8 here
# 2016
#asetup AtlasProduction 20.7.5.3 here
#asetup AtlasProduction 20.7.7.5 here
#asetup AtlasDerivation 20.7.7.5 here
#asetup AthAnalysisBase 2.4.6 here

#setupATLAS
lsetup emi  
lsetup panda
#voms-proxy-init -voms atlas:/atlas/jp -valid 96:00

#asetup AtlasProduction 20.7.7.5 here
#period C
#asetup AtlasProduction 20.7.6.5 here
#period D
#asetup AtlasProduction 20.7.6.7 here
#period K
#asetup AtlasProduction 20.7.8.3 here
#asetup AtlasProduction 20.7.8.2 here
# broken
#asetup AtlasProduction 20.7.6.2 here
#asetup AtlasProduction 21.0.13.1 here 
asetup Athena,21.0.77 here
#asetup AtlasOffline, 21.0.14 here
#asetup AtlasOffline, 21.0.16 here
