#!/bin/bash

run='./run.sh'
eval $run
#for number in 0 4 8 12 16 20 24;do
for number in 4 8 12 16 20 24;do
sleep 50m
sed -i s/"data18_physics_Main_Ztap.root"/"data18_physics_Main_Ztap${number}.root"/g ./efficiencyloop.cpp
sed -i s/"thmin = 0"/"thmin = ${number}"/g ./efficiencyloop.cpp
sed -i s/log.txt/log${number}.txt/g ./run.sh

eval $run
echo $number
sleep 10m

sed -i s/data18_physics_Main_Ztap${number}.root/data18_physics_Main_Ztap.root/g ./efficiencyloop.cpp
sed -i s/"thmin = ${number}"/"thmin = 0"/g ./efficiencyloop.cpp
sed -i s/log${number}.txt/log.txt/g ./run.sh
rm efficiencyloop_cpp*
done

