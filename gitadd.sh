#!/bin/bash
branch="dev_kayamash"
file1="efficiencyloop.cpp"
file2="Efficiency.cpp"
file3="Efficiency.chh"
file4="CalcEff.cpp"
file5="CalcEff.chh"
file7="run.sh"
file8="bsub.sh"
file9="gitadd.sh"
file10="scp.sh"
file11="check.sh"
add="git add "
message="bug fix"
push="git push origin "

#eval $add$file1
#eval $add$file2
#eval $add$file3
#eval $add$file4
#eval $add$file5
#eval $add$file7
#eval $add$file8
#eval $add$file9
#eval $add$file10
#eval $add$file11
git add -A
git commit -m "${message}"
eval $push$branch




