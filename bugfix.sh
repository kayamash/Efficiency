#!/bin/bash
branch="dev_kayamash"
file1="Efficiency.cpp"
file2="Efficiency.chh"
file3="CalcEfficiency.cpp"
file4="CalcEfficiency.chh"
file5="efficiencyloop.cpp"
add="git add "
commit="git commit -m "bug fix""
push="git push origin "

eval $add$file1
eval $add$file2
eval $add$file3
eval $add$file4
eval $add$file5
eval $commit
eval $push$branch
