#!/bin/bash
branch="dev_kayamash"
add="git add ."
commit="git commit -m "
message="bug fix"
push="git push origin "

#eval $add
#eval $commit"${message}"
#eval $push$branch

git add .
git commit -m "${message}"
eval $push$branch




