#!/bin/bash
branch="master"
add="git add ."
message="bug fix"
push="git push origin "

eval $add
git commit -m "${message}"
eval $push$branch




