#!/bin/bash
branch="dev_kayamash"
add="git add ."
commit="git commit -a --allow-emply-message -m ''"
push="git push origin "

eval $add
eval $commit
eval $push$branch
