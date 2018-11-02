#!/bin/bash
branch="dev_kayamash"
add="git add ."
commit="git commit -a"
insert="i"
comment="bug fix"
save=":wq"
push="git push origin "

eval $add
eval $commit
echo $insert
echo $comment
echo $save
eval $push$branch
