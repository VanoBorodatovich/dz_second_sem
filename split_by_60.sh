#!/usr/bin/env bash

# ${#foo} - length of string
# i=$((i + 1)) - add to variable
# "${foo:$i:1}" - print from i to i+1
for f in $(cat ENSM_prots)
do
  if [[ $f == ">"* ]]
  then
    echo $f >> file.txt
  else
    for (( i=0; i<${#f}; i=$((i + 60)) ));
    do
      echo "${f:$i:60}" >> file.txt
    done
  fi
done