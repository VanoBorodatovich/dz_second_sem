#!/usr/bin/env bash

for f in $(cat file.txt)
do
  if [[ $f == ">"* ]]
  then
    echo $f >> elif.txt
  else
    echo -n $f >> elif.txt
  fi
done
