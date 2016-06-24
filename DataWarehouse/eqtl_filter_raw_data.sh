#! /bin/bash

# NOTE: this should be run in the same directory as the raw eQTL files

mkdir ./filtered_sets

for qtl in $( ls *.txt )
do
  echo $file  
  awk '($5 < 0.1)' $file > filtered_sets/$file
done

