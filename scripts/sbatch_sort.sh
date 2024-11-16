#!/bin/bash
#### cd /WORK/hechang/
#### ln -s /data/2022-bioinfo-shared/data/quiz-I/exRNA-long .

awk 'NR > 1 {print $1,$5}' ./exRNA-long/metadata.txt | while read line 
do 
    sbatch ./scripts/sort.sh $line
    #echo $line
done
