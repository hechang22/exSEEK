#!/bin/bash
#SBATCH -J reads_assignment
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --output=./log/sbatch_sort/%j.out
#SBATCH --error=./log/sbatch_sort/%j.err

#py环境
export PATH=/data/2022-bioinfo-shared/softwares/miniconda3/envs/python-env/bin:$PATH

#传入的参数作为$1与$2可用
sample_id=$1
strand=$2
# 按mRNA,lncRNA,snoRNA,snRNA,srpRNA,tRNA,YRNA,intron,pseudogene的优先级把每个fragment assign到特定RNA类型
python3 /data/2022-bioinfo-shared/data/quiz-I/scripts/reads-assignment.py \
--input ./exRNA-long/bam/${sample_id}.bam --strandness ${strand} \
--output ./long/assign/${sample_id}.txt --bed-dir ./exRNA-long/genome/bed
