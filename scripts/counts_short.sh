#!bin/bash

#### ln -s /data/2022-bioinfo-shared/data/quiz-I/exRNA-short .

#配置脚本依赖库的环境变量
export PATH=/data/2022-bioinfo-shared/softwares/miniconda3/envs/python-env/bin:$PATH

touch ./short/sample_ids.txt

awk 'NR > 1 {print $1}' ./exRNA-short/metadata.txt | while read line
do
    echo $line >> ./short/sample_ids.txt

    #miRNA feature count
    python3 /data/2022-bioinfo-shared/data/quiz-I/scripts/count-transcripts.py \
        --bam ./exRNA-short/bam/$line/miRNA.bam \
        --counts ./short/counts/miR/${line}.txt  \
        --stats ./short/counts/miR/${line}_s.txt -s forward
    #piRNA feature count
    python3 /data/2022-bioinfo-shared/data/quiz-I/scripts/count-transcripts.py \
        --bam ./exRNA-short/bam/$line/piRNA.bam \
        --counts ./short/counts/piR/${line}.txt  \
        --stats ./short/counts/piR/${line}_s.txt -s forward
done

# merge
conda activate r-env
Rscript ./scripts/short_RNA_merge.r