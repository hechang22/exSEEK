#!/bin/bash

#### cd ~/hc-22
#### ln -s /data/2022-bioinfo-shared/data/quiz-I/exRNA-long .

# 配置featurecount与py脚本依赖的环境变量
export PATH=/data/2022-bioinfo-shared/softwares/subread/subread-2.0.3-Linux-x86_64/bin:$PATH
export PATH=/data/2022-bioinfo-shared/softwares/miniconda3/envs/python-env/bin:$PATH

# 建立一个id文件为merge做准备
touch ./long/sample_ids.txt

# 根据library属性进行featurecount
awk 'NR > 1 {print $1,$5}' ./exRNA-long/metadata.txt | while read line 
do 
    if [[ -n "$(echo $line | grep "forward")" ]];then
        di="1"
    else
        di="0"
    fi
    sample_id=`echo "$line" | cut -d " " -f 1`
    echo $sample_id >> ./long/sample_ids.txt

    #估计假期没人用教学集群所以稍稍提高了线程数
    featureCounts --countReadPairs -O -M -s $di \
    -p -t exon -g gene_id -a ./exRNA-long/genome/gff/gencode.v38.annotation.gff3 \
    -o ./long/counts/${sample_id}.txt -T 24 ./exRNA-long/bam/${sample_id}.bam
done

python3 /data/2022-bioinfo-shared/data/quiz-I/scripts/summarize-table.py --indir ./long/counts \
 --formatter '{}' --sample-ids ./long/sample_ids.txt \
 --row-field 0 --row-name gene_id --first-line 2 --value-field 6 --fillna \
 --output ./long/count.matrix.txt
