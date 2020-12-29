#!/bin/env zsh -e

cd ..
dir=(`pwd`)

echo 'SRAの名前を入力してください(.sraの前)'
read accession

echo 'リファレンスゲノムを入力して下さい[mm9,mm10]'
read ref

echo 'クラスター数(k-means)'
read cluster

data=(`grep ${accession} data.csv`)

arr=( `echo ${data} | tr -s ',' ' '`)
sra=(${arr[1]})
sample=(${arr[2]})
seq=(${arr[3]})
ends=(${arr[4]})
reserch=(${arr[5]})
spike=(${arr[6]})

mkdir ${dir}/${reserch}/${sample}_${sra}/${ref}/deeptools/individual
cd ${dir}/${reserch}/${sample}_${sra}/${ref}/deeptools/individual

len=(`pytho3.7 /Users/shigenseigyo/Desktop/reference/effectivegenome/x-removed.py ${ref}`)

bamCoverage --bam ${dir}/${reserch}/${sample}_${sra}/${ref}/bowtie2/${sra}.drm.bam -o individual_${sra}.bw \
    --binSize 10 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize ${len} \
    --ignoreForNormalization chrX \
    --ignoreDuplicates \
    --extendReads

computeMatrix scale-regions --skipZeros -S individual_${sra}.bw -R /Users/shigenseigyo/Desktop/reference/Mus_musculus/UCSC/${ref}/Annotation/Genes/genes.bed -b 1000 -o ${sra}_matrix.mat.gz --outFileNameMatrix IndividualValues.tab -p 8

plotHeatmap --colorMap cividis -m ${sra}_matrix.mat.gz -out ${sra}_Heatmap.png
plotHeatmap --kmeans ${cluster} --colorMap cividis --outFileSortedRegions cluster_HeatmapsortedRegions.bed -m ${sra}_matrix.mat.gz -out ${sra}_culusterHeatmap.png

plotProfile -m ${sra}_matrix.mat.gz -out ${sra}_Profile.png --numPlotsPerRow 2