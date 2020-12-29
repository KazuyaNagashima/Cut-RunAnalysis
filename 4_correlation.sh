#!/bin/env zsh -e

cd ..
dir=(`pwd`)

echo 'SRAの名前を入力してください(rep1)(.sraの前)'
read accession1

echo 'SRAの名前を入力してください(rep2)(.sraの前)'
read accession2

echo 'リファレンスゲノムを入力して下さい[mm9,mm10]'
read ref

echo '日付(ex:201223)'
read date

data1=(`grep ${accession1} data.csv`)

arr1=( `echo ${data1} | tr -s ',' ' '`)
sra1=(${arr1[1]})
sample1=(${arr1[2]})
seq1=(${arr1[3]})
ends1=(${arr1[4]})
reserch1=(${arr1[5]})
spike1=(${arr1[6]})


data2=(`grep ${accession2} data.csv`)

arr2=( `echo ${data2} | tr -s ',' ' '`)
sra2=(${arr2[1]})
sample2=(${arr2[2]})
seq2=(${arr2[3]})
ends2=(${arr2[4]})
reserch2=(${arr2[5]})
spike2=(${arr2[6]})


mkdir ${dir}/comparison/${date}_correlation
cd ${dir}/comparison/${date}_correlation

touch analysis_detail.txt

echo '<<rep1>>' | tee -a analysis_detail.txt
echo 'SRA:'${sra1} | tee -a analysis_detail.txt
echo 'SAMPLE:'${sample1} | tee -a analysis_detail.txt
echo 'SEQ:'${seq1} | tee -a analysis_detail.txt
echo 'READ:'${ends1} | tee -a analysis_detail.txt
echo 'RESERCH:'${reserch1} | tee -a analysis_detail.txt
echo 'REFERENCE:'${ref} | tee -a analysis_detail.txt
echo 'SPIKE-IN:'${spike1} | tee -a analysis_detail.txt

echo '\n'| tee -a analysis_detail.txt

echo '<<rep2>>' | tee -a analysis_detail.txt
echo 'SRA:'${sra2} | tee -a analysis_detail.txt
echo 'SAMPLE:'${sample2} | tee -a analysis_detail.txt
echo 'SEQ:'${seq2} | tee -a analysis_detail.txt
echo 'READ:'${ends2} | tee -a analysis_detail.txt
echo 'RESERCH:'${reserch2} | tee -a analysis_detail.txt
echo 'REFERENCE:'${ref} | tee -a analysis_detail.txt
echo 'SPIKE-IN:'${spike2} | tee -a analysis_detail.txt

echo '\n'| tee -a analysis_detail.txt
echo 'binsize:1kb' | tee -a analysis_detail.txt


multiBigwigSummary bins --bwfiles ${dir}/${reserch1}/${sample1}_${sra1}/${ref}/deeptools/individual/individual_${sra1}.bw ${dir}/${reserch2}/${sample2}_${sra2}/${ref}/deeptools/individual/individual_${sra2}.bw -p 7 \
-o results.npz --outRawCounts readCounts.tab \
--binSize 1000 --labels rep1 rep2

#可視化
plotCorrelation -in results.npz -c pearson -p scatterplot -o ${date}_scatter \
--removeOutliers #外れ値を除く
