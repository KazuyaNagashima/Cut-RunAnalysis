#!/bin/env zsh -e

cd ..
dir=(`pwd`)

echo 'SRAの名前を入力してください(.sraの前)'
read accession

echo 'リファレンスゲノムを入力して下さい[mm9,mm10]'
read ref


data=(`grep ${accession} data.csv`)

arr=( `echo ${data} | tr -s ',' ' '`)
sra=(${arr[1]})
sample=(${arr[2]})
seq=(${arr[3]})
ends=(${arr[4]})
reserch=(${arr[5]})
spike=(${arr[6]})

#mkdir ${dir}/${reserch}/${sample}_${sra}/${ref}/homer
cd ${dir}/${reserch}/${sample}_${sra}/${ref}/homer

#bed2pos.pl ${dir}/${reserch}/${sample}_${sra}/${ref}/macs2/${sra}_peaks.gappedPeak >${sra}.mspeaks.hb
#annotatePeaks.pl ${sra}.mspeaks.hb ${ref} -annStats annotate.log >${sra}_peaks.annotated.txt

python3.7 ${dir}/program/annotateTofig.py ${sra}
