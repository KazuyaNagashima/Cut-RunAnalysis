#!/bin/env zsh -e

cd ..

#解析の確認
echo 'data.csvに記入したSRAを全てCut＆Run解析します。？[y/n]'
read qestion
if [ $qestion = 'y' ]; then
  echo 'プログラムを実行します'
else
  echo 'プログラムを終了します'
  exit
fi

#メールの設定
echo '@edu.k.u-tokyo.ac.jpのアドレスを設定すると、プログラム終了時にメールが届きます。メールを設定しますか？[y/n]'
read mail
if [ $mail = 'y' ]; then
  echo '@edu.k.u-tokyo.ac.jpの前を入力してください'
  read address
else
  echo 'メールを設定しません'
fi


dir=(`pwd`)


for ref in 'mm9' 'mm10'
do

 cat data.csv | sed -e '1d' | while read line
 do
  list=( `echo ${line} | tr -s ',' ' '`)

  ls ${dir}/${list[5]}/${list[2]}_${list[1]}/${ref}
     if [ $? = 0 ] ; then
        echo '過去に解析したデータがあるため、スキップします。'
     else
        echo '過去に解析したデータがありません。解析します。'
       
       touch ${list[1]}_analysis.detail

       data=(`grep ${list[1]} data.csv`)

       arr=( `echo ${data} | tr -s ',' ' '`)
       sra=(${arr[1]})
       sample=(${arr[2]})
       seq=(${arr[3]})
       ends=(${arr[4]})
       reserch=(${arr[5]})
       spike=(${arr[6]})
      
       echo '<<sample detail>>' | tee -a ${list[1]}_analysis.detail
       echo 'SRA:'${sra} | tee -a ${list[1]}_analysis.detail
       echo 'SAMPLE:'${sample} | tee -a ${list[1]}_analysis.detail
       echo 'SEQ:'${seq} | tee -a ${list[1]}_analysis.detail
       echo 'READ:'${ends} | tee -a ${list[1]}_analysis.detail
       echo 'RESERCH:'${reserch} | tee -a ${list[1]}_analysis.detail
       echo 'REFERENCE:'${ref} | tee -a ${list[1]}_analysis.detail
       echo 'SPIKE-IN:'${spike} | tee -a ${list[1]}_analysis.detail

       mkdir ${dir}/${reserch}
       mkdir ${dir}/${reserch}/${sample}_${sra}
       mkdir ${dir}/${reserch}/${sample}_${sra}/${ref}
       cd ${dir}/${reserch}/${sample}_${sra}/${ref}
       
       #pfastq-dumpでSRAをfastqに変換
       if [ ${ends} = 'single' ]; then
         pfastq-dump --threads 8  --gzip -s ${dir}/sra/${sra}.sra
       elif [ ${ends} = 'pair' ]; then
         #Pair Endsの時は--split-filesをつける
         pfastq-dump --threads 8 --gzip --split-files -s ${dir}/sra/${sra}.sra
       else
         echo 'readをsingleかpairに設定してください'
         echo 'プログラムを終了します'
         exit
       fi

       mkdir fastq
       if [ ${ends} = 'single' ]; then
        mv ${sra}_.fastq.gz fastq
       else
        mv ${sra}_1.fastq.gz ${sra}_2.fastq.gz fastq
       fi

       #fastpによるQuality Check(20bp>,q>20)
        mkdir fastp_reports
       if [ ${ends} = 'single' ]; then
        fastp -i ./fastq/${sra}_.fastq.gz -o ./fastp_reports/fastp_${sra}.fastq -h ./fastp_reports/report_fastp.html -j ./fastp_reports/report_fastp.json -q 20 --length_required 20
       else
        fastp -i ./fastq/${sra}_1.fastq.gz -I ./fastq/${sra}_2.fastq.gz -o ./fastp_reports/fastp_${sra}R1.fastq -O ./fastp_reports/fastp_${sra}R2.fastq -h ./fastp_reports/report_fastp.html -j ./fastp_reports/report_fastp.json -q 20 --length_required 20
       fi

       #bowtie2によるマッピング
        mkdir bowtie2
        if [ ${ends} = 'single' ]; then
         bowtie2 -p 8  -t -q -N 1 --no-unal -x /Users/shigenseigyo/Desktop/reference/Mus_musculus/UCSC/${ref}/Sequence/Bowtie2Index/genome -U ./fastp_reports/fastp_${sra}.fastq -S ./bowtie2/${sra}.sam
         rm ./fastp_reports/fastp_${sra}.fastq
        else
         bowtie2 -p 8  -t -q -N 1 -I 150 -X 800 --no-unal  --no-mixed --no-discordant -x /Users/shigenseigyo/Desktop/reference/Mus_musculus/UCSC/${ref}/Sequence/Bowtie2Index/genome -1 ./fastp_reports/fastp_${sra}R1.fastq -2 ./fastp_reports/fastp_${sra}R2.fastq -S ./bowtie2/${sra}.sam
         rm ./fastp_reports/fastp_${sra}R1.fastq ./fastp_reports/fastp_${sra}R2.fastq
        fi
         # -t 時間 -q FASTAQ -N ミスマッチ -I 有効ペアエンド最小フラグメント -X 有効ペアエンド最大フラグメント --no-unal 失敗した読み取りのSAMを抑制 --no-mixed 常にペアでアライメント --no-discordant 一致しないとき、-Iや-Xを満たしていない不一致アライメントを探すのを抑制


         #multi mapped reads はbowtie2では"XS"とタグ付されているのでgrepで除く
         grep -v "XS" ./bowtie2/${sra}.sam >./bowtie2/${sra}.uq.sam
         rm ./bowtie2/${sra}.sam

         #samtoolsでbowtie2の出力であるsamをbam(binary)に変換
         samtools view -bS ./bowtie2/${sra}.uq.sam > ./bowtie2/${sra}.bam
         rm ./bowtie2/${sra}.uq.sam

         #samtoolsでbam(binary)のソート
         samtools sort ./bowtie2/${sra}.bam > ./bowtie2/${sra}_sorted.bam


         #picard toolsでPCRバイアスを除く
         picard MarkDuplicates INPUT=./bowtie2/${sra}_sorted.bam OUTPUT=./bowtie2/${sra}.drm.bam METRICS_FILE=./bowtie2/${sra}.metrics AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT


         #samtoolsでbam(binary)のインデックスを作成
         samtools index ./bowtie2/${sra}.drm.bam


         #picard toolsでinsert size を確認
         mkdir picard
         picard CollectInsertSizeMetrics INPUT=./bowtie2/${sra}.drm.bam OUTPUT=./picard/${sra}_insert_size_metrics.txt H=./picard/${sra}_insert_size_metrics.pdf MINIMUM_PCT=0



         if [ $spike = 'none' ] ; then
            #deeptoolsで正規化と可視化用bigwig(またはbedgraph)ファイルの作成
            mkdir deeptools
            bamCoverage -p 6 -b ./bowtie2/${sra}.drm.bam -o ./deeptools/${sra}.bw -of bigwig --binSize 20 -e 250 --scaleFactor 1
            # --binSize bigwig / bedgraphファイルの出力用のビンのサイズ（ベース）
            # -e 読み取りをフラグメントサイズに拡張できます。設定されている場合、例外なく、各読み取りが延長されます。 注：この機能は、RNA-seqなどのスプライスドリードデータではスキップされた領域を超えてリードを拡張するため、通常はお勧めしません。
            # --scaleFactor 計算されたスケーリング係数（該当しない場合は1）にこれが乗算されます 
            # bigWigフォーマットでは決められたbin幅でのリード数の補正値が全ゲノムで計算されます. つまり, ゲノム上の各binでのChIP-seqのsignal強度が計算. 
            # 各ビンのＲＰＫＭ値は、式 "read counts/((bin_length/1000) × (total_reads/106)) "に従って計算。


            #macs2によるピークコール
            mkdir macs2
            if [ ${ends} = 'single' ]; then
               macs2 callpeak -t ./bowtie2/${sra}.drm.bam -n ${sra} -g mm --outdir macs2 -q 0.01 --broad --broad-cutoff 0.1 --nolambda --SPMR --nomodel -B 
               # -f入力ファイルのフォーマット"BAMPE"はbamのペアエンド -g ゲノムサイズ -q q値（最小FDR）カットオフ --broad 広いピークコール --broad-cutoff 広いピークの除去。p値が設定されていない場合、q値で除去
               # --nolambda バックグラウンドラムダをローカルラムダとして使用。ピーク候補領域でのローカルバイアスを考慮しない。--SPMR bedgraph(bdg)の作成 --nomodel 読み取りを5 '-> 3'方向に拡張 -B ラムダをbedGraphファイルに制御
            else
               macs2 callpeak -t ./bowtie2/${sra}.drm.bam -n ${sra} -f BAMPE -g mm --outdir macs2 -q 0.01 --broad --broad-cutoff 0.1 --nolambda --SPMR --nomodel -B 
            fi
         else

            mkdir spike-in

            if [ ${ends} = 'single' ]; then
               bowtie2 -p 8  -t -q -N 1 --no-unal -x /Users/shigenseigyo/Desktop/reference/spike-in/UCSC/${spike}/Sequence/Bowtie2Index/genome -U ./fastp_reports/fastp_${sra}.fastq -S ./spike-in/${spike}_${sra}.sam
            else
               bowtie2 -p 8  -t -q -N 1 -I 150 -X 800 --no-unal  --no-mixed --no-discordant -x /Users/shigenseigyo/Desktop/reference/spike-in/UCSC/${spike}/Sequence/Bowtie2Index/genome -1 ./fastp_reports/fastp_${sra}R1.fastq -2 ./fastp_reports/fastp_${sra}R2.fastq -S ./spike-in/${spike}_${sra}.sam
            fi

            samtools view -bS ./spike-in/${spike}_${sra}.sam > ./spike-in/${spike}_${sra}.bam
            rm ./spike-in/${spike}_${sra}.sam

            samtools sort ./spike-in/${spike}_${sra}.bam > ./spike-in/${spike}_${sra}_sorted.bam

            samtools index ./spike-in/${spike}_${sra}_sorted.bam 

         fi

     fi
 done

done


#終了時メールで送信
if [ $mail = 'y' ]; then
  echo "Process done　メールを送信しました。" | mail -s "Process done" ${address}@edu.k.u-tokyo.ac.jp
else
  echo "Process done"
fi