#####set data#########
filename=("Eed_CTR_morula_B6129xPWK_H3K4me3_rep1")
######################

cd ..

dir="/Users/shigenseigyo/Desktop/analysis/cut_and_run"
cd ${dir}
mkdir ${dir}/${filename}
mv ${filename}_R1.fastq ${filename}
mv ${filename}_R2.fastq ${filename}
cd ${dir}/${filename}


#fastpによるQuality Check
fastp -i ${filename}_R1.fastq -I ${filename}_R2.fastq -o fastp_${filename}R1.fastq -O fastp_${filename}R2.fastq -h report_fastp.html -j report_fastp.json -q 20 --length_required 20

mkdir fastp_report
mv report_fastp.html report_fastp.json fastp_report


#bowtie2によるマッピング
bowtie2 -p 7  -t -q -I 150 -X 800 --no-unal  --no-mixed --no-discordant -x /Users/shigenseigyo/Desktop/reference/snpsplit/PWK_PhJ/bowtie2_index/PWK_PhJ.SNP.nmasked -1 fastp_${filename}R1.fastq -2 fastp_${filename}R2.fastq -S ${filename}.sam
#リファレンスゲノムにおけるSNPをNでマスクする必要がある(SNPsplit)
# -t 時間 -q FASTAQ -I 有効ペアエンド最小フラグメント -X 有効ペアエンド最大フラグメント --no-unal 失敗した読み取りのSAMを抑制 --no-mixed 常にペアでアライメント --no-discordant 一致しないとき、-Iや-Xを満たしていない不一致アライメントを探すのを抑制


#multi mapped reads はbowtie2では"XS"とタグ付されているのでgrepで除く
grep -v "XS" ${filename}.sam >${filename}.uq.sam


#samtoolsでbowtie2の出力であるsamをbam(binary)に変換 
samtools view -bS ${filename}.uq.sam > ${filename}.bam


#samtoolsでbam(binary)のソート
samtools sort ${filename}.bam > ${filename}_sorted.bam


#picard toolsでPCRバイアスを除く
picard MarkDuplicates INPUT=${filename}_sorted.bam OUTPUT=${filename}.drm.bam METRICS_FILE=${filename}.metrics AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT


#samtoolsでbam(binary)のインデックスを作成
samtools index ${filename}.drm.bam


#インサートサイズの可視化
picard CollectInsertSizeMetrics INPUT=${filename}.drm.bam OUTPUT=${filename}_insert_size_metrics.txt H=${filename}_insert_size_metrics.pdf MINIMUM_PCT=0


#snpsplitでアレル特異的なリードを抽出
SNPsplit --paired --snp_file /Users/shigenseigyo/Desktop/reference/snpsplit/PWK_PhJ/PWK_PhJ.SNP.txt ${filename}.drm.bam

samtools sort ${filename}.drm.genome1.bam > ${filename}.drm.paternal_sorted.bam
samtools index ${filename}.drm.paternal_sorted.bam

#deeptoolsで正規化と可視化用bigwig(またはbedgraph)ファイルの作成
bamCoverage -p 6 -b ${filename}.drm.paternal_sorted.bam -o ${filename}.paternal.bw -of bigwig --binSize 20 -e 250 --scaleFactor 1
bamCoverage -p 6 -b ${filename}.drm.bam -o ${filename}.bw -of bigwig --binSize 20 -e 250 --scaleFactor 1
# --binSize bigwig / bedgraphファイルの出力用のビンのサイズ（ベース）
# -e 読み取りをフラグメントサイズに拡張できます。設定されている場合、例外なく、各読み取りが延長されます。 注：この機能は、RNA-seqなどのスプライスドリードデータではスキップされた領域を超えてリードを拡張するため、通常はお勧めしません。
# --scaleFactor 計算されたスケーリング係数（該当しない場合は1）にこれが乗算されます 
# bigWigフォーマットでは決められたbin幅でのリード数の補正値が全ゲノムで計算されます. つまり, ゲノム上の各binでのChIP-seqのsignal強度が計算. 
# 各ビンのＲＰＫＭ値は、式 "read counts/((bin_length/1000) × (total_reads/106)) "に従って計算。


#macs2によるピークコール
#macs2 callpeak -t ${filename}.drm.bam -n ${filename} -f BAMPE -g mm -q 0.01 --broad --broad-cutoff 0.1 --nolambda --SPMR --nomodel -B
#macs2 callpeak -t ${filename}.drm.paternal_sorted.bam -n ${filename} -f BAMPE -g mm -q 0.01 --broad --broad-cutoff 0.1 --nolambda --SPMR --nomodel -B

# -f入力ファイルのフォーマット"BAMPE"はbamのペアエンド -g ゲノムサイズ -q q値（最小FDR）カットオフ --broad 広いピークコール --broad-cutoff 広いピークの除去。p値が設定されていない場合、q値で除去
# --nolambda バックグラウンドラムダをローカルラムダとして使用。ピーク候補領域でのローカルバイアスを考慮しない。--SPMR bedgraph(bdg)の作成 --nomodel 読み取りを5 '-> 3'方向に拡張 -B ラムダをbedGraphファイルに制御


#cat ${filename}_peaks.gappedPeak |perl -e 'while(<>){chomp; @array=split/\t/; print "chr$array[0]\t$array[1]\t$array[2]\t$array[4]\t$array[3]\n"}'>${filename}.mspeaks.bed 
#bed2pos.pl ${filename}.mspeaks.bed >${filename}.mspeaks.hb
#annotatePeaks.pl ${filename}.mspeaks.hb mm9 >${filename}_peaks.annotated.txt


#mkdir trush
#mv fastp_${filename}R1.fastq fastp_${filename}R2.fastq ${filename}.sam ${filename}.uq.sam ${filename}.bam ${filename}_sorted.bam ${filename}.mspeaks.bed ${filename}.mspeaks.hb
#mkdir 

