# bulkATACseq

# pre-analysis

## 1,trim reads/adaptor --*trim_galore* 

```shell
trim_galore --quality 20 --paired ${meta}${i}.1.clean.fq.gz ${meta}${i}.2.clean.fq.gz -o ${trim}
```

## 2,mapping reads to genome -- *bwa mem*

```shell
bwa mem -t 4 -o ${sam}${i}.sam ${Mouse_index} ${trim}${i}.1*fq.gz ${trim}${i}.2*fq.gz
```

## 3,remove duplicates --*picard* 

``` shell
##convert sam file to bam file
samtools view -q 20 -hbS -o ${prebam}${i}.bam ${sam}${i}.sam 
##sort bam
samtools sort -o ${prebam}${i}.sorted.bam ${prebam}${i}.bam
##remove duplicates
java -jar $picard MarkDuplicates M=dupstats REMOVE_DUPLICATES=TRUE I=${prebam}${i}.sorted.bam o=${prebam}${i}.sorted.NODUPS.bam
```

## 4,remove  mitochondria and chrY --*samtools*

 ``` shell
samtools view -h ${prebam}${i}.sorted.NODUPS.bam | awk '$3 != "chrM" {print $0}'| awk '$3 != "chrY" {print $0}'|samtools view -Sb - > ${prebam}${i}.sorted_NODUPS_NOMT.bam
 ```

## 5,get proper paired mapping bam --*samtools*

​	bam mapping situation can be checked by *samtools flagstate*

``` shell
samtools flagstat [bam_file] 
```



​	*samtools flagstate* 输出结果为：![clip_image001](C:\Users\skye\AppData\Local\Temp\msohtmlclip1\01\clip_image001.png)

![](C:\Users\skye\Documents\md\picture\clip_image001.png)

```shell
##first fixmate bam flag
samtools sort -n ${prebam}${i}.sorted_NODUPS_NOMT.bam -o -|samtools fixmate - ${trimb}${i}.fixmate.bam
## get proper paired bam
samtools view -f 2 -bF 12 ${trimb}${i}.fixmate.bam -o ${pair}${i}.pair.bam
```

​	note :	-f 筛选与参数mapping的值（正向选择），2 代表reads正确mapping

​			-b 输出bam格式，默认为sam格式

​			-F 筛选与参数不mapping的值（反向选择），12：拆分为4：该序列没有比对到参考序列上，8：该序列的mate没有比对到参考序列上

## 6, fix adaptor shift in Tn5 enzyme digestion --*deeptools*

```shell
##alignmentSieve need an indexed bam file
##if samtools sort doesn't work ,bedtools sort may help
samtools sort ${pair}${i}.pair.bam -o ${pair}${i}.sort.pair.bam
samtools index ${pair}${i}.sort.pair.bam
alignmentSieve -b ${pair}${i}.sort.pair.bam -o ${shifted}${i}.shifted.bam --ATACshift
samtools sort ${shifted}${i}.shifted.bam -o ${shifted}${i}.shifted.sort.bam
samtools index ${shifted}${i}.shifted.sort.bam
```

after these steps, we finally get  an indexed proper bam file

## 7,call peak by shifted bam file --*Macs*

```shell
macs2 callpeak -t ${shifted}${i}.shifted.bam -n ${peak}${i}_MACS -g mm -f BAM --nomodel --shift -100 --extsize 200
###when generate peak in all replicate,we add all replicate's shifted bam after '-t'  
```



# Quality Control

## 1, replicate correlation

### a,first get reads count by bins --*deeptools*

```shell
multiBamSummary bins --bamfiles [path]/PD10_BF_12_5FL1-5.shifted.sort.bam [path]/PD10_BF_12_5FL2-5.shifted.sort.bam [path]/PD10_BF_12_5FL4-5.shifted.sort.bam --outRawCounts ../correlation/FL12_cor.txt -o ../correlation/FL12.5_results.npz
```

### b,plot correlation map --*scatterplot3d*（10，9）

 ``` R
library(scatterplot3d)
smp<-read.table("data1/result/correlation/txt/NB_cor.txt")
smp<-smp[,4:6]
a<-which(apply(smp,1,sum) == 0)
smp<-smp[-a,]
smp_d<-log2(smp)
scatterplot3d(smp_d,main="NB",color = "blue",pch = 16,angle = 70,type = "p",box = F,xlab = "NB_1",ylab = "NB_2",zlab = "NB_3")
##get correlation score
dat<-matrix(nrow = 1,ncol = 3)
colnames(dat)<-c("1-2","2-3","1-3")
dat[1,1]<-cor(smp$V4,smp$V5,method = "pearson")
dat[1,2]<-cor(smp$V6,smp$V5,method = "pearson")
dat[1,3]<-cor(smp$V4,smp$V6,method = "pearson")
 ```

the output is like:

​	![1560436258812](C:\Users\skye\Documents\md\picture\1560436258812.png)

## 2, fragsize_distribution（6，5）--*ATACseqQC*

``` R
library(ATACseqQC)
fragSize<-fragSizeDist(bamfile,bamFile.labels="bamname")
```

the output is :

​	![1560436822904](C:\Users\skye\Documents\md\picture\1560436822904.png)

## 3, TSS-distribution

### a, peak frequency profile around TSS --*ChIPseeker*

``` R
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
promoter<-getPromoters(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,upstream = 2000,downstream = 2000)
files<-list("data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/PL12_combine_MACS_peaks.narrowPeak","data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/FL12_combine_MACS_peaks.narrowPeak","data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/FL16_combine_MACS_peaks.narrowPeak")
names(files)<-c("PL12","FL12","FL16")
tagmatrixlist<-lapply(files,getTagMatrix,windows=promoter)
plotAvgProf(tagmatrixlist,xlim = c(-2000,2000))

```

the output is :

![1560437274207](C:\Users\skye\Documents\md\picture\1560437274207.png)

note : the output is not like results from deeptools plotprofile ,plotprofile using reads count as ordinates while ChIPseeker using peak count frequency

### b, percent of peaks annotate to promoter vs random peaks distribution --*bedtools*

generate random peaks from Macs peaks

``` shell
out=/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/peak_bed/random_peak/
genome=/home/skye/data1/atac_seq_2018/database/meme/mm10_nochrMY.chrom.sizes
for i in $(ls *narrowPeak |cut -d "." -f 1)
do
	bedtools shuffle -i ${i}.narrowPeak -g $genome -noOverlapping >${out}${i}_random.bed
done
```

get promoter percent around all annotation

``` R
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
files<-list.files(path = "/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/bed",pattern = "bed")
path<-"/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/bed/"
out<-"/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/bed/anno/"
for(i in 1:length(files)){
    peak<-paste(path,files[i],sep = "")
    a<-readPeakFile(peak)
  ##annotate peak
    a<-annotatePeak(a,TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,annoDb = "org.Mm.eg.db")
    a<-as.data.frame(a)
    ##grep promoter peak
    b<-grep(pattern = "Promoter",x = a$annotation)
    ##caculate promoter distribution
    percent<-length(b)/length(a[,1])
  	percent<-as.data.frame(percent)
  	write.table(a,paste(out,files[i],"_anno.txt",sep = ""),sep = "\t",quote = F,col.names = T,row.names = F)
  	write.table(percent,paste(out,files[i],"_percent.txt",sep = ""),sep = "\t",quote = F,col.names = T,row.names = F)
}
```

then using excel to do the bar plot.

# peak reads distribution

## 1, fragment length distribution

```shell
#!/bin/bash
##extract fragsize
bam=/home/skye/data1/testout/cbam/bamlist  
profile=/home/skye/data1/testout/cbam/
frag=/home/skye/data1/testout/fragment2/
par=/home/skye/data1/testout/par/
names=$(cut -d "/" -f 7 $bam | cut -d "." -f 1)
for i in ${names} 
do 
samtools view -bF 12 ${profile}${i}.sorted_NODUPS_NOMT.bam | tee ${par}${i}.paired.bam | awk '{print $9}' > ${frag}${i}_frag.txt
done
fragfile=/home/skye/data1/testout/fragment2/
frag2=/home/skye/data1/testout/frag2/
ls /home/skye/data1/testout/fragment2/ > /home/skye/data1/testout/frag.list
frag=/home/skye/data1/testout/frag.list
names=$(cut -d '/' -f 7 $frag | cut -d '.' -f 1)
for i in $names
do
cat ${fragfile}${i}.txt | awk '{if($1>0 && $1 <= 100) print $1}' | tee ${frag2}${i}_free.txt |wc -l >> ${frag2}${i}_12345.txt 
cat ${fragfile}${i}.txt | awk '{if($1>180 && $1 <= 247) print $1}' | tee ${frag2}${i}_mono.txt | wc -l >> ${frag2}${i}_12345.txt 
cat ${fragfile}${i}.txt | awk '{if($1>315 && $1 <= 473) print $1}' | tee ${frag2}${i}_di.txt |wc -l >> ${frag2}${i}_12345.txt
cat ${fragfile}${i}.txt | awk '{if($1>558 && $1 <= 615) print $1}' | tee ${frag2}${i}_tri.txt |wc -l >>${frag2}${i}_12345.txt
cat ${fragfile}${i}.txt | awk '{if($1>615) print $1}' | tee ${frag2}${i}_distal.txt | wc -l >>${frag2}${i}_12345.txt
done

```

the fragment length distribution is separate to ,

​	nucleosome free: 0<fragment length  <=100

​	mononucleosome :180<fragment length  <=247

​	dinucleosome:315<fragment length  <=473

​	tri-nucleosome:558<fragment length  <=615

​	more than tri-nuclleosome: >615



## ***2, reads count profile --*ChIPseeker*

### a, profile around TSS 

```R
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
promoter<-getPromoters(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,upstream = 2000,downstream = 2000)
##files is a list of all narrowPeak paths
files<-list("data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/PL12_combine_MACS_peaks.narrowPeak","data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/FL12_combine_MACS_peaks.narrowPeak","data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/FL16_combine_MACS_peaks.narrowPeak")
names(files)<-c("PL12","FL12","FL16")
tagmatrixlist<-lapply(files,getTagMatrix,windows=promoter)
plotAvgProf(tagmatrixlist,xlim = c(-2000,2000))
```

the output is like:

![1560584274644](C:\Users\skye\Documents\md\picture\1560584274644.png)

### b, profile around peak summit(atlas peak summit)

```R
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
atlas_summit_bed<-read.table("~/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/atlas_file/all_sample_atlas/all_atlas")
atlas_summit_bed[,4]<-round((atlas_summit_bed[,2]+atlas_summit_bed[,3])/2)-1000
atlas_summit_bed[,5]<-round((atlas_summit_bed[,2]+atlas_summit_bed[,3])/2)+1000
atlas_summit_bed<-atlas_summit_bed[,c(1,4,5)]
colnames(atlas_summit_bed)<-c("chr","start","end")
atlas_summit_1kb<-makeGRangesFromDataFrame(atlas_summit_bed)
data<-list("/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/PL12_combine_MACS_peaks.narrowPeak","/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/FL12_combine_MACS_peaks.narrowPeak","/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/Macs_peak/FL16_combine_MACS_peaks.narrowPeak")
names(data)<-c("PL12","FL12","FL16")
tagmatrixlist<-lapply(data, getTagMatrix,windows=atlas_summit_1kb)
plotAvgProf(tagmatrixlist,xlim = c(-1000,1000))
```

the output is:

​				![1560585518342](C:\Users\skye\Documents\md\picture\1560585518342.png)

​	

# atlas of all ATAC-seq peaks (PCA)

we use *diffbind* to make PCA analysis

first we need a sample sheet

<table>
   <tr>
      <td>SampleID</td>
      <td>Tissue</td>
      <td>Replicate</td>
      <td>bamReads</td>
      <td>Peaks</td>
      <td>PeakCaller</td>
      <td>ScoreCol</td>
      <td>LowerBetter</td>
   </tr>
   <tr>
      <td>PL12_1</td>
      <td>PL12</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_12_5PL1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_12_5PL1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>PL12_2</td>
      <td>PL12</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_12_5PL2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_12_5PL2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>PL12_3</td>
      <td>PL12</td>
      <td>3</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_12_5PL3-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_12_5PL3-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>FL12_1</td>
      <td>FL12</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_12_5FL1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_12_5FL1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>FL12_2</td>
      <td>FL12</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_12_5FL2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_12_5FL2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>FL12_3</td>
      <td>FL12</td>
      <td>3</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_12_5FL4-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_12_5FL4-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>FL16_1</td>
      <td>FL16</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_16_5FL1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_16_5FL1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>FL16_2</td>
      <td>FL16</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_16_5FL2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_16_5FL2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>FL16_3</td>
      <td>FL16</td>
      <td>3</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/shiftedbam_pff/PD10_BF_16_5FL3-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/PL_FL/peak/PD10_BF_16_5FL3-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>NB_1</td>
      <td>NB</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_NB_BM1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_NB_BM1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>NB_2</td>
      <td>NB</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_NB_BM2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_NB_BM2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>NB_3</td>
      <td>NB</td>
      <td>3</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_NB_BM3-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_NB_BM3-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W1_1</td>
      <td>W1</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_1W_BM1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_1W_BM1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W1_2</td>
      <td>W1</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_1W_BM2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_1W_BM2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W1_3</td>
      <td>W1</td>
      <td>3</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_1W_BM3-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_1W_BM3-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W3_1</td>
      <td>W3</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_3W_BM1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_3W_BM1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W3_2</td>
      <td>W3</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_3W_BM2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_3W_BM2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W4_1</td>
      <td>W4</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_4W_BM1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_4W_BM1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W4_2</td>
      <td>W4</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_4W_BM2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_4W_BM2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W4_3</td>
      <td>W4</td>
      <td>3</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_4W_BM3-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_4W_BM3-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W6_1</td>
      <td>W6</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_6W_BM1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_6W_BM1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W6_2</td>
      <td>W6</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_6W_BM2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_6W_BM2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W6_3</td>
      <td>W6</td>
      <td>3</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_6W_BM3-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_6W_BM3-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W8_1</td>
      <td>W8</td>
      <td>1</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_8W_BM1-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_8W_BM1-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W8_2</td>
      <td>W8</td>
      <td>2</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_8W_BM2-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_8W_BM2-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td>W8_3</td>
      <td>W8</td>
      <td>3</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/shiftedbam/PD10_BF_8W_BM3-5.shifted.sort.bam</td>
      <td>/home/skye/data1/atac_seq_2018/pre_data/BM/peak/PD10_BF_8W_BM3-5_MACS_peaks.narrowPeak</td>
      <td>narrow</td>
      <td>8</td>
      <td>FALSE</td>
   </tr>
   <tr>
      <td></td>
   </tr>
</table>


then we use diffbind to make PCA

reference : https://www.jianshu.com/p/f849bd55ac27 sample<-read.csv("data1/result/diff/diffbind_sample.csv")

a<-dba(sampleSheet = sample)

b<-dba.count(a,bUseSummarizeOverlaps = TRUE)

dba.plotPCA(b,attributes = DBA_TISSUE,label = DBA_ID)

## delete 'label= DBA_ID' to remove labels in dot

the output is like:

​	![1560434411947](C:\Users\skye\Documents\md\picture\1560434411947.png)



# find differential peaks



## 1, make atlas

### a, +-250bp peak summit

​	call peak by all replicates

- get peak summit +- 250bp bed file

  ``` R
  a<-read.table([macs_summit.bed])
  a[2]<-a[2]-250
  a[3]<-a[3]+250
  write.table(a[,1:3],[summit_250bp.bed],sep = "\t",col.names = F,row.names = F,quote = F)
  ```

  peak summit file is from macs2 call peak in pre-analysis

### b, merge peaks --*bedtools*

``` shell
## cat all sample peaks together
cat [sample1_peak_bed] [sample2_peak_bed] ... >all_sample_peak.bed
## sort peak file
##if summit -250 bp <0 there is an error, set negative value to 0
sortBed -i all_sample_peak.bed >all_sample_peak_sorted.bed
## merge peaks when overlap >0
bedtools merge -i all_sample_peak_sorted.bed >all_sample_peak_atlas
```

### c, get sample reads count --*deeptools*

```shell
multiBamSummary BED-file --BED ${atlas_file} --bamfiles ${pff_bamfile}*sort.bam ${BM_bamfile}*sort.bam -o ${out}all_count.npz --outRawCounts ${out}all_count.txt
##convert 'chr start end' to 'chr_start_end' as a unique peak name 
cat ${out}all_count.txt |awk '{print $1 "_" $2 "_" $3 }' >a
paste a ${out}all_count.txt >b
## exchange first 3 column to chr_start_end
#cat b |awk '{print $1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 }' >${out}all_countnames.txt
##update
cut -f2,3,4 --complement b >${out}all_countnames.txt
```

## 2, differential analysis

### a, *deseq2*

​	deseq2 needs raw counts to do the job ,do not use the normalized data

```R
###all_atlas peak_diff
mycounts<-read.delim("data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/all_countnames.txt") ###counts from multibamsummary
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
colnames(mycounts)<-c("FL12_1","FL12_2","FL12_3","PL12_1","PL12_2","PL12_3","FL16_1","FL16_2","FL16_3","W1_1","W1_2","W1_3","W3_1",
                      "W3_2","W4_1","W4_2","W4_3","W6_1","W6_2","W6_3","W8_1","W8_2","W8_3","NB_1","NB_2","NB_3")
coldata<-data.frame(row.names = colnames(mycounts),condition=c(rep("FL12",3),rep("PL12",3),rep("FL16",3),rep("W1",3),rep("W3",2),rep("W4",3),rep("W6",3),rep("W8",3),rep("NB",3)))
all(rownames(coldata) == colnames(mycounts))
library(DESeq2)
dds<-DESeqDataSetFromMatrix(mycounts,coldata,design = ~condition)
dds <- DESeq(dds)

res<-results(dds,contrast = c("condition","FL12","PL12"))
FL12_PL12<-as.data.frame(res)
FL12_PL12_up<-FL12_PL12[which(FL12_PL12$log2FoldChange > 1),]
FL12_PL12_down<-FL12_PL12[which(FL12_PL12$log2FoldChange < -1),]

res<-results(dds,contrast = c("condition","FL16","FL12"))
FL16_FL12<-as.data.frame(res)
FL16_FL12_up<-FL16_FL12[which(FL16_FL12$log2FoldChange > 1),]
FL16_FL12_down<-FL16_FL12[which(FL16_FL12$log2FoldChange < -1),]

res<-results(dds,contrast = c("condition","NB","FL16"))
NB_FL16<-as.data.frame(res)
NB_FL16_up<-NB_FL16[which(NB_FL16$log2FoldChange > 1),]
NB_FL16_down<-NB_FL16[which(NB_FL16$log2FoldChange < -1),]

res_bm<-results(dds,contrast = c("condition","W1","NB"))
W1_NB<-as.data.frame(res_bm)
W1_NB_UP<-W1_NB[which(W1_NB$log2FoldChange > 1),]
W1_NB_DOWN<-W1_NB[which(W1_NB$log2FoldChange < -1),]

res_bm<-results(dds,contrast = c("condition","W3","W1"))
W3_NB1<-as.data.frame(res_bm)
W3_NB1_UP<-W3_NB1[which(W3_NB1$log2FoldChange > 1),]
W3_NB1_DOWN<-W3_NB1[which(W3_NB1$log2FoldChange < -1),]

res_bm<-results(dds,contrast = c("condition","W4","W3"))
W4_W3<-as.data.frame(res_bm)
W4_W3_UP<-W4_W3[which(W4_W3$log2FoldChange > 1),]
W4_W3_DOWN<-W4_W3[which(W4_W3$log2FoldChange < -1),]

res_bm<-results(dds,contrast = c("condition","W6","W4"))
W6_W4<-as.data.frame(res_bm)
W6_W4_UP<-W6_W4[which(W6_W4$log2FoldChange > 1),]
W6_W4_DOWN<-W6_W4[which(W6_W4$log2FoldChange < -1),]

res_bm<-results(dds,contrast = c("condition","W8","W6"))
W8_W6<-as.data.frame(res_bm)
W8_W6_UP<-W8_W6[which(W8_W6$log2FoldChange > 1),]
W8_W6_DOWN<-W8_W6[which(W8_W6$log2FoldChange < -1),]


write.table(NB_FL16_up,"data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/NB_FL16_up.txt",sep = "\t",quote = F,col.names = T,row.names = T)
write.table(NB_FL16_down,"data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/NB_FL16_down.txt",sep = "\t",quote = F,col.names = T,row.names = T)
write.table(FL16_FL12_down,"data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/FL16_FL12_down.txt",sep = "\t",quote = F,col.names = T,row.names = T)
write.table(FL16_FL12_up,"data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/FL16_FL12_up.txt",sep = "\t",quote = F,col.names = T,row.names = T)
write.table(FL12_PL12_up,"data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/FL12_PL12_up.txt",sep = "\t",quote = F,col.names = T,row.names = T)
write.table(FL12_PL12_down,"data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/FL12_PL12_down.txt",sep = "\t",quote = F,col.names = T,row.names = T)

```





# regulatory of TF

## 1,motif enrichment --*homer*

```shell
##first we need to convert 'chr_start_end' to 'chr start end' bed format
cat diff_peak.txt |awk -F "[_\t]" '{print $1 "\t" $2 "\t" $3}' >diff_peak.bed
## delete colnames in each bed file
sed -i '1d' diff_peak.bed
## homer find motif in each sample
findMotifsGenome.pl diff_peak.bed mm10 [out dir]
```



## 2, find significant motif

first we extract motif and it's significant motif P-value

$1 is motif name ;$3 is motif enrich P-value

the $2 in output file is "-log10P-value"

```shell
cat motif_known.txt |awk '{print $1 "\t" $3}' |awk -F "[\t][1][e][-]" '{print $1 "\t" $2}' >motif_known_pvalue.txt
```



### a, diff_peak motif heatmap --*pheatmap*

``` R
FL16_FL12_UP<-read.table("data1/atac_seq_2018/pre_data/PL_FL/deeptools_qu/multibamsummary/peak_diff/motif/all_knowmotif/motif_log10p/FL12_FL16_down_lgp")
FL16_FL12_DOWN<-read.table("data1/atac_seq_2018/pre_data/PL_FL/deeptools_qu/multibamsummary/peak_diff/motif/all_knowmotif/motif_log10p/FL12_FL16_up_lgp")
colnames(FL16_FL12_DOWN)<-c("name","value")
colnames(FL16_FL12_UP)<-c("name","value")
FL16_FL12<-merge(FL16_FL12_UP,FL16_FL12_DOWN,by="name",all = T)
colnames(FL16_FL12)<-c("name","FL16_FL12_up","FL16_FL12_down")
FL16_FL12<-_FL16_FL12[!duplicated(FL16_FL12$name),]
rownames(FL16_FL12)<-FL16_FL12$name
FL16_FL12<-FL16_FL12[,-1]
library(pheatmap)
pheatmap(FL16_FL12,scale = "column")
```

we extract significant motif in each "column" by manual (refer to color depth)

```R
## significant motif heatmap plot
##select most significant motif name in each column
sig<-read.table("data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/motif/FL16_FL12_sort.txt")
FL16_FL12$name<-rownames(FL16_FL12)
##the format in sig file and FL16_FL12 is different
FL16_FL12$name<-gsub("-","−",FL16_FL12$name)
##select sig motif order in p-value file 
b=1
c=c()
for(i in sig$V1){
  a<-which(FL16_FL12$name == i)
  c[b]=a
  b=b+1
}
##reorder motif rank ,in p-value file ,for adjusting heatmap cluster to make heatmap more beauty
FL16_FL12_select<-rbind(FL16_FL12[c,],FL16_FL12[-c,])
pheatmap(FL16_FL12[c,1:2],cluster_rows = F,cluster_cols = F,show_rownames = F)
```

the output is like

![1560564950012](C:\Users\skye\Documents\md\picture\1560564950012.png)

### b, diff_peak motif dot plot --*ggplot2*（5，3）

```R
fl12_pl12_down<-read.table("data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/motif/all_p/FL12_PL12_down_known_p.txt")
fl12_pl12_up<-read.table("data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/motif/all_p/FL12_PL12_up_known_p.txt")
library(pheatmap)
library(ggplot2)
fl12_pl12_down<-fl12_pl12_down[-1,]
fl12_pl12_down$V2<-as.numeric(fl12_pl12_down$V2)
##first no color scale ,just check top ranking motif's type
ggplot(data = fl12_pl12_down,aes(x=order(fl12_pl12_down$V2,decreasing = T),y=fl12_pl12_down$V2))+geom_point()
##select by manual (refer to the plot above)
##add color level
fl12_pl12_down$color<-factor(c(rep("bZIP",8),rep("Zf",3),"bZIP",rep("Zf",2),rep("other",length(fl12_pl12_down$V1)-14)),levels = c("bZIP","Zf","other"))
###4 color
##ETS:"#f8766d",Zf:"#00b0f6",Runt:"#5684E9",bZIP:"#FF9900"/"#00b0f6"（rescue）,other:"black"
ggplot(data = a,aes(x=order(a$V2,decreasing = T),y=a$V2))+geom_point(position = position_jitter(width = 1,height = 0.5),aes(color=a$color))+scale_color_manual(values = c("#FF9900","#f8766d","#5684E9","black"))+
  theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
###3 color
ggplot(data = a,aes(x=order(a$V2,decreasing = T),y=a$V2))+geom_point(position = position_jitter(width = 2,height = 2),aes(color=a$color))+scale_color_manual(values = c("#FF9900","#00b0f6","black"))+
  theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
```

the output is:

![1560565022157](C:\Users\skye\Documents\md\picture\1560565022157.png)



## 2, significant motif target gene

### a, scan for target site --*homer*

homer has motif pwm database in it's install file path

```shell
fa_file=/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/diff_peak_fa/
motif_database=/home/skye/Downloads/biosoft/homer/motifs/
significant_motif=/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/motif/significant_motif/
target=/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/motif/significant_motif/target/

for motif_sample in $(ls ${significant_motif}*txt)
do
	name=$(echo ${motif_sample} |awk -F "/" '{print $NF}' |cut -d '.' -f1 |cut -d '_' -f2,3,4)
	mkdir ${target}${name}
	for motif in $(cat ${motif_sample} |awk -F "(" '{print $1}' |tr 'A-Z' 'a-z')
	do
		scanMotifGenomeWide.pl ${motif_database}${motif}.motif ${fa_file}${name}.fa -bed > ${target}${name}/${motif}.bed
	done
done 
```

note : some motif name is not proper match motif name in database ,needs manual adjust

### b, annotate target site to genes and make GO/KEGG analysis --*ChIPseeker*

``` shell
## this is a loop script for batch operation
significant_motif=/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/motif/significant_motif/
target=/home/skye/data1/atac_seq_2018/pre_data/PL_FL/replicate_peak/qu_250bp_bed/all_atlas_diffpeak/diff_peak/motif/significant_motif/target/
for motif_sample in $(ls ${significant_motif}*txt) 
do
	name=$(echo ${motif_sample} |awk -F "/" '{print $NF}' |cut -d '.' -f1 |cut -d '_' -f2,3,4)
	for motif in $(cat ${motif_sample} |awk -F "(" '{print $1}' |tr 'A-Z' 'a-z')
	do
		cat ${target}${name}/${motif}.bed |awk -F "[:-\t]" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >${target}${name}/${motif}_peak.bed
		mkdir ${target}${name}/anno_cluster
		mkdir ${target}${name}/anno_cluster/${motif}
		Rscript /home/skye/data1/atac_seq_2018/script/anno_cluste_inside.R ${motif} ${target}${name}/ ${target}${name}/anno_cluster/${motif}/
	done 
done

```

the Rscript in this shell script is the main function to make annotation

```R
##argv[1]:motif name;argv[2]:peak file path ;argv[3]:outfile_path
argv<-commandArgs(TRUE)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
peak<-read.table(paste(argv[2],argv[1],"_peak.bed",sep = ""))
peak$v6<-apply(peak[,c(2,4)],1,sum)
peak$v7<-apply(peak[,c(2,5)],1,sum)
peak<-peak[,c(1,6,7)]
write.table(peak,paste(argv[2],argv[1],"_seq.bed",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
bed<-readPeakFile(peakfile = paste(argv[2],argv[1],"_seq.bed",sep = ""))
bed<-annotatePeak(peak = bed,TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,annoDb = "org.Mm.eg.db" )
pdf(paste(argv[3],argv[1],"anno_bar",".pdf",sep = ""),width = 6,height = 3)
plotAnnoBar(bed)
dev.off()
bed<-as.data.frame(bed)
write.table(bed,paste(argv[3],argv[1],"_annoda.txt",sep = ""),sep = "\t",quote = F)
##GO
GO<-enrichGO(gene = bed$SYMBOL[!duplicated(bed$SYMBOL)],OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont="ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
pdf(paste(argv[3],argv[1],"go_25_bar",".pdf",sep = ""),width = 10,height = 6)
barplot(GO,showCategory = 25)
dev.off()
GO<-as.data.frame(GO)
write.table(GO,paste(argv[3],argv[1],"_anno_GO.txt",sep = ""),sep = "\t",quote = F)
##kegg
kegg<-enrichKEGG(gene = bed$geneId[!duplicated(bed$geneId)],organism = "mmu",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
pdf(paste(argv[3],argv[1],"kegg_25_bar",".pdf",sep = ""),width = 10,height = 6)
barplot(kegg,showCategory = 25)
dev.off()
kegg<-as.data.frame(kegg)
write.table(kegg,paste(argv[3],argv[1],"_anno_kegg.txt",sep = ""),sep = "\t",quote = F)
```



# connect motif target gene to RNA-seq expression --*ChIPseeker*

```R
###program to get target gene fpkm and GO/KEGG CLASSIFICATION
##
##argvs[1]:name,argvs[2]:compare sample ,argvs[3]:anno file path
argvs<-commandArgs(TRUE)
fl12_pl12_fc<-read.csv("/media/skye/W-DATA/motif_target/PLA125vsFL125 result.csv") ##fl12_pl12_foldchange-RNA fpkm
fl12_pl12_fc$logFC<- -fl12_pl12_fc$logFC
fl16_fl12_fc<-read.csv("/media/skye/W-DATA/motif_target/FL165vsFL125 result.csv")##fl12_fl16_foldchange-RNA fpkm
##input
##TF_family is made by cat all TF family target site together(selected by manual)
TF_family<-read.delim(paste(argvs[3],argvs[1],sep = ""),header = F)
TF_family<-TF_family[which(TF_family$V17 != "GENENAME"),]
TF_family<-as.data.frame(as.character(TF_family$V17[!duplicated(TF_family$V17)]))
###extract TARGET gene fpkm 
c<-c()
a=1
for(i in TF_family[,1]){
  if(argvs[2] == "fl12_pl12_fc"){
    if(0< length(which(fl12_pl12_fc$X == i) <2)){
    c[a]<-which(fl12_pl12_fc$X == i)
    a=a+1
    }
  }
  else if (argvs[2] == "fl16_fl12_fc"){
    if(0< length(which(fl16_fl12_fc$X == i) <2)){
      c[a]<-which(fl16_fl12_fc$X == i)
      a=a+1
    }
  }
}
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
if (argvs[2] == "fl12_pl12_fc"){
  all<-fl12_pl12_fc[c,]
  up<-fl12_pl12_fc[which(fl12_pl12_fc$logFC >0 & fl12_pl12_fc$adj.P.Val <= 0.05),]
  down<-fl12_pl12_fc[which(fl12_pl12_fc$logFC <0 & fl12_pl12_fc$adj.P.Val <= 0.05),]
  write.table(up,paste(argvs[3],argvs[1],"_up_gene_fc",sep = ""),quote = F,sep = "\t",row.names = F)
  write.table(down,paste(argvs[3],argvs[1],"_down_gene_fc",sep = ""),quote = F,sep = "\t",row.names = F)
  up_go<-enrichGO(up[,1],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  down_go<-enrichGO(down[,1],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  write.table(as.data.frame(up_go),paste(argvs[3],argvs[1],"_up_gene_GO",sep = ""),quote = F,sep = "\t",row.names = F)
  write.table(as.data.frame(down_go),paste(argvs[3],argvs[1],"_down_gene_GO",sep = ""),quote = F,sep = "\t",row.names = F)
  up_entrez<-select(org.Mm.eg.db,keys = as.character(up[,1]),columns = "ENTREZID",keytype = "SYMBOL")
  down_entrez<-select(org.Mm.eg.db,keys = as.character(down[,1]),columns = "ENTREZID",keytype = "SYMBOL")
  up_kegg<-enrichKEGG(up_entrez[!is.na(up_entrez[,2]),2],keyType = "kegg",organism = "mmu",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  down_kegg<-enrichKEGG(down_entrez[!is.na(down_entrez[,2]),2],keyType = "kegg",organism = "mmu",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  write.table(as.data.frame(up_kegg),paste(argvs[3],argvs[1],"_up_gene_kegg",sep = ""),quote = F,sep = "\t",row.names = F)
  write.table(as.data.frame(down_kegg),paste(argvs[3],argvs[1],"_down_gene_kegg",sep = ""),quote = F,sep = "\t",row.names = F)
}
if(argvs[2] == "fl16_fl12_fc"){
  all<-fl16_fl12_fc[c,]
  up<-fl16_fl12_fc[which(fl16_fl12_fc$logFC >0 & fl16_fl12_fc$adj.P.Val <= 0.05),]
  down<-fl16_fl12_fc[which(fl16_fl12_fc$logFC <0 & fl16_fl12_fc$adj.P.Val <= 0.05),]
  write.table(up,paste(argvs[3],argvs[1],"_up_gene_fc",sep = ""),quote = F,sep = "\t",row.names = F)
  write.table(down,paste(argvs[3],argvs[1],"_down_gene_fc",sep = ""),quote = F,sep = "\t",row.names = F)
  up_go<-enrichGO(up[,1],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  down_go<-enrichGO(down[,1],OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  write.table(as.data.frame(up_go),paste(argvs[3],argvs[1],"_up_gene_GO",sep = ""),quote = F,sep = "\t",row.names = F)
  write.table(as.data.frame(down_go),paste(argvs[3],argvs[1],"_down_gene_GO",sep = ""),quote = F,sep = "\t",row.names = F)
  up_entrez<-select(org.Mm.eg.db,keys = as.character(up[,1]),columns = "ENTREZID",keytype = "SYMBOL")
  down_entrez<-select(org.Mm.eg.db,keys = as.character(down[,1]),columns = "ENTREZID",keytype = "SYMBOL")
  up_kegg<-enrichKEGG(up_entrez[!is.na(up_entrez[,2]),2],keyType = "kegg",organism = "mmu",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  down_kegg<-enrichKEGG(down_entrez[!is.na(down_entrez[,2]),2],keyType = "kegg",organism = "mmu",pvalueCutoff = 0.05,qvalueCutoff = 0.2)
  write.table(as.data.frame(up_kegg),paste(argvs[3],argvs[1],"_up_gene_kegg",sep = ""),quote = F,sep = "\t",row.names = F)
  write.table(as.data.frame(down_kegg),paste(argvs[3],argvs[1],"_down_gene_kegg",sep = ""),quote = F,sep = "\t",row.names = F)
}
## the plot in the loop may get some trouble
png(paste(argvs[3],argvs[1],"_up_GO.png"),width = 1000,height = 600)
barplot(up_go,showCategory=25)
dev.off()
png(paste(argvs[3],argvs[1],"_down_GO.png"),width = 1000,height = 600)
barplot(down_go,showCategory=25)
dev.off()
png(paste(argvs[3],argvs[1],"_up_kegg.png"),width = 1000,height = 600)
barplot(up_kegg,showCategory=25)
dev.off()
png(paste(argvs[3],argvs[1],"_down_kegg.png"),width = 1000,height = 600)
barplot(down_kegg,showCategory=25)
dev.off()
```





## more

https://www.jianshu.com/p/32b2fab75c24

https://nucleoatac.readthedocs.io/en/latest/nucleoatac/

 

