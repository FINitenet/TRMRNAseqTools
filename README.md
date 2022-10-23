# 0_workflow.sh

## v1.0.1

--删除ShortStack预测sRNA的部分，重新编写该模块

## v1.0.2

--删除featureCounts --fraction参数，无该参数时total reads增加1M左右，添加该参数则少0.5M左右，总体不影响。
--增加rawdata QC
--默认clean data长度在10-35nt，ShortStack reads分配以及RNA分型图均使用18-28nt

## v2.0(2021-7-9)

--使用getopt代替getopts，长选项现在可用
--优化分析流程
1.原始下机数据QC
--rawreads number
2.去除3'接头，低质量reads (-q 20 --trim-n)(reads length > 10nt)(map2genome)
--reads numebr
--reads length distribution(/genome_len_dist/)
--extract reads with high count number (top100) (/uniq_reads_count/)
3.保留自定义长度reads (目前的需求是18-28nt) (后续分析则继续使用该长度测序数据）
--reads numebr
4.有效数据（validata），即能够回帖到参考基因组的数据
--reads numebr
--percentage of validata (map number/Step3 number）
5.干净数据（cleandata），即去除trsnoRNA后保留的数据
--reads numebr
--percentage of trsnoRNA (map number/Step4 number）
6.RNA分型 (data source:Step3)
--ShortStack解决multi-mapping reads分配的问题
--featureCounts对reads type进行注释及统计，使用自定义的注释文件(/Annotation-all.reads/)
--各类型reads长度分布(/Annotation-type_len_dis/)
7.pri-miRNA count (data source:Step5)
--featureCounts -O -M --fraction -s 1 --largestOverlap -t miRNA_primary_transcript -g Name 使用miRBase注释文件
-O          Assign reads to all their overlapping meta-features (or features if -f is specified).
-M          Multi-mapping reads will also be counted. For a multi-mapping read, all its reported alignments will be counted. The 'NH' tag in BAM/SAM input is used to detect multi-mapping reads.
--fraction  Assign fractional counts to features. When both '-M' and '-O' are specified, each alignment will carry a fractional count of 1/(x*y).
8.数据可视化
--ShortStack bam2bw
--JBrowse
9.新sRNA loci预测
--不完善，参考文章 [Integrated annotations and analyses of small RNA–producing loci from 47 diverse plants]

## v2.0.1(2021-7-22)

--增加 miRNA hairpin 和 mature 的计数，参考文章 [Bioinformatics Analysis of Small RNAs in Plants Using Next Generation Sequencing Technologies]
--增加 featureCounts --fraction 参数

## v2.0.2(2021-8-2)

--优化mapping流程
1.删除validata
2.去接头仅保留大于10-nt reads
3.map2genome为全长reads map(percentage of validata)
4.map2trsnoRNA，去除ncRNA(percentage of trsnoRNA)
5.re-map2genome，cleandata，保留为bowtie mapping genome methods

## v2.0.2b(2021-9-1)

--本次更新仅做说明，没有功能及参数的改变
--fraction参数除了使得多重注释reads的被按比例分配以外，会导致输出的结果文件中 "XT:Z:" 出现多个注释。虽然对type_len_dist的结果并没有影响(grep时未加-w，已确认过，只是total reads会有浮动)，但是对reads分类时仍然需要注意该问题。

## v2.0.3(2021-9-30)

--优化长度分布统计逻辑，删除不必要统计参数

## V2.0.5(2021-11-17)

--对reads分类重新划分
1.ncRNA --> otherRNA
2.others --> NoFeatures --> otherRNA(only reads distribution used)
(otherRNA=nontranslating_CDS+ncRNA+NoFeatures(others), results same as before just name changed)
3.去除了在Airport11中的TAS基因

## V2.0.5a(2021-12-19)

--增加miRNA re-annotation(mature miRNA merge by same sequence),仅保留该注释结果的BAM文件
--优化多个结果文件，使其更符合R语言作图习惯

## V2.0.5b(2021-12-22)

--使用pigz -p 8 代替gzip

## V2.1.0(2022-1-10)

--在Wang et.al, 2019, PNAS对miRNA/miRNA重新梳理（*Reference_Source/Arabidopsis_Reference/ath_mature_bowtie_index/PNAS_miRNA.fa*）的基础上，修正了部分错误的序列信息，构建了新的miRNA/miRNA*序列文件 （*Reference_Source/Arabidopsis_Reference/ath_mature_bowtie_index/PNAS_miRNA_ver.211213.fa*）。

--基于重新构建的miRNA/miRNA*序列文件，通过回帖到基因组，以及参照miRBase（*Reference_Source/Arabidopsis_Reference/ath.gff3*）hairpin的注释，构建了新的miRNA annotation GFF3文件（*Reference_Source/Arabidopsis_Reference/ath.reannotation.gff3*）。

--在对reads计数时，除对miRBase mature 以及hairpin的计数结果外，增加miRNA re-annotation的计数结果，并以miRNA family的形式重新求和计算。
