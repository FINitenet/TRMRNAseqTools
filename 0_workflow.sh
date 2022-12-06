#!/bin/bash

set -euo pipefail

helpdoc() {
	cat <<EOF
Description:

    This shellscript is used to run the pipeline to analysis small RNA sequence data

Usage:

    $0 -p <INT> -n <INT> -m <INT> -g <genome> -i <path> -o <path> -b <path> 

Global options:

    -p/--thread         Threads used. [Default:"1"; Integer].
    -i/--input          Input file directory. [Default:".fastq.gz"].
    -o/--output         Output directory. [Default:"./out"]
    -g/--genome	        Currently supported genome: ath, osa, zea, sly [Default:"ath"].
    -b/--batch          List file of input files (Required).
    -n/--min_length     Discard reads that became shorter than length INT because of either quality or adapter trimming. A value of '0' effectively disables this behaviour.[default: 20 bp].
    -m/--max_length     Discard reads that are longer than <INT> bp after trimming.
    -h/--help           Prints this usage information. [Flag]

Alignment options:

    --align_only        If this switch is present, the pipeline run will terminate after the alignment phase with no analysis performed.
    
EOF
}

Settings() {
	cat <<EOF

Basic Informations:

    hostname: $(hostname)
    script: $0
    Version: V1.0.0.20221120_Beta
    Last update: 2022-11-20

Settings:

    indir: $input
    outdir: $output
    readfile: ${list_arr[*]}
    required_cores: $thread
    genome: $genome
    annotation source: /bios-store1/chenyc/scripts/source/$genome
    keep reads min length: $min
    keep reads max length: $max
EOF
}

###  default parameters  ###
min=18
max=28
output=./out
scriptDir=$(dirname $0)
min_thread=1
max_thread=48
genome=ath
tag=fastq.gz
flag=2

if [ $# = 0 ]; then
	helpdoc
	exit 1
fi

###  get input  ###

ARGS=$(getopt -o p:i:o:g:b:n:m:f:h -l thread:,input:,output:,genome:,batch:,min_length:,max_length:,help,align_only -n "$0" -- "$@")
if [ $? -eq 0 ]; then
	eval set -- "${ARGS}"
else
	printf "\033[31mWarning: Please recheck erro information above\033[0m\n"
	exit 1
fi

while true; do
	case "$1" in
	-h | --help)
		helpdoc
		exit 0
		shift
		;;
	-p | --thread)
		thread=$2
		if ((thread != min_thread && thread > max_thread)); then
			printf "\033[31mWarning: thread must be (0,24], default: 1\033[0m\n"
			exit 1
		fi
		shift 2
		;;
	-m | --max_length)
		max=$2
		shift 2
		;;
	-n | --min_length)
		min=$2
		shift 2
		;;
	-i | --input)
		input=$2
		if [ ! -d $input ]; then
			printf "\033[31mWarning: There is no such folder coresponding to $input\033[0m\n"
			exit 1
		fi
		shift 2
		;;
	-o | --output)
		output=$2
		if [ ! -d $output ]; then
			printf "\033[31mWarning: Output directory $output doesn't exist, creating it for you...\033[0m\n"
			mkdir -p $output
		fi
		shift 2
		;;
	-b | --batch)
		LIST=$2
		if [ ! -f $LIST ]; then
			printf "\033[31mWarning: There is no such file coresponding to $LIST\033[0m\n"
			exit 1
		else
			list=$(cat $LIST)
			list_arr=(${list//,/})
		fi
		shift 2
		;;
	-g | --genome)
		genome=$2
		shift 2
		;;
	-f)
		tag=$2
		shift 2
		;;
	--align_only)
		flag=1
		shift
		;;
	--)
		shift
		break
		;;
	*)
		printf "\033[31mWarning: Unknown option: $1 $2\033[0m\n"
		helpdoc
		exit 1
		;;

	esac
done

source $scriptDir/source/$genome

Settings

StartTime=$(date +"%Y-%m-%d %H:%M:%S %Z")

echo
echo "Run Progress and Messages:"
echo '-------------------------------------------------------'
echo "Program starting at $StartTime"
echo '-------------------------------------------------------'
echo

dir1=$output/01.QC
dir2=$output/02.Mapping
dir3=$output/03.Quant
dir4=$output/04.VIs

###### PART1 QC ######
#rawdata qc

if [ ! -d "$dir1/rawdata_qc/" ]; then
	echo
	echo
	echo "[ $(date) ] rawdata QC"
	echo '-----------------------------------------------'
	mkdir -p $dir1/rawdata_qc
	fastqc -t 16 --nogroup -o $dir1/rawdata_qc $input/*.$tag
	multiqc $dir1/rawdata_qc -o $dir1/rawdata_qc/multiqc_result/ -n rawdata_QC
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

#remove the 3′adapter
if [ ! -d "$dir1/trim_adapter/" ]; then
	echo
	echo
	echo "[ $(date) ] Remove the 3′adapter with trim_galore and rename sample name"
	echo '-----------------------------------------------'
	mkdir -p $dir1/trim_adapter
	myvar=0
	for i in ${list}; do
		echo "trim_galore -q 20 --length 10 --trim-n --basename ${i}"
		trim_galore -j 8 -q 20 --basename "$i" --length 10 --consider_already_trimmed 10 --trim-n --fastqc --fastqc_args "--nogroup" --gzip $input/"$i"*.$tag -o $dir1/trim_adapter/ &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait

	multiqc $dir1/trim_adapter/ -o $dir1/trim_adapter/multiqc_result -n trim_adapter_QC
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

#Keep reads length in $min-$maxnt
if [ ! -d "$dir1/$min-$max/" ]; then
	echo
	echo
	echo "[ $(date) ] filter reads -- keep reads length in ${min}-${max} nt"
	echo '-----------------------------------------------'
	mkdir -p $dir1/$min-$max
	myvar=0
	for i in ${list}; do
		echo "seqkit seq -j $thread -w 0 -g -m $min -M $max ${i}_trimmed.fq.gz"
		seqkit seq -j $thread -w 0 -g -m $min -M $max $dir1/trim_adapter/"$i"_trimmed.fq.gz -o $dir1/$min-$max/"$i"_trimmed.fq && pigz -p 8 $dir1/$min-$max/"$i"_trimmed.fq &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait

	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

########################################################################

###### PART2 Mapping ######
#preprocess & mapping

if [ ! -d "$dir2/map2genome/" ]; then
	echo
	echo
	echo "[ $(date) ] filter reads -- mapping to ${genome} genome"
	echo '-----------------------------------------------'
	mkdir -p $dir2/map2genome
	myvar=0
	for i in ${list}; do
		echo "bowtie map2genome ${i}"
		bowtie -p $thread -v 0 --no-unal --best --strata -a \
			-x ${ath[genome_bowtie_index]} $dir1/trim_adapter/"$i"_trimmed.fq.gz \
			-S $dir2/map2genome/"$i"_aligned.sam --al $dir2/map2genome/"$i"_aligned.fastq --un $dir2/map2genome/"$i"_unaligned.fastq >$dir2/map2genome/"$i".mapresults.txt 2>&1 &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait

	fastqc -t 16 $dir2/map2genome/*_aligned.fastq && multiqc $dir2/map2genome/ -o $dir2/map2genome/multiqc_result -n map2genome_QC
	bash $scriptDir/module/7_summary_bowtie.sh $dir2/map2genome $dir2/map2genome $LIST
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

if [ ! -d "$dir2/map2trsnoRNA/" ]; then
	echo
	echo
	echo "[ $(date) ] filter reads -- mapping to ncRNA seq (rRNA, tRNA, etc)"
	echo '-----------------------------------------------'
	mkdir -p $dir2/map2trsnoRNA
	myvar=0
	for i in ${list}; do
		echo "bowtie map2trsnoRNA ${i}"
		bowtie -p $thread -v 0 --no-unal --best --strata -a -x ${ath[trsno_bowtie_index]} $dir2/map2genome/"$i"_aligned.fastq \
			-S $dir2/map2trsnoRNA/"$i".trsno.sam \
			--un $dir2/map2trsnoRNA/"$i"_unaligned.fastq >$dir2/map2trsnoRNA/"$i".mapresults.txt 2>&1 &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

if [ ! -d "$dir2/cleandata/" ]; then
	echo
	echo
	echo "[ $(date) ] filter reads -- re-mapping to ${genome} genome"
	echo '-----------------------------------------------'
	mkdir -p $dir2/cleandata
	myvar=0
	for i in ${list}; do
		echo "bowtie re-mapping to ${genome} genome ${i}"
		bowtie -p $thread -v 0 --no-unal --best --strata -a -m 50 \
			-x ${ath[genome_bowtie_index]} $dir2/map2trsnoRNA/"$i"_unaligned.fastq \
			-S $dir2/cleandata/"$i".remap.sam --al $dir2/cleandata/"$i"_aligned.fastq >$dir2/cleandata/"$i".mapresults.txt 2>&1 && pigz -p 8 $dir2/cleandata/"$i"_aligned.fastq &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait

	fastqc -t 16 $dir2/cleandata/*_aligned.fastq.gz && multiqc $dir2/cleandata/ -o $dir2/cleandata/multiqc_result -n cleandata_QC
	bash $scriptDir/module/7_summary_bowtie.sh $dir2/cleandata $dir2/cleandata $LIST
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

#ShortStack for multi-mapping reads 0-mismatch
if [ ! -d "$dir2/ShortStack/" ]; then
	echo
	echo
	echo "[ $(date) ] Use ShortStack to assign reads"
	echo '-----------------------------------------------'
	mkdir -p $dir2/ShortStack
	myvar=0
	for i in ${list}; do
		echo "${i} aligned to the genome using ShortStack"
		ShortStack --genomefile ${ath[genomefile]} --outdir $dir2/ShortStack/"$i"_ShortStack --align_only --nohp --keep_quals \
			--mmap u --bowtie_m 1000 --ranmax 50 --mismatches 0 \
			--bowtie_cores $thread --readfile $dir1/$min-$max/"$i"_trimmed.fq.gz && mv $dir2/ShortStack/"$i"_ShortStack/Log.txt $dir2/ShortStack/"$i".log &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait
	mv ${output}/02.Mapping/ShortStack/*/*.bam ${output}/02.Mapping/ShortStack/
	rm -rf $dir2/ShortStack/*_ShortStack

	bash $scriptDir/module/7_summary_shortstack.sh $dir2/ShortStack $dir2/ShortStack $LIST
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

#reads length distribution
if [ ! -d "$dir2/genome_len_dist/" ]; then
	echo
	echo
	echo "[ $(date) ] Analysis reads length distribution"
	echo '-----------------------------------------------'
	mkdir -p $dir2/genome_len_dist
	myvar=0
	for i in ${list}; do
		echo "length_distribution ${i}"
		# cat $output/map2genome/"$i"_aligned.fastq | perl -ne '$s=<>;<>;<>;chomp($s);print length($s)."\n";' | LC_ALL=C sort - | LC_ALL=C uniq -c | awk -v sample=$i '{a[NR]=$2;b[NR]=$1;sum+=$1}END{for(i in a){printf "%s\t%d\t%d\t%0.2f\n", sample,a[i],b[i],b[i]/sum*100}}' | LC_ALL=C sort -k 2 -n > $output/genome_len_dist/"$i"_len_dist.txt && pigz -p 8 $dir2/map2genome/"$i"_aligned.fastq &
		python3 $scriptDir/module/size_dist.py $dir2/map2genome/"$i"_aligned.fastq >$dir2/genome_len_dist/"$i"_len_dist.txt && pigz -p 8 $dir2/map2genome/"$i"_aligned.fastq && sed -i "s%$dir2/map2genome/%%g ; s%_aligned.fastq%%g" $dir2/genome_len_dist/"$i"_len_dist.txt &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait

	cat $dir2/genome_len_dist/*_len_dist.txt >$dir2/genome_len_dist/len_dist_summary
	# awk 'BEGIN{OFS="\t";print "sample","length","count","percentage"}{print $0}' $output/genome_len_dist/len_dist_summary.tmp > $output/genome_len_dist/len_dist_summary.txt
	# cat $output/genome_len_dist/${list_arr[0]}_len_dist.txt > $output/genome_len_dist/tmp1
	# for j in ${list_arr[@]:1}
	#         do
	#         join $output/genome_len_dist/${j}_len_dist.txt $output/genome_len_dist/tmp1 > $output/genome_len_dist/tmp2
	#         mv $output/genome_len_dist/tmp2 $output/genome_len_dist/tmp1
	#         done
	#         wait
	# awk -v OFS='\t' '{$1=$1; print}' $output/genome_len_dist/tmp1 > $output/genome_len_dist/len_dist_summary.txt
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

#extract high abundance reads
if [ ! -d "$dir2/uniq_reads_count/" ]; then
	echo
	echo
	echo "[ $(date) ] filter reads -- extract high abundance reads"
	echo '-----------------------------------------------'
	mkdir -p $dir2/uniq_reads_count
	myvar=0
	for i in ${list}; do
		echo "fastx_collapser -v -i $dir2/map2genome/${i}_aligned.fastq.gz -o $dir2/uniq_reads_count/${i}.sort.fasta"
		zcat $dir2/map2genome/"$i"_aligned.fastq.gz | fastx_collapser -v -o $dir2/uniq_reads_count/"$i".sort.fasta &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait
	echo '-----------------------------------------------'
	myvar=0
	for i in ${list}; do
		echo "Converting ${i}.sort.fasta to oneline.fasta"
		sed -n '1{x;d;x};${H;x;s/\n/ /1;s/\n//g;p;b};/^>/{x;s/\n/ /1;s/\n//g;p;b};H' $dir2/uniq_reads_count/"$i".sort.fasta | cut -f 2 -d"-" >$dir2/uniq_reads_count/"$i".temp
		awk '{OFS="\t";print $2,$1}' $dir2/uniq_reads_count/"$i".temp >$dir2/uniq_reads_count/"$i".uniqseq.count.txt

		echo "Split fasta to $min-$max nt"
		for ((a = $min; a < $max + 1; a++)); do awk -v len=$a '{if(length($2)==len) print $2"\t"$1"\t"len}' $dir2/uniq_reads_count/"$i".temp | awk '{if(NR<101) print}' >$dir2/uniq_reads_count/${i}_${a}.temp; done

		echo "Merge table and named "$i".merge.txt"
		cat $dir2/uniq_reads_count/${i}_* >$dir2/uniq_reads_count/${i}.top100.${min}-${max}.txt

	done
	wait
	rm $dir2/uniq_reads_count/*.temp
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

########################################################################

###### Switch ######
if [ $flag == 1 ]; then
	echo '-----------------------------------------------'
	echo "[ $(date) ] Aligment Run complete. END THE PIPELINE"
	echo '-----------------------------------------------'
	exit
else
	echo '-----------------------------------------------'
	echo "[ $(date) ] Aligment Run complete. Continue the PIPELINE for Annotation and Visualization"
	echo '-----------------------------------------------'
fi

########################################################################

###### PART3 Quant ######
#annotation
#featureCounts v2.0.1

if [ ! -d "$dir3/Annotation-all.reads" ]; then
	echo
	echo
	echo "[ $(date) ] Annotating sRNA reads"
	echo '-----------------------------------------------'
	mkdir -p $dir3/Annotation-all.reads

	echo "featureCounts -M -O --largestOverlap --fraction -t gene,ncRNA_gene -g biotype -T $thread -R BAM -a ${ath[all_type_annotation]} -o $dir3/Annotation-all.reads/all.type.annotation"

	featureCounts -M -O --largestOverlap --fraction -t gene,ncRNA_gene -g biotype -T $thread -R BAM -a ${ath[all_type_annotation]} -o $dir3/Annotation-all.reads/all.type.annotation $dir2/ShortStack/*.bam >$dir3/Annotation-all.reads/log.txt 2>&1

	echo "Editing Summary Files"

	awk 'NR>1{$2=$3=$4=$5=$6=null;print}' $dir3/Annotation-all.reads/all.type.annotation | awk '{for(i=1;i<NF;i++)printf("%s\t",$i);print $NF}' >$dir3/Annotation-all.reads/summary.txt

	sed -i "1 s%$dir2/ShortStack/%%g ; s%_trimmed.bam%%g" $dir3/Annotation-all.reads/summary.txt

	echo "Calculating the number of others reads"

	grep 'Unassigned_NoFeatures' $dir3/Annotation-all.reads/all.type.annotation.summary | awk '{OFS="\t";gsub("Unassigned_NoFeatures","NoFeatures",$1);print $0}' >>$dir3/Annotation-all.reads/summary.txt

	#       awk '{for(i=2;i<=NF;i++)a[i]+=$i;print}END{OFS="\t";printf "part\t";for(j=2;j<=NF;j++)printf "%.2f\t", a[j]; print""}' $output/Annotation-all.reads/summary.txt | grep '^Geneid\|^part' > $output/Annotation-all.reads/others.txt

	#       awk '{for(i=2;i<=NF;i++)a[i]+=$i;print}END{OFS="\t";printf "total\t";for(j=2;j<=NF;j++)printf "%.2f\t", a[j]; print""}' $output/Annotation-all.reads/all.type.annotation.summary |  grep '^total' >>  $output/Annotation-all.reads/others.txt

	#       awk -F'\t' 'BEGIN{printf "others\t"} NR == 2{ for (i = 2; i < NF; i++) hash[i] = $i;}NR > 2{ for (i = 2; i < NF; i++){printf("%.2f\t", $i - hash[i]);hash[i] = $i;}printf("\n");}' $output/Annotation-all.reads/others.txt >> $output/Annotation-all.reads/summary.txt
	myvar=0
	for i in ${list}; do
		samtools index -@ $thread $dir3/Annotation-all.reads/"$i"_trimmed.bam.featureCounts.bam &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done

	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

if [ ! -d "$dir3/Annotation-pri-miRNA" ]; then
	echo
	echo
	echo "[ $(date) ] Annotating miRNA reads"
	echo '-----------------------------------------------'
	mkdir -p $dir3/Annotation-pri-miRNA
	echo "featureCounts -O -M --largestOverlap -t miRNA_primary_transcript -g Name -T $thread -a ${ath[miRNA_annotation]} -o $dir3/Annotation-pri-miRNA/miRNA.hairpin.annotation"
	featureCounts -O -M --largestOverlap -t miRNA_primary_transcript -g Name -T $thread -a ${ath[miRNA_annotation]} -o $dir3/Annotation-pri-miRNA/miRNA.hairpin.annotation $dir2/ShortStack/*.bam >$dir3/Annotation-pri-miRNA/hairpin.log 2>&1
	echo "featureCounts -O -M --largestOverlap -t miRNA -g Name -T $thread -a ${ath[miRNA_annotation]} -o $dir3/Annotation-pri-miRNA/miRNA.mature.annotation"
	featureCounts -O -M --largestOverlap -t miRNA -g Name -T $thread -a ${ath[miRNA_annotation]} -o $dir3/Annotation-pri-miRNA/miRNA.mature.annotation $dir2/ShortStack/*.bam >$dir3/Annotation-pri-miRNA/mature.log 2>&1
	echo "featureCounts -O -M --largestOverlap -t miRNA -g Name -T $thread -R BAM -a ${ath[miRNA_re-annotation]} -o $dir3/Annotation-pri-miRNA/miRNA.re.annotation"
	featureCounts -O -M --largestOverlap -t miRNA -g Name -T $thread -R BAM -a ${ath[miRNA_re-annotation]} -o $dir3/Annotation-pri-miRNA/miRNA.re.annotation $dir2/ShortStack/*.bam >$dir3/Annotation-pri-miRNA/re.log 2>&1

	sed -i "s%$dir2/ShortStack/%%g ; s%_trimmed.bam%%g" $dir3/Annotation-pri-miRNA/miRNA.mature.annotation
	sed -i "s%$dir2/ShortStack/%%g ; s%_trimmed.bam%%g" $dir3/Annotation-pri-miRNA/miRNA.hairpin.annotation
	sed -i "s%$dir2/ShortStack/%%g ; s%_trimmed.bam%%g" $dir3/Annotation-pri-miRNA/miRNA.re.annotation

	python $scriptDir/module/combine_sRNA_by_id_from_featureCounts.py -i $dir3/Annotation-pri-miRNA/miRNA.re.annotation -f /bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_mature_bowtie_index/PNAS_miRNA_ver.211213.fa -o $dir3/Annotation-pri-miRNA/miRNA.re.combine
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

if [ ! -f "$scriptDir/source/${genome}_biotype.txt" ]; then
	awk 'NR>1{gsub("others","Unassigned_NoFeatures",$1);print $1}' $dir3/Annotation-all.reads/summary.txt >$scriptDir/source/${genome}_biotype.txt && type=$(cat $scriptDir/source/${genome}_biotype.txt)
else
	type=$(cat $scriptDir/source/${genome}_biotype.txt)
fi

if [ ! -d "$dir3/Annotation-type_len_dis/" ]; then
	echo
	echo
	echo '-----------------------------------------------'
	mkdir -p $dir3/Annotation-type_len_dis
	myvar=0
	for i in ${list}; do
		echo "[ $(date) ] Editing Summary Files for barplot"
		mkdir -p $output/Annotation-type_len_dis/"$i"
		for ((a=$min;a<$max+1;a++))
			do
			samtools view $output/Annotation-all.reads/"$i"_trimmed.bam.featureCounts.bam | awk -v len=$a '{gsub("M","",$6);if($6 == len) print $20,$22}' | sed -e 's/XS:Z:Assigned//g ; s/XT:Z:\|XS:Z://g ; s/Unassigned_NoFeatures/otherRNA/g ; s/ //g ; s/,/\n/g' > $output/Annotation-type_len_dis/"$i"/"$a".txt
			LC_ALL=C sort $output/Annotation-type_len_dis/"$i"/"$a".txt | LC_ALL=C uniq -c | awk -v len=$a '{printf("%s\t%s\t%s\n",len,$2,$1)}' >> $output/Annotation-type_len_dis/"$i"/"$i".list
			done
		awk '{n+=$3}END{print n}' $output/Annotation-type_len_dis/"$i"/"$i".list > $output/Annotation-type_len_dis/"$i"/total.reads && total=$(cat $output/Annotation-type_len_dis/"$i"/total.reads)
		awk -v n=$total '{printf("%s\t%s\t%s\t%.4f\n",$1,$2,$3,($3/n*100))}' $output/Annotation-type_len_dis/"$i"/"$i".list > $output/Annotation-type_len_dis/"$i"/"$i".summary
		python3 $scriptDir/module/calc_feature_count_by_featureCounts_tags.py -i $dir3/Annotation-all.reads/"$i"_trimmed.bam.featureCounts.bam -o $dir3/Annotation-type_len_dis/"$i".avg.summary &
	done
	rm -rf $output/Annotation-type_len_dis/*/*.txt
	rm -rf $output/Annotation-type_len_dis/*/*.list
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

########################################################################

###### PART4 VIS ######
#Data visualization

if [ ! -d "$dir4/deeptools_bam2bw/" ]; then
	echo
	echo
	echo "[ $(date) ] Converting BAM to bigweg for igv"
	echo '-----------------------------------------------'
	mkdir -p $dir4/deeptools_bam2bw
	myvar=0
	for i in ${list}; do
		echo "samtools index -@ $thread $dir2/ShortStack/"$i"_trimmed.bam"
		samtools index -@ $thread $dir2/ShortStack/"$i"_trimmed.bam &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait
	myvar=0
	for i in ${list}; do
		echo "bamCoverage -b $dir2/ShortStack/"$i"_trimmed.bam --filterRNAstrand forward --binSize 1 --normalizeUsing CPM -o $dir4/deeptools_bam2bw/"$i"_aligned.R.bw"
		bamCoverage --numberOfProcessors max/2 -b $dir2/ShortStack/"$i"_trimmed.bam --filterRNAstrand forward --binSize 1 --normalizeUsing CPM -o $dir4/deeptools_bam2bw/"$i"_aligned.R.bw &
		echo "bamCoverage -b $dir2/ShortStack/"$i"_trimmed.bam --filterRNAstrand reverse --binSize 1 --normalizeUsing CPM -o $dir4/deeptools_bam2bw/"$i"_aligned.F.bw"
		bamCoverage --numberOfProcessors max/2 -b $dir2/ShortStack/"$i"_trimmed.bam --filterRNAstrand reverse --binSize 1 --normalizeUsing CPM -o $dir4/deeptools_bam2bw/"$i"_aligned.F.bw &
		myvar=$(($myvar + 1))
		if [ "$myvar" = "6" ]; then
			myvar=0
			wait
		fi
	done
	wait
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

if [ ! -d "$dir4/deeptools_plot/" ]; then
	echo
	echo
	echo "[ $(date) ] Data QC visualization"
	echo '-----------------------------------------------'
	mkdir -p $dir4/deeptools_plot
	echo "[ $(date) ] multiBamSummary bins --binSize 100"
	multiBamSummary bins --bamfiles $dir2/ShortStack/*_trimmed.bam \
		--minMappingQuality 20 --smartLabels --binSize 100 --numberOfProcessors max \
		-o $dir4/deeptools_plot/readCounts.npz --outRawCounts $dir4/deeptools_plot/readCounts.tab
	echo "[ $(date) ] plotCorrelation --corMethod spearman --whatToPlot scatterplot"
	plotCorrelation -in $dir4/deeptools_plot/readCounts.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot scatterplot -o $dir4/deeptools_plot/scatterplot_SpearmanCorr.png --outFileCorMatrix $dir4/deeptools_plot/scatterplot_SpearmanCorr.tab
	echo "[ $(date) ] plotCorrelation --corMethod spearman --whatToPlot heatmap"
	plotCorrelation -in $dir4/deeptools_plot/readCounts.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap -o $dir4/deeptools_plot/heatmap_SpearmanCorr.png --outFileCorMatrix $dir4/deeptools_plot/heatmap_SpearmanCorr.tab --colorMap RdYlBu --plotNumbers
	echo "[ $(date) ] plotCorrelation --corMethod pearson --whatToPlot scatterplot"
	plotCorrelation -in $dir4/deeptools_plot/readCounts.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" --whatToPlot scatterplot -o $dir4/deeptools_plot/scatterplot_pearsonCorr.png --outFileCorMatrix $dir4/deeptools_plot/scatterplot_pearsonCorr.tab
	echo "[ $(date) ] plotCorrelation --corMethod pearson --whatToPlot heatmap"
	plotCorrelation -in $dir4/deeptools_plot/readCounts.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" --whatToPlot heatmap -o $dir4/deeptools_plot/heatmap_pearsonCorr.png --outFileCorMatrix $dir4/deeptools_plot/heatmap_pearsonCorr.tab --colorMap RdYlBu --plotNumbers
	echo "[ $(date) ] Run complete"
	echo '-----------------------------------------------'
fi

EndTime=$(date +"%Y-%m-%d %H:%M:%S %Z")
echo "Completed at $EndTime"
st=$(date -d "$StartTime" +%s)
et=$(date -d "$EndTime" +%s)
sumTime=$(($et - $st))
echo "Total time is : $sumTime seconds."
