#!/bin/bash

genomefile=/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_bowtie_index/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
loci_file=/bios-store1/chenyc/scripts/sRNA-seq_Analysis_Pipeline/loci/ath_hairpin.loci

input=$1
output=$2
list=$(cat $3)
tag=fastq.gz

#rawdata qc
if [ ! -d "$output/rawdata_qc/" ]; then
        echo
        echo
        echo "[`date`] rawdata QC"
        echo '-----------------------------------------------'
        mkdir -p $output/rawdata_qc
	    fastqc -t 16 --nogroup -o $output/rawdata_qc $input/*."$tag"
        multiqc $output/rawdata_qc -o $output/rawdata_qc/multiqc_result/
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

# read -p "Please check multiqc results and Provide the 3' apapter sequence: " sequence
# wait

sequence=AGATCGGAAGAG
#remove the 3' adapter with cutadapt
if [ ! -d "$output/remove_adapter_3/" ]; then
        echo
        echo
        echo "[`date`] Remove the 3' adapter with cutadapt"
        echo '-----------------------------------------------'
        mkdir -p $output/remove_adapter_3
        myvar=0
        for i in ${list}
            do
            echo "cutadapt -j 24 -a "$sequence" --discard-untrimmed -m 20 -M 24 -o $output/remove_adapter_3/{$i}_trimmed.fq $input/${i}.${tag}"
            cutadapt -j 24 -a "$sequence" --discard-untrimmed -m 20 -M 24 -o $output/remove_adapter_3/"$i"_trimmed.fq $input/"$i"*."$tag" &
			#trim_galore -j 24 -q 20 --length 15 --trim-n --dont_gzip $input/"$i"."$tag" -o $output/remove_adapter_3/ &            
            myvar=$(($myvar + 1 ))
            if [ "$myvar" = "6" ]
                then
                    myvar=0
                    wait
                fi
            done
            wait
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

#remove the 5' adapter with cutadapt
# if [ ! -d "$output/remove_adapter_5/" ]; then
#         echo
#         echo
#         echo "[`date`] Remove the 5' adapter with cutadapt"
#         echo '-----------------------------------------------'
#         mkdir -p $output/remove_adapter_5
#         myvar=0
#         for i in ${list}
#             do
#             echo "cutadapt -j 24 -g ${g} --discard-untrimmed -m 15 -o $output/remove_adapter_5/{$i}_trimmed.fq $input/${i}.${tag}"
#             cutadapt -j 24 -g "$g" --discard-untrimmed -m 15 -o $output/remove_adapter_5/"$i"_trimmed.fq $input/"$i"."$tag" &
#             myvar=$(($myvar + 1 ))
#             if [ "$myvar" = "6" ]
#                 then
#                     myvar=0
#                     wait
#                 fi
#             done
#             wait
#         echo "[`date`] Run complete"
#         echo '-----------------------------------------------'
# fi

#discard low quality reads 
if [ ! -d "$output/filter_lq_reads/" ]; then
        echo
        echo
        echo "[`date`] Discard low quality reads  with FASTX-Toolkit"
        echo '-----------------------------------------------'
        mkdir -p $output/filter_lq_reads
        myvar=0
        for i in ${list}
            do
            echo "fastq_quality_filter -q 20 -p 85 -Q 33 -v -i $output/remove_adapter_3/${i}_trimmed.fq -o $output/filter_lq_reads/${i}_filter.fq"
            fastq_quality_filter -q 20 -p 85 -Q 33 -v -i $output/remove_adapter_3/"$i"_trimmed.fq -o $output/filter_lq_reads/"$i"_filter.fq &
            myvar=$(($myvar + 1 ))
            if [ "$myvar" = "6" ]
                then
                    myvar=0
                    wait
                fi
            done
            wait
	        fastqc -t 16 --nogroup -o $output/filter_lq_reads/*_filter.fq
            multiqc $output/filter_lq_reads -o $output/filter_lq_reads/multiqc_result/
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

#Aligned to the genome
if [ ! -d "$output/map2genome/" ]; then
        echo
        echo
        echo "[`date`] Aligned to the genome with ShortStack"
        echo '-----------------------------------------------'
        mkdir -p $output/map2genome
        myvar=0
        for i in ${list}
            do
            ShortStack --genomefile $genomefile --outdir $output/map2genome/"$i" --mismatches 0 --align_only --nohp --keep_quals --mmap u --bowtie_m 50 --bowtie_cores 24 --readfile $output/filter_lq_reads/"$i"_filter.fq &
            myvar=$(($myvar + 1 ))
            if [ "$myvar" = "2" ]
                then
                    myvar=0
                    wait
                fi
            done
            wait
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/denovo_identification/" ]; then
        echo
        echo
        echo "[ `date` ] Clusters of sRNAs were de novo identified in each library with ShortStack"
        echo '-----------------------------------------------'
        mkdir -p $output/denovo_identification
        myvar=0
        for i in ${list}
                do
                        ShortStack --genomefile $genomefile --outdir $output/denovo_identification/"$i" --mismatches 0 --keep_quals --mincov 2rpm --bowtie_cores 24 --bamfile $output/map2genome/"$i"/"$i"_filter.bam &          
                        myvar=$(($myvar + 1 ))
                        if [ "$myvar" = "2" ]
                         then
                                myvar=0
                                wait
                        fi
                done
                wait
        echo "[ `date` ] Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/loci_drived/" ]; then
        echo
        echo
        echo "[ `date` ] Clusters of sRNAs were identified by known-loci file in each library with ShortStack"
        echo '-----------------------------------------------'
        mkdir -p $output/loci_drived
        myvar=0
        for i in ${list}
                do
                        ShortStack --genomefile $genomefile --locifile $loci_file --outdir $output/loci_drived/"$i" --mismatches 0 --keep_quals --mincov 2rpm --bowtie_cores 24 --bamfile $output/map2genome/"$i"/"$i"_filter.bam &          
                        myvar=$(($myvar + 1 ))
                        if [ "$myvar" = "2" ]
                         then
                                myvar=0
                                wait
                        fi
                done
                wait
        echo "[ `date` ] Run complete"
        echo '-----------------------------------------------'
fi