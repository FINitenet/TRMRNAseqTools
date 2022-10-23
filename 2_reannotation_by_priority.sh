#!/bin/bash
### Create:     
###             2022-5-14
### Update:
###             Null
### Usage:
###             this script [input] [output] [filelist] [priority file] 
###
help() {
    sed -rn 's/^### ?//;T;p' "$0"
}

if [[ $1 = "--help" ]] || [[ $1 = "-h" ]]; then
    help
    exit 1
fi

input=$1
output=$2
list=$(cat $3)
priority=$4

if [ ! -d "$output/Annotation-all.reads" ]; then
    echo
    echo
    echo "[ $(date) ] Annotation reads use all sRNA type"
    mkdir -p "$output/Annotation-all.reads"
    featureCounts -M -O --largestOverlap --fraction -t gene,ncRNA_gene -g biotype -T 24 -R BAM \
        -a /bios-store1/chenyc/scripts/sRNA-seq_Analysis_Pipeline/reference/Arabidopsis_thaliana.TAIR10.53.md.gff3 \
        -o $output/Annotation-all.reads/all.type.annotation $input/*.bam >$output/Annotation-all.reads/log.txt 2>&1
fi

if [ ! -d "$output/Annotation-level" ]; then
    echo
    echo
    echo "[ $(date) ] Sorting by priority files [ filepath: ath_biotype_level.txt] "
    mkdir -p "$output/Annotation-level"
    for i in $list; do
        samtools sort -@ 24 $output/Annotation-all.reads/"$i"_trimmed.bam.featureCounts.bam >$output/Annotation-all.reads/"$i"_sorted.bam && samtools index -@ 24 $output/Annotation-all.reads/"$i"_sorted.bam
        split_featurecounts_bam_by_feature_tags.py -i $output/Annotation-all.reads/"$i"_sorted.bam -o $output/Annotation-level/"$i" -p $priority
    done
    wait
fi

if [ ! -d "$output/doc" ]; then
    echo
    echo
    echo "[ $(date) ] Summary for Plot "
    mkdir -p "$output/doc"
    for i in $(cat $priority); do
        for j in $list; do
            samtools view  $output/Annotation-level/${j}.${i}.bam | awk -v type=$i '{gsub("M","",$6);print $6,type}' | LC_ALL=C sort | LC_ALL=C uniq -c | awk '{OFS="\t";print $3,$2,$1}' >>$output/doc/${j}_summary.txt
        done
    done
    wait
fi

## miRNA
# featureCounts -O --largestOverlap -t miRNA -g Name -s 1 -T 48 \
#     -a /bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath.reannotation.gff3 \
#     -o $output/Annotation-EI/miRNA $output/Annotation-level/*miRNA_primary_transcript.sam >  $output/Annotation-EI/log.miRNA.txt 2>&1

## rRNA
# bioawk -c sam '{if($flag == 16){print ">"$1"\n"revcomp($seq)}else{print ">"$1"\n"$seq}}' $output/Annotation-level/r-T-OX.tRNA.bam > $output/Annotation-level/r-T-OX.tRNA.fa
# bowtie -p 12 -f -v 0 --no-unal -a -x /bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/tRNA_bowtie_index/tRNA $output/Annotation-level/r-T-OX.tRNA.fa -S $output/Annotation-level/r-T-OX.tRNA.sam && \
#     samtools sort -@ 12 $output/Annotation-level/r-T-OX.tRNA.sam > $output/Annotation-level/r-T-OX.tRNA.bam
# samtools index $output/Annotation-level/r-T-OX.tRNA.bam && rm -rf $output/Annotation-level/r-T-OX.tRNA.sam
# samtools idxstats $output/Annotation-level/r-T-OX.tRNA.bam | cut -f 1,3 | sed -e "1 i\ID\tr-T-OX" | sed -e '$d' > $output/Annotation-EI/r-T-OX.tRNA
# paste $output/Annotation-EI/C-T-OX.tRNA $output/Annotation-EI/r-T-OX.tRNA | cut -f 1,2,4 > $output/Annotation-EI/tRNA

if [ ! -d "$output/Annotation-EI" ]; then
    echo
    echo
    echo "[ $(date) ] Expression Quantification use featureCounts "
    mkdir -p $output/Annotation-EI

    array_trsnoRNA=(
        snoRNA
        snRNA
        rRNA
        tRNA
    )

    # array_otherRNA=(
    #     lncRNA
    #     otherRNA
    # )

    for i in $(cat $priority); do
        echo
        echo
        echo "[ $(date) ] Processing: $i"
        if [[ "${array_trsnoRNA[@]}" =~ "${i}" ]]; then

            featureCounts -O --largestOverlap -t ncRNA_gene -g ID -T 48 \
                -a /bios-store1/chenyc/scripts/sRNA-seq_Analysis_Pipeline/reference/ath_trsnoRNA.gff3 \
                -o $output/Annotation-EI/"$i" $output/Annotation-level/*."$i".bam >$output/Annotation-EI/log."$i".txt 2>&1

            sed -i "s/$output\/Annotation-level\///g;s/.$i.bam//g" $output/Annotation-EI/"$i"

        elif [[ "$i" = "miRNA_primary_transcript" ]]; then

            featureCounts -O --largestOverlap -t miRNA -g Name -s 1 -T 48 \
                -a /bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath.reannotation.gff3 \
                -o $output/Annotation-EI/miRNA $output/Annotation-level/*miRNA_primary_transcript.bam >$output/Annotation-EI/log.miRNA.txt 2>&1

            Z-combine_sRNA_by_id_from_featureCounts.py -i $output/Annotation-EI/miRNA -f /bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_mature_bowtie_index/PNAS_miRNA_ver.211213.fa -o $output/Annotation-EI/mature.combined

            sed -i "s/$output\/Annotation-level\///g;s/."$i".miRNA_primary_transcript.bam//g" $output/Annotation-EI/miRNA
            sed -i "s/$output\/Annotation-level\///g;s/."$i".miRNA_primary_transcript.bam//g" $output/Annotation-EI/mature.combined

        elif [[ "$i" = "protein_coding" ]]; then

            featureCounts -O --largestOverlap -t gene -g ID -T 48 \
                -a /bios-store1/chenyc/scripts/sRNA-seq_Analysis_Pipeline/reference/ath_"$i".gff3 \
                -o $output/Annotation-EI/"$i" $output/Annotation-level/*."$i".bam >$output/Annotation-EI/log."$i".txt 2>&1

            sed -i "s/$output\/Annotation-level\///g;s/.$i.bam//g" $output/Annotation-EI/"$i"

        elif [[ "$i" = "otherRNA" ]]; then

            featureCounts -O --largestOverlap -t ncRNA_gene -g Parent -T 48 \
                -a /bios-store1/chenyc/scripts/sRNA-seq_Analysis_Pipeline/reference/ath_"$i".gff3 \
                -o $output/Annotation-EI/"$i" $output/Annotation-level/*."$i".bam >$output/Annotation-EI/log."$i".txt 2>&1

            sed -i "s/$output\/Annotation-level\///g;s/.$i.bam//g" $output/Annotation-EI/"$i"

        elif [[ "$i" = "Unassigned_NoFeatures" ]]; then
            a=($list)
            for ((m = 0; m < ${#a[@]}; m++)); do
                bedtools bamtobed -i $output/Annotation-level/${a[$m]}.Unassigned_NoFeatures.bam -cigar >$output/Annotation-level/${a[$m]}.unknown.bed
            done

            multiIntersectBed -i $output/Annotation-level/*.unknown.bed >$output/Annotation-level/multi.region.noheader.bed

            mergeBed -d 25 -i $output/Annotation-level/multi.region.noheader.bed >$output/Annotation-level/merge.unknown.tmp
            awk '{OFS="\t";print $1,$2,$3,"Cluster_"NR}' $output/Annotation-level/merge.unknown.tmp >$output/Annotation-EI/merge.unknown.bed

            for ((m = 0; m < ${#a[@]}; m++)); do
                bedtools intersect -a $output/Annotation-EI/merge.unknown.bed -b $output/Annotation-level/${a[$m]}.Unassigned_NoFeatures.bam -wa -c >$output/Annotation-EI/${a[$m]}.unknown.summary
            done
        else

            featureCounts -O --largestOverlap -t ncRNA_gene -g ID -T 48 \
                -a /bios-store1/chenyc/scripts/sRNA-seq_Analysis_Pipeline/reference/ath_"$i".gff3 \
                -o $output/Annotation-EI/"$i" $output/Annotation-level/*."$i".bam >$output/Annotation-EI/log."$i".txt 2>&1

            sed -i "s/$output\/Annotation-level\///g;s/.$i.bam//g" $output/Annotation-EI/"$i"

        fi

    done
    wait

fi
