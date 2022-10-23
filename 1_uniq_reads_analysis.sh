list=$(cat filenames.txt)
LIST=filenames.txt
min=10
max=35
genome_bowtie_index=/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_bowtie_index/Arabidopsis_thaliana.TAIR10.dna.toplevel
all_type_annotation=/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/Arabidopsis_thaliana.TAIR10.51.all.gff3
scriptDir=/bios-store1/chenyc/scripts/sRNA-seq_Analysis_Pipeline
output=/bios-store1/chenyc/DataBase_TRM-sRNA-seq/rdr/Project_chensusu_220225N/uniq_analysis





if [ ! -d "$output/map2genome/" ]; then
        echo
        echo
        echo "[`date`] filter reads -- mapping to ${genome} genome"
        echo '-----------------------------------------------'
        mkdir -p $output/map2genome
        myvar=0
        for i in ${list}
                do
                        echo "bowtie map2genome ${i}"
                        bowtie -p 12 -v 0 -f --no-unal --best --strata -m 50\
                                -x $genome_bowtie_index uniq_reads_count/"$i".sort.fasta \
                                -S $output/map2genome/"$i"_aligned.sam --al $output/map2genome/"$i"_aligned.fasta --un $output/map2genome/"$i"_unaligned.fasta  > $output/map2genome/"$i".mapresults.txt 2>&1 &
                        myvar=$(($myvar + 1 ))
                        if [ "$myvar" = "6" ]
                        then
                                myvar=0
                                wait
                        fi
                done
                wait
                bash $scriptDir/module/7_summary_bowtie.sh $output/map2genome $output/map2genome $LIST
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

if  [ ! -d "$output/ShortStack/" ]; then
        echo
        echo
        echo "[ `date` ] Use ShortStack to assign reads"
        echo '-----------------------------------------------'
        mkdir -p $output/ShortStack
        myvar=0
        for i in ${list}
                do
                        echo "${i} aligned to the genome using ShortStack"
                        ShortStack --genomefile ${genome_bowtie_index}.fa --outdir $output/ShortStack/"$i"_ShortStack --align_only --nohp --keep_quals \
                            --mmap u --bowtie_m 50 --ranmax 50 --mismatches 0 \
                            --bowtie_cores 12 --readfile uniq_reads_count/"$i".sort.fasta && mv $output/ShortStack/"$i"_ShortStack/Log.txt $output/ShortStack/"$i".log &
                        myvar=$(($myvar + 1 ))
                        if [ "$myvar" = "6" ]
                        then
                              myvar=0
                              wait
                        fi
                done
                wait
                mv ${output}/ShortStack/*/*.bam ${output}/ShortStack/
                rm -rf $output/ShortStack/*_ShortStack
                bash $scriptDir/module/7_summary_shortstack.sh $output/ShortStack $output/ShortStack $LIST
        echo "[ `date` ] Run complete"
        echo '-----------------------------------------------'
fi


if [ ! -d "$output/genome_len_dist/" ]; then
        echo
        echo
        echo "[`date`] Analysis reads length distribution"
        echo '-----------------------------------------------'
        mkdir -p $output/genome_len_dist
        myvar=0
        for i in ${list}
                do
                        echo "length_distribution ${i}"
                        seqkit fx2tab -j 12 -l -n $output/map2genome/"$i"_aligned.fasta | cut -f 2  | sort - | uniq -c | awk -v sample=$i '{a[NR]=$2;b[NR]=$1;sum+=$1}END{for(i in a){printf "%s\t%d\t%d\t%0.2f\n", sample,a[i],b[i],b[i]/sum*100}}' | sort -k 1 -n > $output/genome_len_dist/"$i"_len_dist.txt &
                        myvar=$(($myvar + 1 ))
                        if [ "$myvar" = "6" ]
                        then
                                myvar=0
                                wait
                        fi
                done
                wait
        cat $output/genome_len_dist/*_len_dist.txt > $output/genome_len_dist/len_dist_summary.tmp
        awk 'BEGIN{OFS="\t";print "sample","length","count","percentage"}{print $0}' $output/genome_len_dist/len_dist_summary.tmp > $output/genome_len_dist/len_dist_summary.txt
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi

if [ ! -d "$output/Annotation-type_len_dis/" ]; then
        echo
        echo
        echo '-----------------------------------------------'
        mkdir -p $output/Annotation-type_len_dis
        myvar=0
        for i in ${list}
                do
                echo "[`date`] Editing Summary Files for barplot"
                mkdir -p $output/Annotation-type_len_dis/"$i"
                featureCounts -O --largestOverlap --fraction -t gene,ncRNA_gene -g biotype -T 12 -R BAM -a $all_type_annotation -o $output/Annotation-type_len_dis/all.type.annotation $output/ShortStack/"$i".sort.bam >> $output/Annotation-type_len_dis/log.txt 2>&1
                for ((a=$min;a<$max+1;a++))
                        do
                        samtools view $output/Annotation-type_len_dis/"$i".sort.bam.featureCounts.bam | awk -v len=$a '{gsub("M","",$6);if($6 == len) print $20,$22}' | sed -e 's/XS:Z:Assigned//g ; s/XT:Z:\|XS:Z://g ; s/Unassigned_NoFeatures/otherRNA/g ; s/ //g ; s/,/\n/g' > $output/Annotation-type_len_dis/"$i"/"$a".txt
                        sort $output/Annotation-type_len_dis/"$i"/"$a".txt | uniq -c | awk -v len=$a '{printf("%s\t%s\t%s\n",len,$2,$1)}' >> $output/Annotation-type_len_dis/"$i"/"$i".list
                        done
                awk '{n+=$3}END{print n}' $output/Annotation-type_len_dis/"$i"/"$i".list > $output/Annotation-type_len_dis/"$i"/total.reads && total=$(cat $output/Annotation-type_len_dis/"$i"/total.reads)
                awk -v n=$total '{printf("%s\t%s\t%s\t%.4f\n",$1,$2,$3,($3/n*100))}' $output/Annotation-type_len_dis/"$i"/"$i".list > $output/Annotation-type_len_dis/"$i"/"$i".summary 
        done
        echo "[`date`] Run complete"
        echo '-----------------------------------------------'
fi