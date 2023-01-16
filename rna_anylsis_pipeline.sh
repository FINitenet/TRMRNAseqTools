mkdir 2_rawdata_qc
fastqc -t 16 -o 2_rawdata_qc 1_rawdata/*.fastq.gz

mkdir 3_trimmed
for i in se1 se2 sh1 sh2 ss1 ss2 wt1 wt2
    do
    trim_galore -j 12 --gzip -o 3_trimmed --paired 1_rawdata/"$i"*_R1_*.fastq.gz 1_rawdata/"$i"*_R2_*.fastq.gz
    done

mkdir 4_mapping
for i in C-F C-S d-F -d-S
    do
    hisat2 -p 6 -t -x /bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_hisat_index/Arabidopsis_thaliana.TAIR10.dna.toplevel -1 3_trimmed/"$i"*_R1_*.fq.gz -2 3_trimmed/"$i"*_R2_*.fq.gz -S 4_mapping/"$i".sam --summary-file 4_mapping/"$i"_mapping_summary.txt --new-summary --no-unal --dta &
    done

mkdir 5_featureCounts
featureCounts -p -T 12 -t gene,ncRNA_gene -g ID --extraAttributes biotype -a /bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/Arabidopsis_thaliana.TAIR10.53.chr.gff3 -o 5_featureCounts/featureCounts.count 4_mapping/*.sam
sed -i 's/4_mapping\///g; s/\.sam//g; s/^gene://g' 5_featureCounts/featureCounts.count