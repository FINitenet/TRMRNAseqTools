#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :run.py
# @Time      :2023/07/25 10:51:19
# @Author    :Yuchen@rlab

import glob
import os
import shutil
import subprocess
import sys
from collections import Counter

import pandas as pd
import pysam
from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def data_trim():
    adapter_search = subprocess.run(["/bios-store1/chenyc/scripts/Github_scripts/DNApi/dnapi.py", input_path + "/" + key + ".fastq.gz"], check=True, stdout=subprocess.PIPE, text=True)
    adapter = adapter_search.stdout.split()[0]
    trim_cmd = ["trim_galore", "--fastqc", 
                "-a", adapter,
                "--gzip", "--length", "10", "--max_length", "35", "--consider_already_trimmed", "10",
                "--trim-n", "--suppress_warn", "-j", "12", "-o", trim_path, input_path + "/" + key + ".fastq.gz"]
    subprocess.run(trim_cmd, check=True)

def umi_extract():
    fastq_file = os.path.join(trim_path, key + "_trimmed.fq.gz")
    insert_umi = Counter()
    with pysam.FastqFile(fastq_file) as fh:
        for record in fh:
            # 查询子序列在序列中出现的次数
            if record.sequence.count("AACTGTAGGCACCATCAAT") > 0:
                iu = record.sequence.split("AACTGTAGGCACCATCAAT", 1)
                if 30 >= len(iu[0]) >= 10 and len(iu[1]) == 12 and 'N' not in set(iu[1]):
                    insert_umi.update(["_".join(iu)])

    with bgzf.BgzfWriter(os.path.join(umi_path, key + "_trimmed.fa.gz"), "wb") as out:
        SeqIO.write(
            (SeqRecord(Seq(i_u.split("_")[0]), id="_".join([i_u, str(number)]), description="")
            for i_u, number in insert_umi.items()),
            out,
            "fasta"
        )

def data_length_dist():
    # Read FASTA file and filter sequences based on length
    # must check there is no fai files in the same directory
    if os.path.exists(umi_path + "/" + key + "_trimmed.fa.gz.fai"):
        os.remove(umi_path + "/" + key + "_trimmed.fa.gz.fai")
    fasta_file = os.path.join(umi_path, key + "_trimmed.fa.gz")
    sequence_lengths = []
    with pysam.FastaFile(fasta_file) as fh:
        for sequence_name in fh.references:
            sequence_length = fh.get_reference_length(sequence_name)
            sequence_lengths.append(sequence_length)

        total_count = len(sequence_lengths)
        counts = pd.Series(sequence_lengths).value_counts().sort_index()
        percentages = (counts / total_count) * 100

        table = pd.DataFrame({"Length": counts.index, "Count": counts.values, "Percentage": percentages.values})
        table = table.sort_values(by="Length").reset_index(drop=True)
        table.to_csv(doc_path + "/" + key + "_rawreads_dist.csv", index=False, sep="\t")

def data_mapping():
    bowtie_fq_cmd = [
    "bowtie",
    "-p", "12",
    "-m", "50",
    "-v", "0",
    "--best",
    "--strata",
    "-a",
    "--no-unal",
    "--shmem",
    "-x", genomeindex, 
    umi_path + "/" + key + '_trimmed.fq.gz',
    "-S", mapping_path + '/' + key + '_aligned.sam',
    "--al", mapping_path + '/' + key + '_aligned.fa'
    ]
    bowtie_log = open(doc_path + '/' + key + '.mapresults.txt', 'w')
    log = subprocess.run(bowtie_fq_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    bowtie_log.write(log.stderr)
    bowtie_log.close()

def convert_bam2bigwig():
    samfile = mapping_path + '/' + key + '_aligned.sam'
    bamfile = mapping_path + '/' + key + '_aligned.bam'
    pysam.sort("-o", bamfile, samfile)
    pysam.index(bamfile)
    os.remove(samfile)

    # bw_cmd = [
    # "ShortTracks",
    # "--bamfile", bamfile,
    # "--mode", "readlength",
    # "--stranded"
    # ]
    # subprocess.run(bw_cmd, check=True)
    # bw_files = glob.glob(mapping_path + '/*.bw')
    # for bw_file in bw_files:
    #     shutil.move(bw_file, doc_path)

if __name__ == "__main__":

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    doc_path = sys.argv[2]
    flag = False
    fastq_files = glob.glob(sys.argv[1] + "/*.fastq.gz", recursive=True)
    trim_path = os.path.join(output_path, "2_trim_adapter/")
    umi_path = os.path.join(output_path, "3_extract_umi/")
    mapping_path = os.path.join(output_path, "4_mapping/")
    for key in fastq_files:
        key = os.path.basename(key).replace(".fastq.gz", "")
        print(key)
        # step1: trim adapter
        if not os.path.exists(trim_path):
            os.mkdir(trim_path)
        data_trim()

        # step2: extract umi
        if flag == True:
            if not os.path.exists(umi_path):
                os.mkdir(umi_path)
            umi_extract()
        else:
            umi_path = trim_path

        # step3: mapping
        if not os.path.exists(mapping_path):
            os.mkdir(mapping_path)
        genomeindex = "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_bowtie_index/Arabidopsis_thaliana.TAIR10.dna.toplevel"
        data_mapping()

        # step4: data length distribution
        # data_length_dist()
        convert_bam2bigwig()
