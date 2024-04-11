#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :run.py
# @Time      :2023/07/25 10:51:19
# @Author    :Yuchen@rlab
# @Version   :2.0.20240325

import glob
import os
import shutil
import subprocess
import sys
from collections import Counter
import argparse

import pandas as pd
import pysam
from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def data_trim(ADAPTER, TRIM_FLAG, TRIM_PATH, fastq_file):
    """
    Trim adapter and low quality reads,
    flag = 1 means trim adapter use self adapter, 
    flag = 2 means auto detect adapter by trim_galore,(default) 
    flag = 3 means auto detect adapter by dnapi.py
    """
    
    if ADAPTER and TRIM_FLAG == 1:
        trim_cmd = ["trim_galore", "--fastqc", 
                "-a", ADAPTER,
                "--gzip", "--length", "10", "--consider_already_trimmed", "10",
                "--trim-n", "--suppress_warn", "-j", "12", "-o", TRIM_PATH, fastq_file]
    elif ADAPTER == None and TRIM_FLAG == 2:
        trim_cmd = ["trim_galore", "--fastqc", 
                "--length", "18", "--max_length", "28", 
                "--gzip", "--consider_already_trimmed", "10",
                "--trim-n", "--suppress_warn", "-j", "12", "-o", TRIM_PATH, fastq_file]
    elif ADAPTER == None and TRIM_FLAG == 3:
        adapter_search = subprocess.run(["/bios-store1/chenyc/scripts/Github_scripts/DNApi/dnapi.py", fastq_file], check=True, stdout=subprocess.PIPE, text=True)
        adapter = adapter_search.stdout.split()[0]
        trim_cmd = ["trim_galore", "--fastqc", 
                    "-a", adapter,
                    "--gzip", "--length", "10", "--consider_already_trimmed", "10",
                    "--trim-n", "--suppress_warn", "-j", "12", "-o", TRIM_PATH, fastq_file]
    subprocess.run(trim_cmd, check=True)

def data_umi_extract(BASENAME, TRIM_PATH, UMI_PATH):
    fastq_file = os.path.join(TRIM_PATH, BASENAME + "_trimmed.fq.gz")
    insert_umi = Counter()
    with pysam.FastqFile(fastq_file) as fh:
        for record in fh:
            # 查询子序列在序列中出现的次数
            if record.sequence.count("AACTGTAGGCACCATCAAT") > 0:
                iu = record.sequence.split("AACTGTAGGCACCATCAAT", 1)
                if 30 >= len(iu[0]) >= 10 and len(iu[1]) == 12 and 'N' not in set(iu[1]):
                    insert_umi.update(["_".join(iu)])

    with bgzf.BgzfWriter(os.path.join(UMI_PATH, BASENAME + "_trimmed.fa.gz"), "wb") as out:
        SeqIO.write(
            (SeqRecord(Seq(i_u.split("_")[0]), id="_".join([i_u, str(number)]), description="")
            for i_u, number in insert_umi.items()),
            out,
            "fasta"
        )
    
    bowtie_fa_cmd = [
    "bowtie",
    "-f",
    "-p", "12",
    "-m", "50",
    "-v", "0",
    "--best",
    "--strata",
    "-a",
    "--no-unal",
    "--shmem",
    "-x", genomeindex, 
    UMI_PATH + "/" + BASENAME + '_trimmed.fa.gz',
    "-S", UMI_PATH + '/' + BASENAME + '_aligned.sam',
    "--al", UMI_PATH + '/' + BASENAME + '_aligned.fa'
    ]
    bowtie_log = open(UMI_PATH + '/' + BASENAME + '.mapresults.txt', 'w')
    log = subprocess.run(bowtie_fa_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    bowtie_log.write(log.stderr)
    bowtie_log.close()

def data_mapping(BASENAME, TRIM_PATH, MAPPING_PATH, DOC_PATH):
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
    TRIM_PATH + "/" + BASENAME + '_trimmed.fq.gz',
    "-S", MAPPING_PATH + '/' + BASENAME + '_aligned.sam',
    "--al", MAPPING_PATH + '/' + BASENAME + '_aligned.fq'
    ]
    bowtie_log = open(MAPPING_PATH + '/' + BASENAME + '.mapresults.txt', 'w')
    log = subprocess.run(bowtie_fq_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    bowtie_log.write(log.stderr)
    bowtie_log.close()

    sequence_lengths = []
    with pysam.FastaFile(MAPPING_PATH + '/' + BASENAME + '_aligned.fq') as fh:
        for sequence_name in fh.references:
            sequence_length = fh.get_reference_length(sequence_name)
            sequence_lengths.append(sequence_length)

        total_count = len(sequence_lengths)
        counts = pd.Series(sequence_lengths).value_counts().sort_index()
        percentages = (counts / total_count) * 100

        table = pd.DataFrame({"Length": counts.index, "Count": counts.values, "Percentage": percentages.values})
        table = table.sort_values(by="Length").reset_index(drop=True)
        table.to_csv(DOC_PATH + "/" + BASENAME + "_map2genome_dist.csv", index=False, sep="\t")


def convert_bam2bigwig(MAPPING_PATH, BW_PATH):
    if not os.path.exists(BW_PATH):
        os.mkdir(BW_PATH)
    samfile = MAPPING_PATH + '/' + BASENAME + '_aligned.sam'
    bamfile = MAPPING_PATH + '/' + BASENAME + '_aligned.bam'
    pysam.sort("-o", bamfile, samfile)
    pysam.index(bamfile)
    os.remove(samfile)

    bw_cmd = [
    "/home/chenyc/anaconda3/envs/ShortStack4/bin/ShortTracks",
    "--bamfile", bamfile,
    "--mode", "readlength",
    "--stranded"
    ]
    subprocess.run(bw_cmd, check=True)
    bw_files = glob.glob(MAPPING_PATH + '/*.bw')
    for bw_file in bw_files:
        shutil.move(bw_file, BW_PATH)

def data_length_dist(INPUT_FILE, DOC_PATH):
    # Read FASTA file and filter sequences based on length
    # must check there is no fai files in the same directory
    
    sequence_lengths = []
    with pysam.FastaFile(INPUT_FILE) as fh:
        for sequence_name in fh.references:
            sequence_length = fh.get_reference_length(sequence_name)
            sequence_lengths.append(sequence_length)

        total_count = len(sequence_lengths)
        counts = pd.Series(sequence_lengths).value_counts().sort_index()
        percentages = (counts / total_count) * 100

        table = pd.DataFrame({"Length": counts.index, "Count": counts.values, "Percentage": percentages.values})
        table = table.sort_values(by="Length").reset_index(drop=True)
        table.to_csv(DOC_PATH + "/" + BASENAME + "_rawreads_dist.csv", index=False, sep="\t")


# Argument parsing / help message / version
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))
parser.add_argument("-v", "--version", action="version",
                    version='%(prog)s v2.0.20240325')
parser.add_argument("-i", "--inputdir", type= str, default= os.getcwd(),
                    help="path to input directory, default is current directory")
parser.add_argument("-o", "--outdir", type= str, default= os.getcwd(),
                    help="path to output directory, default is current directory")
parser.add_argument("-g", "--genome_index", type= str, default= "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_bowtie_index/Arabidopsis_thaliana.TAIR10.dna.toplevel",
                    help="path to genome index file")
parser.add_argument("-a", "--adapter", type= str, default= None)
parser.add_argument("-t", "--trim_flag", type= int, default= 3)
parser.add_argument("-u", "--umi_flag", type= bool, default= True)
parser.add_argument("-s", "--suffix", type= str, default= ".fastq.gz")


args = parser.parse_args()


if __name__ == "__main__":

    INPUTDIR = args.inputdir
    OUTDIR = args.outdir
    DOC_PATH = args.outdir

    TRIM_FLAG = 3
    UMI_FLAG = args.umi_flag
    ADAPTER = None
    suffix = args.suffix

    genomeindex = args.genome_index

    TRIM_PATH = os.path.join(OUTDIR, "2_trim_adapter/")
    UMI_PATH = os.path.join(OUTDIR, "3_extract_umi/")
    MAPPING_PATH = os.path.join(OUTDIR, "4_mapping/")
    BW_PATH = os.path.join(OUTDIR, "5_bam2bigwig/")

    fastq_files = glob.glob(INPUTDIR + "/*" + suffix, recursive=True)

    for fastq_file in fastq_files:
        BASENAME = os.path.basename(fastq_file).replace(suffix, "")
        print(BASENAME)
        # step1: trim adapter
        if not os.path.exists(TRIM_PATH):
            os.mkdir(TRIM_PATH)
        data_trim(ADAPTER, TRIM_FLAG, TRIM_PATH, fastq_file)

        # step2: extract umi
        if UMI_FLAG == True:
            if not os.path.exists(UMI_PATH):
                os.mkdir(UMI_PATH)
            data_umi_extract(BASENAME, TRIM_PATH, UMI_PATH)

        # step3: mapping to genome
        if not os.path.exists(MAPPING_PATH):
            os.mkdir(MAPPING_PATH)
        data_mapping(BASENAME, TRIM_PATH, MAPPING_PATH, DOC_PATH)
        convert_bam2bigwig(MAPPING_PATH, BW_PATH)
