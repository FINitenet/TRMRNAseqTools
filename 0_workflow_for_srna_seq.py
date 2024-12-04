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
    flag = 1 means trim adapter use self adapter
    flag = 2 means auto detect adapter by trim_galore(default) 
    flag = 3 means auto detect adapter by dnapi.py
    """
    
    if ADAPTER and TRIM_FLAG == 1:
        trim_cmd = ["trim_galore", "--fastqc", "--fastqc_args", "-t 16 --nogroup",
                "-a", ADAPTER, "--basename", BASENAME,
                "--gzip", "--length", "10", "--max_length", "30", "--consider_already_trimmed", "10",
                "--trim-n", "--suppress_warn", "-j", "12", "-o", TRIM_PATH, fastq_file]
    elif ADAPTER == None and TRIM_FLAG == 2:
        trim_cmd = ["trim_galore", "--fastqc", "--fastqc_args", "-t 16 --nogroup",
                "--length", "10", "--basename", BASENAME,
                "--gzip", "--consider_already_trimmed", "10",
                "--trim-n", "--suppress_warn", "-j", "12", "-o", TRIM_PATH, fastq_file]
    elif ADAPTER == None and TRIM_FLAG == 3:
        adapter_search = subprocess.run(["/bios-store1/chenyc/scripts/Github_scripts/DNApi/dnapi.py", fastq_file], check=True, stdout=subprocess.PIPE, text=True)
        adapter = adapter_search.stdout.split()[0]
        trim_cmd = ["trim_galore", "--fastqc", "--fastqc_args", "-t 16 --nogroup",
                    "-a", adapter, "--basename", BASENAME, 
                    "--gzip", "--length", "10", "--consider_already_trimmed", "10",
                    "--trim-n", "--suppress_warn", "-j", "12", "-o", TRIM_PATH, fastq_file]
    subprocess.run(trim_cmd, check=True)

def find_adapter(fqfile, prefix, n):
    """Given path to fastq/fastq.gz and a prefix, find best 3' adapter.

    Inputs:
    - fqfile : path to fastq or gzip-compressed fastq file. 
       Can also be .fa, .fa.gz, .fasta, .fasta.gz
    - prefix : sequence to use a 'prefix' .. a common microRNA
    - n : Integer for use by tqdm for positioning progress tracker

    Output:
    - adapter : Up to 20nt long best adapter sequence.
    
    Depends on FastqGeneralIterator from Bio.SeqIO.QualityIO
    """
    prefixcount = 0
    seqcount = 0
    aseqs = {}
    format = ''

    # Check to see if it is fasta. If not, it must be fastq.
    possible_fa_exts = ['.fa', '.fasta', '.fa.gz', '.fasta.gz']
    # compile a regex of the valid exts for fasta files
    pattern = '$|'.join(possible_fa_exts) + '$'
    pattern = pattern.replace('.', '\.')
    fa_exts = re.compile(pattern)
    if fa_exts.search(fqfile):
        format = 'fasta'
    else:
        format = 'fastq'
    
    if(re.search("gz$", fqfile)):
        fqhandle = gzip.open(fqfile, "rt")
    else:
        fqhandle = open(fqfile, "rt")
    
    if format == 'fastq':
        for title, seq, qual in tqdm(FastqGeneralIterator(fqhandle),
            position=n, leave=None, unit_scale=0.25, desc=fqfile,
            unit='reads', mininterval=1, disable=None):
            seqcount += 1
            if(seq.startswith(prefix)) :
                prefixcount += 1
                aseq = seq[len(prefix):]
                if(len(aseq) > 20) :
                    aseq = aseq[0:20]
                if aseq in aseqs:
                    aseqs[aseq] += 1
                else:
                    aseqs[aseq] = 1
    elif format == 'fasta':
        for record in tqdm(SeqIO.parse(fqhandle, "fasta"),
            position=n, leave=None, unit_scale=0.25, desc=fqfile,
            unit='reads', mininterval=1, disable=None
        ):
            seqcount += 1
            faseq = str(record.seq)
            if faseq.startswith(prefix):
                prefixcount += 1
                aseq = faseq[len(prefix):]
                if(len(aseq) > 20) :
                    aseq = aseq[0:20]
                if aseq in aseqs:
                    aseqs[aseq] += 1
                else:
                    aseqs[aseq] = 1

    fqhandle.close()
    #print("\t\tPrefix", prefix, "found", prefixcount, "times out of", seqcount,
    #    " records.")
    if prefixcount > 0:
        bestaseq = sorted(aseqs, key=aseqs.get, reverse=True)[0]
        bestperc = round((aseqs[bestaseq] / prefixcount) * 100, 1)
        #print("\t\tMost frequent adapter (first 20nts):", bestaseq, "found",
        #    aseqs[bestaseq], "times (", bestperc, '%)')
        ###print (f'For {fqfile} adapter {bestaseq} will be used for trimming.')
        return (fqfile, bestaseq) 
        
    else:
        #print(f'Autotrim failed to find an adapter for {fqfile}!')
        #print(f' Consider a different autotrim_key!') 
        # We don't sys.exit() here because this function may be running
        #  as one of several threads. Main script should abort instead.
        return (fqfile, None)

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

    with bgzf.BgzfWriter(os.path.join(UMI_PATH, BASENAME + "_umi.fa.gz"), "wb") as out:
        SeqIO.write(
            (SeqRecord(Seq(i_u.split("_")[0]), id="_".join([i_u, str(number)]), description="")
            for i_u, number in insert_umi.items()),
            out,
            "fasta"
        )
    
    bowtie_fa_cmd = [
    "bowtie",
    "-f",
    "-p", "24",
    "-m", "50",
    "-v", "0",
    "--best",
    "--strata",
    "-a",
    "--no-unal",
    "-x", genomeindex, 
    UMI_PATH + "/" + BASENAME + '_umi.fa.gz',
    "-S", UMI_PATH + '/' + BASENAME + '_aligned.sam',
    "--al", UMI_PATH + '/' + BASENAME + '_aligned.fa'
    ]
    bowtie_log = open(UMI_PATH + '/' + BASENAME + '.mapresults.txt', 'w')
    log = subprocess.run(bowtie_fa_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    bowtie_log.write(log.stderr)
    bowtie_log.close()
    
def data_mapping(BASENAME, TRIM_PATH, MAPPING_PATH):
    bowtie_fq_cmd = [
    "bowtie",
    "-p", "24",
    "-m", "50",
    "-v", "0",
    "--best",
    "--strata",
    "-a",
    "--no-unal",
    "-x", genomeindex, 
    TRIM_PATH + "/" + BASENAME + '_trimmed.fq.gz',
    "-S", MAPPING_PATH + '/' + BASENAME + '_aligned.sam',
    "--al", MAPPING_PATH + '/' + BASENAME + '_aligned.fq'
    ]
    bowtie_log = open(MAPPING_PATH + '/' + BASENAME + '.mapresults.txt', 'w')
    log = subprocess.run(bowtie_fq_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    bowtie_log.write(log.stderr)
    bowtie_log.close()

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


def umi_kit():
    if not os.path.exists(UMI_PATH):
        os.mkdir(UMI_PATH)
    TRIM_FLAG = 2
    ADAPTER = None
    data_trim(ADAPTER, TRIM_FLAG, UMI_PATH, fastq_file)
    data_umi_extract(BASENAME, UMI_PATH, UMI_PATH)
    data_length_dist(UMI_PATH + "/" + BASENAME + "_aligned.fa", DOC_PATH)

    TRIM_FLAG = 1
    ADAPTER = "AACTGTAGGCACCATCAAT"
    data_trim(ADAPTER, TRIM_FLAG, TRIM_PATH, fastq_file)
    data_mapping(BASENAME, TRIM_PATH, MAPPING_PATH)

def default_kit():
    TRIM_FLAG = 2
    ADAPTER = None
    data_trim(ADAPTER, TRIM_FLAG, TRIM_PATH, fastq_file)
    data_mapping(BASENAME, TRIM_PATH, MAPPING_PATH)

# Argument parsing / help message / version
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))
parser.add_argument("-v", "--version", action="version",
                    version='%(prog)s v2.0.20240507')
parser.add_argument("-i", "--inputdir", type= str, default= os.getcwd(),
                    help="path to input directory, default is current directory")
parser.add_argument("-o", "--outdir", type= str, default= os.getcwd(),
                    help="path to output directory, default is current directory")
parser.add_argument("-g", "--genome_index", type= str, default= "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_bowtie_index/Arabidopsis_thaliana.TAIR10.dna.toplevel",
                    help="path to genome index file")
parser.add_argument("-u", "--umi_flag", type= int, default= 1)
parser.add_argument("-s", "--suffix", type= str, default= "_R1.fastq.gz")


args = parser.parse_args()


if __name__ == "__main__":

    INPUTDIR = args.inputdir
    OUTDIR = args.outdir
    DOC_PATH = args.outdir
    UMI_FLAG = args.umi_flag
    
    suffix = args.suffix
    genomeindex = args.genome_index

    TRIM_PATH = os.path.join(OUTDIR, "2_trim_adapter/")
    UMI_PATH = os.path.join(OUTDIR, "3_extract_umi/")
    MAPPING_PATH = os.path.join(OUTDIR, "4_mapping/")
    BW_PATH = os.path.join(OUTDIR, "5_bam2bigwig/")

    fastq_files = glob.glob(INPUTDIR + "/*" + suffix, recursive=True)

    if not os.path.exists(TRIM_PATH):
        os.mkdir(TRIM_PATH)
    if not os.path.exists(MAPPING_PATH):
        os.mkdir(MAPPING_PATH)

    for fastq_file in fastq_files:
        BASENAME = os.path.basename(fastq_file).replace(suffix, "")
        print(BASENAME)
        if UMI_FLAG == 1:
            umi_kit()
        else:
            default_kit()
        convert_bam2bigwig(MAPPING_PATH, BW_PATH)
