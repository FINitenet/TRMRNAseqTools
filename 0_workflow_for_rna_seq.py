#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :0_workflow_for_rna_seq.py
# @Time      :2023/03/16 16:28:47
# @Author    :Yuchen@rlab

import os
import glob
import subprocess
import pysam
import argparse
import datetime
from multiprocessing import Process

def build_meta_info(meta_file):
    with open(meta_file, "r") as f:
        file_lines = f.readlines()
        meta = {}
        for line in file_lines:
            values = line.strip().split(" ")
            key = values[0]
            meta[key] = values[1:]
        return(meta)

def rawdata_qc(input_path, output_path):
    qc_path = os.path.join(output_path, "2_rawdata_qc/")
    if not os.path.exists(qc_path):
        os.mkdir(qc_path)
        fastq_files = os.path.join(input_path, "*.fastq.gz")
        # fastqc_cmd = ["fastqc", "-t", "24", fastq_files,"-o", qc_path]
        # subprocess.call(fastqc_cmd)
        fastqc_cmd = "fastqc -t 24 -o " + qc_path + " " + fastq_files
        os.system(fastqc_cmd)
    else:
        print("2_rawdata_qc already exists, skip this step.")


# %%
def data_trim(meta, input_path, output_path):
    trim_path = os.path.join(output_path, "3_trimmed_data/")
    if not os.path.exists(trim_path):
        os.mkdir(trim_path)
        for key, values in meta.items():
            library_type = values[2]
            if library_type == "paired":
                fastq1 = os.path.join(input_path, key + "_1.fastq.gz")
                fastq2 = os.path.join(input_path, key + "_2.fastq.gz")
                trim_cmd = ["trim_galore", "--paired", "--fastqc", "--gzip", "--suppress_warn", "-j", "24", "-o", trim_path, fastq1, fastq2]
                subprocess.run(trim_cmd, check=True)
            elif library_type == "single":
                fastq = os.path.join(input_path, key + ".fastq.gz")
                trim_cmd = ["trim_galore", "--fastqc", "--gzip", "--suppress_warn", "-j", "24", "-o", trim_path, fastq]
                subprocess.run(trim_cmd, check=True)
            else:
                print("Error: %s" % library_type)
    return(trim_path)

# %%
def data_mapping(meta, input_path, output_path):
    mapping_path = os.path.join(output_path, "4_mapping/")
    if not os.path.exists(mapping_path):
        os.mkdir(mapping_path)
        for key, values in meta.items():
            library_type = values[2]
            mapping_result = os.path.join(mapping_path, key + "_aligned.sam")
            mapping_summary = os.path.join(mapping_path, key + "_mapping_summary.txt")
            if library_type == "paired":
                fastq1 = os.path.join(input_path, key + "_val_1.fq.gz")
                fastq2 = os.path.join(input_path, key + "_val_2.fq.gz")
                mapping_cmd = [
                    "hisat2", "-p", "24", "-t",
                    "-x", bowtie2_index,
                    "-1", fastq1, "-2", fastq2,
                    "-S", mapping_result,
                    "--summary-file", mapping_summary,
                    "--new-summary", "--no-unal", "--dta"
                ]
                subprocess.run(mapping_cmd, check=True)
            elif library_type == "single":
                fastq = os.path.join(input_path, key + "_trimmed.fq.gz")
                mapping_cmd = [
                    "hisat2", "-p", "24", "-t",
                    "-x", bowtie2_index,
                    "-U", fastq,
                    "-S", mapping_result,
                    "--summary-file", mapping_summary,
                    "--new-summary", "--no-unal", "--dta"
                ]
                subprocess.run(mapping_cmd, check=True)
            else:
                print("Error: %s" % library_type)
    return(mapping_path)

# %%
def gene_annotation(meta, input_path, output_path):
    anno_path = os.path.join(output_path, "5_gene_annotation/")
    if not os.path.exists(anno_path):
        os.mkdir(anno_path)
        for key, values in meta.items():
            library_type = values[2]
            mapping_file = os.path.join(input_path, key + "_aligned.sam")
            
            if library_type == "paired":
                anno_cmd = [
                    "featureCounts",
                    "-a", gff_file,
                    "-o", anno_path + key + ".gene.counts.txt",
                    mapping_file,
                    "-t", "gene,ncRNA_gene",
                    "-g", "ID",
                    "-T", "48",
                    "-R", "BAM",
                    "-p",
                    "--largestOverlap"
                ]
                with open(anno_path + key + ".log", "wb") as log_file:
                    res = subprocess.run(anno_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    log_file.write(res.stdout) 
            elif library_type == "single":
                anno_cmd = [
                    "featureCounts",
                    "-a", gff_file,
                    "-o", anno_path + key + ".gene.counts",
                    mapping_file,
                    "-t", "gene,ncRNA_gene",
                    "-g", "ID",
                    "-T", "48",
                    "-R", "BAM",
                    "--largestOverlap"
                ]
                with open(anno_path + key + ".log", "wb") as log_file:
                    res = subprocess.run(anno_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    log_file.write(res.stdout)

def data_compress(sam_file, bam_file):
    pysam.sort("-o", bam_file, sam_file)
    pysam.index(bam_file)
    os.remove(sam_file)

def rna_preprocessing():
    meta = build_meta_info(meta_file)
    rawdata_qc(input_path, output_path)
    trim_path = data_trim(meta, input_path, args.outdir)
    mapping_path = data_mapping(meta, trim_path, args.outdir)
    gene_annotation(meta, mapping_path, args.outdir)
    

# Argument parsing / help message / version
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))
parser.add_argument("-v", "--version", action="version",
                    version='%(prog)s v1.0.20230321')
parser.add_argument("-i", "--inputdir", type= str, default= os.getcwd(),
                    help="path to input directory, default is current directory")
parser.add_argument("-o", "--outdir", type= str, default= os.getcwd(),
                    help="path to output directory, default is current directory")
parser.add_argument("-b", "--batch_file",
                    help="path to batch file, which contains sample information")

args = parser.parse_args()



if __name__ == "__main__":

    startTime = datetime.datetime.now()
    print('Start Time:', startTime)
    
    # global parameters
    input_path = args.inputdir
    output_path = args.outdir
    meta_file = args.batch_file
    bowtie2_index = "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_hisat_index/Arabidopsis_thaliana.TAIR10.dna.toplevel"
    gff_file = "/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/Arabidopsis_thaliana.TAIR10.53.chr.ridsiRNA.gff3"
    
    # main pipeline
    rna_preprocessing()
    
    # compress sam files to bam files
    sam_files = glob.glob(os.path.join(output_path + "/4_mapping/", "*.sam"))
    progress = []
    for sam_file in sam_files:
        bam_file = sam_file.replace(".sam", ".sorted.bam")
        p = Process(target=data_compress, args=(sam_file, bam_file))
        p.start()
        progress.append(p)
    for p in progress:
        p.join()

    print("All done!")
    endTime = datetime.datetime.now()
    time = (endTime - startTime).seconds
    print('End Time:', endTime)
    print("This programme run: %s s" % (time))