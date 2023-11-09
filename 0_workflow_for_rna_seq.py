#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :0_workflow_for_rna_seq.py
# @Time      :2023/10/08 23:22:50
# @Author    :Yuchen@rlab

import os
import glob
import subprocess
import pysam
import argparse
import datetime
import pandas as pd
from multiprocessing import Process


def build_meta_dict(input_path):
    meta_dict = {}
    # Define the allowed library layouts
    allowed_layouts = ['PAIRED', 'SINGLE']
    # Check if the library layout is valid
    if library_layouts not in allowed_layouts:
        print('ERROR: library_layouts is not PAIRED or SINGLE')
    else:
    # Loop over all files in the directory
        for filename in os.listdir(input_path):
        # Loop over all run accessions
            for run in run_accessions:
            # Check if the filename starts with the run accession
                if filename.startswith(run):
                # Check the library layout
                    if library_layouts == 'PAIRED':
                    # Check if the filename ends with _1.fastq.gz or _2.fastq.gz
                        if filename.endswith(('_1.fastq.gz', '_2.fastq.gz')):
                            filename = filename.replace('.fastq.gz', '')
                            meta_dict.setdefault(run, []).append(filename)
                    elif library_layouts == 'SINGLE':
                    # Check if the filename ends with .fastq.gz
                        if filename.endswith('.fastq.gz'):
                            filename = filename.replace('.fastq.gz', '')
                            meta_dict.setdefault(run, []).append(filename)
    return meta_dict
        
def data_trim(meta_dict, trim_path):
    if not os.path.exists(trim_path):
        os.mkdir(trim_path)
    for key, values in meta_dict.items():
        if library_layouts == 'PAIRED':
            fastq1 = os.path.join(input_path, values[0] + ".fastq.gz")
            fastq2 = os.path.join(input_path, values[1] + ".fastq.gz")
            trim_cmd = ["trim_galore", "--paired", "--fastqc", "--gzip", "--suppress_warn", "-j", "8", "-o", trim_path, fastq1, fastq2]
            subprocess.run(trim_cmd, check=True)
        elif library_layouts == 'SINGLE':
            fastq = os.path.join(input_path, values[0] + ".fastq.gz")
            trim_cmd = ["trim_galore", "--fastqc", "--gzip", "--suppress_warn", "-j", "8", "-o", trim_path, fastq]
            subprocess.run(trim_cmd, check=True)
        else:
            print("ERROR: library_layouts is not PAIRED or SINGLE")
    print("Trimming is done!")


def data_mapping(meta_dict, mapping_path):
    if not os.path.exists(mapping_path):
        os.mkdir(mapping_path)
    for key, values in meta_dict.items():
        mapping_result = os.path.join(mapping_path, key + "_aligned.sam")
        mapping_summary = os.path.join(mapping_path, key + "_mapping_summary.txt")
        if library_layouts == 'PAIRED':
            fastq1 = os.path.join(trim_path, values[0] + "_val_1.fq.gz")
            fastq2 = os.path.join(trim_path, values[1] + "_val_2.fq.gz")
            mapping_cmd = [
                "hisat2", "-p", "24", "-t",
                "-x", hisat_index,
                "-1", fastq1, "-2", fastq2,
                "-S", mapping_result,
                "--summary-file", mapping_summary,
                "--new-summary", "--no-unal", "--dta"
            ]
            subprocess.run(mapping_cmd, check=True)
        elif library_layouts == 'SINGLE':
            fastq = os.path.join(trim_path, values[0] + "_trimmed.fq.gz")
            mapping_cmd = [
                "hisat2", "-p", "24", "-t",
                "-x", hisat_index,
                "-U", fastq,
                "-S", mapping_result,
                "--summary-file", mapping_summary,
                "--new-summary", "--no-unal", "--dta"
            ]
            subprocess.run(mapping_cmd, check=True)
        else:
            print("ERROR: library_layouts is not PAIRED or SINGLE")
    return(mapping_path)

def gene_annotation(meta, anno_path):
    if not os.path.exists(anno_path):
        os.mkdir(anno_path)
    for key, values in meta.items():
        mapping_file = os.path.join(mapping_path, key + "_aligned.sam")
        if library_layouts == 'PAIRED':
            anno_cmd = [
                "featureCounts",
                "-a", gff_file,
                "-o", anno_path + key + ".gene.counts",
                mapping_file,
                "-t", "gene,ncRNA_gene",
                "-g", "ID",
                "--extraAttributes", "biotype",
                "-T", "48",
                "-p",
                "-O",
                "-M"
            ]
            with open(anno_path + key + ".log", "wb") as log_file:
                res = subprocess.run(anno_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                log_file.write(res.stdout) 
        elif library_layouts == 'SINGLE':
            anno_cmd = [
                "featureCounts",
                "-a", gff_file,
                "-o", anno_path + key + ".gene.counts",
                mapping_file,
                "-t", "gene,ncRNA_gene",
                "-g", "ID",
                "--extraAttributes", "biotype",
                "-T", "48",
                "-O",
                "-M"
            ]
            with open(anno_path + key + ".log", "wb") as log_file:
                res = subprocess.run(anno_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                log_file.write(res.stdout)

        df = pd.read_csv(anno_path + key + ".gene.counts", sep="\t", skiprows=[0])
        # 获取样本名列表
        sample_names = df.columns[1:]
        # 去掉路径，并更新样本名
        new_sample_names = [name.split('/')[-1].replace(".sam","") for name in sample_names]
        df.columns = ['Geneid'] + new_sample_names
        # 保存回文件
        df.to_csv(anno_path + key + ".gene.counts", sep='\t', index=False)


def data_compress(sam_file, bam_file):
    pysam.sort("-o", bam_file, sam_file)
    pysam.index(bam_file)
    os.remove(sam_file)

def merge_gene_matrix(anno_path):
    # 要查找的目录路径
    path = anno_path

    # 合并结果保存的文件路径
    output_file = anno_path + "merged_gene_matrix.txt"

    columns_to_merge = 7  # 前多少列需要合并，根据不同的定量方法进行修改，这里是featureCounts的结果，需要合并前6列

    # 定义用于保存数据的字典
    data_dict = {}

    # 遍历指定目录下的所有文件
    for file in os.listdir(path):
        # 如果文件名以 "gene.counts" 结尾
        if file.endswith("gene.counts"):
            # 打开文件并逐行读取
            with open(os.path.join(path, file), "r") as f:
                for line in f:
                    line_data = tuple(line.strip().split("\t")[:columns_to_merge])
                    # 将前 columns_to_merge 列作为键将行数据添加到 data_dict 中
                    if line_data not in data_dict:
                        data_dict[line_data] = line.strip().split("\t")[columns_to_merge:]
                    else:
                        data_dict[line_data].extend(line.strip().split("\t")[columns_to_merge:])
            
    # 将结果写入输出文件中
    with open(output_file, "w") as out_f:
        for key, values in data_dict.items():
            out_f.write("\t".join(list(key) + values) + "\n")

# Argument parsing / help message / version
parser = argparse.ArgumentParser(prog=os.path.basename(__file__))
parser.add_argument("-v", "--version", action="version",
                    version='%(prog)s v1.0.20231008')
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
    # input_path = '/home/chenyc/Bioinformatics/chenyc/DataBase/NCBI/DataBase_RNA/PRJNA893215'
    # output_path = '/home/chenyc/Bioinformatics/chenyc/test/test4rna'
    input_path = args.inputdir
    output_path = args.outdir

    hisat_index = '/bios-store1/chenyc/Reference_Source/Arabidopsis_Reference/ath_chr_hisat_index/Arabidopsis_thaliana.TAIR10.dna.toplevel'
    gff_file = '/home/chenyc/Bioinformatics/chenyc/Reference_Source/Arabidopsis_Reference/Arabidopsis_thaliana.TAIR10.53.chr.ridsiRNA.gff3'

    trim_path = os.path.join(output_path, "1_trimmed_data/")
    mapping_path = os.path.join(output_path, "2_map2genome/")
    anno_path = os.path.join(output_path, "3_gene_annotation/")

    # Read the SraRunTable.txt file as a DataFrame
    meta = pd.read_csv(args.batch_file)
    library_layouts = meta['LibraryLayout'].unique()
    run_accessions = meta['Run'].unique()
    
    # Step1: data trimming
    # build meta dict
    meta_dict = build_meta_dict(input_path)
    # data trimming
    data_trim(meta_dict, trim_path)
    # Step2: data mapping
    # data mapping
    mapping_path = data_mapping(meta_dict, mapping_path)
    # Step3: gene annotation
    # gene annotation
    gene_annotation(meta_dict, anno_path)
    # merge gene matrix
    merge_gene_matrix(anno_path)
    # Step4: compress sam files to bam files
    sam_files = glob.glob(os.path.join(output_path + "/2_map2genome/", "*.sam"))
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