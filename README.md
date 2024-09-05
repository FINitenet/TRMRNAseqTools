# TRMRNAseqTools

Integrated TRM sRNA Sequencing Data Analysis for Plant

Author: Yuchen Cai

Current Version: 1.0_Beta

Latest updata: 11/04/2024

## Introduction

The TRMRNASeqTools is a SHELL and Python based pipeline designed for automatic general analysis for TRM sRNA sequencing data in supported plants.

Currently it is able to process small RNA-seq and generally analyze them.

If you have any questions or comments, please submit an issue in the GitHub or directly email to Yuchen Cai.

## Function
![image](https://github.com/user-attachments/assets/2dc7966c-12fa-4ef9-bf41-1ae2863acb7c)

## Upate Summary

### V2.1.0(2022-1-10)

--在Wang et.al, 2019, PNAS对miRNA/miRNA重新梳理（*Reference_Source/Arabidopsis_Reference/ath_mature_bowtie_index/PNAS_miRNA.fa*）的基础上，修正了部分错误的序列信息，构建了新的miRNA/miRNA*序列文件 （*Reference_Source/Arabidopsis_Reference/ath_mature_bowtie_index/PNAS_miRNA_ver.211213.fa*）

--基于重新构建的miRNA/miRNA*序列文件，通过回帖到基因组，以及参照miRBase（*Reference_Source/Arabidopsis_Reference/ath.gff3*）hairpin的注释，构建了新的miRNA annotation GFF3文件（*Reference_Source/Arabidopsis_Reference/ath.reannotation.gff3*）

--在对reads计数时，除对miRBase mature 以及hairpin的计数结果外，增加miRNA re-annotation的计数结果，并以miRNA family的形式重新求和计算。

### V1.0.0.20220927_Beta(2022-9-27)

--使用更加规范的版本号

--添加--aligned-only参数

### V1.0.0.20221111_Beta(2022-11-11)

--使用python scripts 统计 genome RNA length distribution（*size_dist.py*)

### V1.0.0.20221120_Beta(2022-11-20)

--使用python scripts 统计各个长度上RNA的种类分布 (*calc_feature_count_by_featureCounts_tags.py*)

--修改了一些文字描述

### V1.0.0.20230313_Beta(2023-3-13)

--完全使用python scripts 统计各个长度上RNA的种类分布 (*calc_feature_count_by_featureCounts_tags.py*)，替换原代码 (*0_workflow_Annotation-type_len_dist.sh*)

    1.与原代码相同的统计方法，重复统计；

    2.multiple-annotation reads count平均分(1/len(feature));

    3.按照RNA type优先级对multiple-annotation reads进行注释。

--更换注释文件 (*TRMRNAseqTools/reference/Arabidopsis_thaliana.TAIR10.53.md.gff3*)
