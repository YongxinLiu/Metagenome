SHELL:=/bin/bash

	# 宏基因组第一版流程配置文件，请根据项目具体情况修改 v1.1 2019/4/23
	# Config file of Metagenome pipeline version 1, please modify according to experiment design
  # v1.0 2019/4/23 建议标准流程
  # v1.1 2020/12/28 开展测试完整流程

# 1. 有参分析流程参数 Parameters of reference-based pipeline

	# 工作目录 Working directory
	# 修改wd为当前工作目录pwd
	wd=`pwd`
	# 设置j最大运行任务/p线程数/p1非并行任务线程数，超过CPU数量效率反而会降低
	j=3
	p=24
	p1=72
	# make init # 建立分析所需子目录
	# 准备实验设计(result/design.txt)和测序数据(seq/*.fq.gz)和数据库(修改如下参数)


## 1.1. qc 质控和去宿主
	
	# 质控软件trimmomatic安装目录
	trimmomatic_path=/conda/share/trimmomatic/
	# 宿主基因组bowtie2索引，如人human/Homo_sapiens, 拟南芥ath/水稻rice/水麦wheat/苜蓿med/bt2
	host_bt2=/db/host/med/bt2
  # 接头文件：通常cleandata已经去除接头，因无须指定接头文件。检查multiqc中Adapter Content 如有接头，则查找接头文件/conda/share/trimmomatic/adapters/中，手动指定参数，如 ILLUMINACLIP:/conda/share/trimmomatic/adapters/TruSeq2/3-PE.fa:2:40:15

## 1.2. 物种和功能组成定量 humman2


## 1.3. 整理物种组成表和基本绘图 Summary metaphlan2 and plot
	
	# 绘制热图高丰度菌数据
	tax_top=25

## 1.4 kraken2
	
	kraken2_db=/db/kraken2
	kraken2_header=`tail -n+2 result/design.txt | cut -f 1|head -n1`



# 2. 无参分析流程 De novo assemble pipeline

## 2.1. khmer质控

	# 去除覆盖度>20区域，<3的低频kmer
	khmer_high=20
	khmer_low=2
	# 内存上限，最小8G，可根据内存调整，如1TB内存8个任务运行，可最大设置64G，保证内存使用不过半
	khmer_memory=64G

## 2.2. 组装 Assemble

	# megahit参数
	# 最小kmin，默认21，越小计算量越大，土壤推荐27，25-31范围合理
	kmin=27
	# 最长kmer，默认141适合PE150数据，PE100可设为91
	kmax=91
	# kmer步长，默认12，最小10低覆盖度组装更好，最大28，越大精度越差但计算更少
	kstep=12

	# 组装软件选择 Choose assemble method: megahit / metaspades
	assemble_method=megahit
	# 组装模式 single/group/all
	# 单样品组装，适合大数据megahit_single；按组最合理megahit_group；全部适合小样本megahit_all
	
	# 评估长阈值，默认500
	quast_len=200

	# salmon参数
	salmon_kmer=31

## 2.3. 基因组注释 Genome annotation

	# 基因注释软件选择 Choose assemble method: prokka / prodigal(基因存在多转录本，暂不考虑)
	annotation_method=prokka

### 2.3.3 构建非冗余基集 Non-redundancy gene set(大数据可选)

	cdhit_coverage=0.9
	cdhit_similarity=0.95
	cdhit_mem=900000


## 2.3. 物种注释 kraken2
	
	# 物种注释分三个级别：reads、contig、gene


## 2.4 功能数据库注释

### 2.4.1 eggNOG
	
	# 拆分行2000000，1M行
	# split_line=10000
	# eggnog数据库位置 483m
	eggnog_db=/home/meta/db/eggnog2
	# 复制数据库至内存22G，加速检索(前提内存足够大) cp /home/meta/db/eggnog2 /dev/shm
	# eggnog_db=/dev/shm

### 2.4.2 KEGG

	kegg_dmnd=/db/kegg/kegg_76
	# 标准化单位，默认为一百万tpm/rpm，可选1、100或1000
	unit=1000000

### 2.4.3 CAZyome 碳水化合物数据库

	# http://cys.bios.niu.edu/dbCAN2/download/
	dbcan2_dmnd=/db/protein/dbcan2/CAZyDB.07312018
	dbcan2_anno=/db/protein/dbcan2/fam_description.txt

### 2.4.4 CAZyome 碳水化合物数据库

	# http://www.dantaslab.org/resfams
	resfams_dmnd=/mnt/zhou/yongxin/db/ResFams/Resfams-proteins
	resfams_anno=/mnt/bai/yongxin/data/db/ResFams/Resfams-proteins_class.tsv

include /home/meta/soft/Metagenome/denovo1/pipeline.md
