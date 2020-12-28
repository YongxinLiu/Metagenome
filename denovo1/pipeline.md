#SHELL:=/bin/bash

	# 宏基因组分析流程第一版 
	# Metagenome pipeline version 1



help:

	# 帮助文档: 流程所需的软件、脚本及数据库版本
	# Help: Pipeline dependency version of softwares, scripts and databases
	
	# 软件 Softwares
	fastqc -v # 0.11.5 测序数据质量评估
	multiqc --version # 1.4
		multiqc -h # 质控合并
	kneaddata --version # v0.6.1 质控和去宿主
		kneaddata_read_count_table -h
	parallel --version # 并行/多线程任务管理 20141022
	humann2 --version # v0.11.1
		humann2_join_tables
		humann2_renorm_table
		humann2_split_stratified_table
	metaphlan2.py --version # 2.7.6 (5 March 2018)
		merge_metaphlan_tables.py -h
		metaphlan_hclust_heatmap.py -h
	# microbiome_helper # 2017, https://github.com/mlangill/microbiome_helper
		# metaphlan_to_stamp.pl
	graphlan.py -h # 0.9.7 (21 July 2014)
		export2graphlan.py
		graphlan_annotate.py
	run_lefse.py -h # 1.0
		lefse-format_input.py
		lefse-plot_cladogram.py
		lefse-plot_res.py
		lefse-plot_features.py
	kraken2 --version # 2.0.7-beta
	khmer # kmer质控软件
		interleave-reads.py # 双端序列交叉合并为单文件
		trim-low-abund.py # 低丰度kmer过滤
		split-paired-reads.py # 双端单文件还原为双端和非配对
	megahit -h # v1.1.3
	quast.py # 5.0.0
	salmon -v # 非比对基因定量软件 0.11.2
	prokka --version # 基因组注释 1.13.3
	cd-hit # 去冗余 
		cd-hit-est # 核酸水平去冗余
	diamond help # 类blast比对工具 v0.8.22.84
	hmmsearch -h # 结构域界定 3.1b2
		hmmbuild -h # hmm构建

	fastp -v # version 0.12.5 fastq文件质控、质量类型转换
	usearch10 --help # v10.0.240_i86linux64 扩增子分析软件
	clustalo --version # 1.2.1 多序列对齐
	biom --version # 2.1.5， OTU表格式转换

	# 脚本 Scripts
	#alpha_boxplot.sh -h # 1.0 基于usearch alpha_div绘制箱线图

	# 数据库 Databases

	# 更新日志 Update log
	# 2018-09-12 Add kneaddata, humann2 for reference based pipeline
	# 2018-09-23 Add kraken2 for taxonomy in reads levels, salmon genes quaitity
	# 2019-06-06 Update



# 1. 有参分析流程 Reference-based pipeline

	# 准备实验设计、测序数据和参数数据库
	# Prepare experiment design (result/metadata.txt), sequecing data (seq/*.fq.gz) and database

10init:

	# 清理零字节文件重启项目 Clean zero-byte files, restart project
	find . -name "*" -type f -size 0c | xargs -n 1 rm -f
	# 建立程序必须目录 Create basic directory
	mkdir -p seq temp result
	touch $@


## 1.1. 质控并移除宿主 Quality control & Remove host

### 1.1.1 质量评估 Quality access
	
	# 此步选做，一般原始序列合并时已经统计过数据量

11qa: 10init

	touch $@
	# 质量评估 Quality access
	time fastqc -t ${p1} seq/*1.fq.gz
	time fastqc -t ${p1} seq/*2.fq.gz
	# 质量报告汇总
	multiqc -d seq/ -o result/

### 1.1.2 质控并移除宿主 Quality control & Remove host

11qc: 10init

	touch $@
	mkdir -p temp/11qc
	# 并行质控去宿主 kneaddata call trimmomatic quality control and bowtie2 remove host
	time parallel --xapply -j ${j} \
		"kneaddata -i seq/{1}_1.fq.gz -i seq/{1}_2.fq.gz \
		-o temp/11qc -v -t ${p} --remove-intermediate-output \
		--trimmomatic ${trimmomatic_path} --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:100' \
		--bowtie2-options '--very-sensitive --dovetail' -db ${host_bt2}" \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
	kneaddata_read_count_table --input temp/11qc --output temp/11kneaddata_stat.txt
	cut -f 1,2,4,12 temp/11kneaddata_stat.txt | awk 'BEGIN{OFS=FS="\t"} {print $$0,$$3/$$2*100,$$4/$$3*100}' \
		| sed 's/_1_kneaddata//' | sed '1 s/-nan/Hi-Q%/;s/-nan/rm_host%/' > result/11kneaddata_stat.txt
	cat result/11kneaddata_stat.txt

### 1.13 提交数据准备 Submit clean data

	# 保存质控和去宿主后数据至seq/clean目录，连同保存宿主比例

11gsa: 11qc

	touch $@
	mkdir -p submit
	# hard link to clean
	ln temp/11qc/*data_paired* submit/
	# rename accroding to ID
	rename 's/_1_kneaddata_paired//;s/fastq/fq/' submit/*
	# compress for reducing space
	pigz submit/*
	ln result/11kneaddata_stat.txt submit/stat.txt
	ln result/metadata.txt submit/metadata.txt
	md5sum submit/*_1.fq.gz > submit/md5sum1.txt
	md5sum submit/*_2.fq.gz > submit/md5sum2.txt
	paste submit/md5sum1.txt submit/md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' > submit/md5sum.txt
	sed -i 's/submit\///g' submit/md5sum.txt
	cat submit/md5sum.txt


## 1.2. 物种和功能组成定量 humann2

### 1.2.1 humann2输入文件准备：双端文件cat合并

12humann2_concat: 11qc

	touch $@
	# 生成humann2输入要求的合并文件 cat pair-end for humann2
	mkdir -p temp/12concat
	parallel --xapply -j ${j} \
		"cat temp/11qc/{1}*data_paired* > temp/12concat/{1}.fq" \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
	# 查看样品数量和大小 show samples size
	ls -l temp/12concat/*.fq

### 1.2.2 humann2计算，包括metaphlan2

12humann2: 12humann2_concat

	touch $@
	mkdir -p temp/12humann2
	time parallel -j ${j} \
		'humann2 --input {} --threads ${p} \
		--output temp/12humann2/ ' \
		::: temp/12concat/*.fq 

### 1.2.3 功能组成整理 humann2_sum

12humann2_sum: 12humann2

	touch $@
	mkdir -p result/12humann2
	# 合并所有样品，通路丰度`pathabundance`包括各功能和具体的物种组成，还有基因家族`genefamilies`(太多)，通路覆盖度`pathcoverage`层面可以进分析
	humann2_join_tables --input temp/12humann2/ --file_name pathabundance --output result/12humann2/pathabundance.tsv
	sed -i 's/_Abundance//g' result/12humann2/pathabundance.tsv
	# 标准化为相对丰度relab或百万分数cpm
	humann2_renorm_table --input result/12humann2/pathabundance.tsv --units relab \
		--output result/12humann2/pathabundance_relab.tsv
	# 分层结果，结果stratified(每个菌的功能组成)和unstratified(功能组成)两个
	humann2_split_stratified_table --input result/12humann2/pathabundance_relab.tsv \
		--output result/12humann2/
	sed -i 's/# Pathway/MetaCyc_pathway/' result/12humann2/pathabundance_relab_*stratified.tsv


## 1.3. 整理物种组成表和基本绘图 Summary metaphlan2 and plot

### 1.3.1 整理物种组成表 Summary metaphlan2

13metaphaln2_sum: 12humann2_sum

	touch $@
	# metaphlan2功能组成
	mkdir -p result/13metaphlan2
	# 合并单样品为表
	merge_metaphlan_tables.py temp/12humann2/*_humann2_temp/*_metaphlan_bugs_list.tsv | \
		sed 's/_metaphlan_bugs_list//g' | sed '/^#/d' > result/13metaphlan2/taxonomy.tsv
	# 转换为stamp的多级格式，株水不完整，不去掉列无法对齐
	metaphlan_to_stamp.pl result/13metaphlan2/taxonomy.tsv > result/13metaphlan2/taxonomy.spf
	# 绘制热图
	metaphlan_hclust_heatmap.py --in result/13metaphlan2/taxonomy.tsv \
		--out result/13metaphlan2/taxonomy_heatmap_top.pdf \
		-c bbcry --top ${tax_top} --minv 0.1 -s log 

### 1.3.2 GraPhlAn图

13metaphaln2_graphlan: 13metaphaln2_sum
	touch $@
	# metaphlan2 to graphlan
	export2graphlan.py --skip_rows 1,2 -i result/13metaphlan2/taxonomy.tsv \
		--tree temp/13ref_taxonomy.tree --annotation temp/13ref_taxonomy.annot \
		--most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 \
		--annotations 5,6 --external_annotations 7 --min_clade_size 1
	# graphlan annotation
	graphlan_annotate.py --annot temp/13ref_taxonomy.annot temp/13ref_taxonomy.tree temp/13ref_taxonomy.xml
	# output PDF figure, annoat and legend
	graphlan.py temp/13ref_taxonomy.xml result/13metaphlan2/taxonomy_graphlan.pdf --external_legends --dpi 300 

### 1.3.3 物种组成LEfSe差异分析(可选)

	# 需要有完整的实验分组，且样本末尾数字为样本重复编号才可分析

13metaphaln2_lefse: 13metaphaln2_graphlan
	touch $@
	mkdir -p result/13metaphlan2_lefse
	# LEfSe差异分析和Cladogram
	# 修改样本品为组名，要求重复为结尾数值
	sed '1 s/[0-9]$$//g' result/13metaphlan2/taxonomy.tsv | grep -v '#' > result/13metaphlan2_lefse/lefse.txt
	# 格式转换为lefse内部格式
	lefse-format_input.py  result/13metaphlan2_lefse/lefse.txt temp/13lefse_input.in -c 1 -o 1000000
	# 运行lefse
	run_lefse.py temp/13lefse_input.in temp/13lefse_input.res
	# 绘制物种树注释差异
	lefse-plot_cladogram.py temp/13lefse_input.res result/13metaphlan2_lefse/lefse_cladogram.pdf --format pdf --dpi 600 
	# 绘制所有差异features柱状图
	lefse-plot_res.py temp/13lefse_input.res result/13metaphlan2_lefse/lefse_res.pdf --format pdf --dpi 600
	# 批量绘制所有差异features柱状图
	lefse-plot_features.py -f diff --archive none --format pdf \
		temp/13lefse_input.in temp/13lefse_input.res result/13metaphlan2_lefse/t_


## 1.4. kraken2物种组成(可选)

### 1.4.1 基于NCBI完整基因组数据库的k-mer物种注释 Taxonomy assign by k-mer and based on NCBI database

14kraken2_reads: 11qc
	touch $@
	mkdir -p temp/14kraken2_reads
	# Based on submited clean data
	time parallel -j ${j} \
		'kraken2 --db ${kraken2_db} --paired submit/{1}*.fq.gz \
		--threads ${p} --use-names --use-mpa-style --report-zero-counts \
		--report temp/14kraken2_reads/{1}_report \
		--output temp/14kraken2_reads/{1}_output' >> temp/14kraken2_reads/kraken2.log 2>&1 \
		::: `tail -n+2 result/metadata.txt | cut -f 1`

### 1.4.2 合并为矩阵 merge into matrix

14kraken2_reads_sum: 14kraken2_reads
	touch $@
	mkdir -p result/14kraken2_reads
	parallel -j ${j} \
		'sort temp/14kraken2_reads/{1}_report | cut -f 2 | sed "1 s/^/{1}\n/" > temp/14kraken2_reads/{1}_count ' \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
	sort temp/14kraken2_reads/${kraken2_header}_report | cut -f 1 | sed "1 s/^/Taxonomy\n/" > temp/14kraken2_reads/0header_count
	paste temp/14kraken2_reads/*count > result/14kraken2_reads/taxonomy_count.txt

## 1.5 Clean 清理有参过程临时文件

15clean_ref: 
	touch $@
	# 1.1 质控qc质控后归档为submit，删除qc临时文件
	rm -rf temp/11qc
	# 1.2 humman2筛选文件
	rm -rf temp/12concat
	# 备份临时文件中的metaphlan2重要结果
	mkdir temp/12metaphlan2
	cp -f temp/12humann2/*_humann2_temp/*_metaphlan_bugs_list.tsv temp/12metaphlan2/
	# 删除humann2筛选文件
	rm -rf temp/12humann2/*temp
	# 1.4 kraken2比对文件
	rm -r temp/14kraken2_reads/*output


# 2. 无参分析流程 De novo assemble pipeline

## 2.1. khmer质控(可选)

21khmer: 11qc
	touch $@
	mkdir -p temp/21khmer
	# -V is metagenome, not single genome. -M set max memory for save server. -f force write
	time parallel -j ${j} \
		'interleave-reads.py temp/11qc/{1}_1_kneaddata_paired_1.fastq temp/11qc/{1}_1_kneaddata_paired_2.fastq | \
		trim-low-abund.py -V -M ${khmer_memory} -Z ${khmer_high} -C ${khmer_low} - -o temp/21khmer/{1}.fq; \
		split-paired-reads.py -f -0 temp/21khmer/{1}_0.fq -1 temp/21khmer/{1}_1.fq -2 temp/21khmer/{1}_2.fq temp/21khmer/{1}.fq' \
		::: `tail -n+2 result/metadata.txt | cut -f 1`


## 2.2. Assemble 组装

### 2.2.0 基于qc质控序列单样本组装(大项目)

	# ifseq, else, endif必须顶格写
22assemble_single: 11gsa
	touch $@
	echo -e "Assemble method: ${assemble_method}\nAssemble mode: assemble single\n"
ifeq (${assemble_method}, megahit)
	# 采用metahit单样品拼接
	# rm -rf temp/22megahit
	mkdir -p temp/22megahit
	time parallel -j ${j} \
		'megahit -t ${p1} --k-min ${kmin} --k-max ${kmax} --k-step ${kstep} \
		-1 submit/{1}_1.fq.gz -2 submit/{1}_2.fq.gz \
		-o temp/22megahit/{1} >> temp/22megahit/megahit.log 2>&1 ' \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
else ifeq (${assemble_method}, metaspades)
	# metaspades单样品拼接
	parallel --xapply -j ${j} \
		"metaspades.py -t ${p} -m 500 \
		-1 temp/11qc/{1}_1_kneaddata_paired_1.fastq \
		-2 temp/11qc/{1}_1_kneaddata_paired_2.fastq \
		-o temp/22metaspades/{1}" \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in megahit or metaspades")
endif
	echo -ne 'Assemble finished!!!\n' 

### 2.2.1 基于khmer质控后序列拼接(可选)

22megahit_all_k: 21khmer
	touch $@
	# 采用metahit拼接所有样本
	# rm -rf temp/22megahit_all
	time megahit -t ${p1} --k-min ${kmin} --k-max ${kmax} --k-step ${kstep} \
	-1 `ls temp/21khmer/*_1.fq|tr '\n' ','|sed 's/,$$//'` \
	-2 `ls temp/21khmer/*_1.fq|tr '\n' ','|sed 's/,$$//'` \
	-o temp/22megahit_all

### 2.2.2 基于qc质控后序列混合组装(6-30个样品小项目)

22megahit_all: 11qc
	touch $@
	# 采用metahit拼接所有样本，理论上双端序列文件大小完全相关，为何有差异
	# rm -rf temp/22megahit_all
	# --k-min ${kmin} --k-max ${kmax} --k-step ${kstep} 
	time megahit -t ${p1} \
		-1 `ls temp/11qc/*_paired_1.fastq|tr '\n' ','|sed 's/,$$//'` \
		-2 `ls temp/11qc/*_paired_2.fastq|tr '\n' ','|sed 's/,$$//'` \
		-o temp/22megahit_all --continue

## 2.2.3 megahit_all_quast评估

	# 依赖 22megahit_all | 22megahit_all_k
22megahit_all_quast: 22megahit_all
	touch $@
	# 拼接结果评估
	mkdir -p result/22megahit_all
	ln -f temp/22megahit_all/final.contigs.fa result/22megahit_all/
	time quast.py -o result/22megahit_all/ -m ${quast_len} -t ${p} result/22megahit_all/final.contigs.fa

### 2.2.4 Contig定量salmon

22megahit_all_salmon: 22megahit_all_quast
	touch $@
	mkdir -p temp/22salmon_contig
	# 1. 建索引
	# -t 转录本序列，--type 类型fmd/quasi，-k kmer长度默认31, -i 索引
	salmon index -t temp/22megahit_all/final.contigs.fa -p ${p} \
		-i temp/22salmon_contig/index --type quasi -k ${salmon_kmer}
	# 2. 定量
	parallel -j ${j} \
		'salmon quant -i temp/22salmon_contig/index -l A -p ${p} --meta \
		-1 temp/11qc/{1}_1_kneaddata_paired_1.fastq -2 temp/11qc/{1}_1_kneaddata_paired_2.fastq \
		-o temp/22salmon_contig/{1}.quant' \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
	# 3. 合并
	salmon quantmerge --quants temp/22salmon_contig/*.quant -o result/22megahit_all/contig.TPM
	salmon quantmerge --quants temp/22salmon_contig/*.quant --column NumReads -o result/22megahit_all/contig.count
	sed -i '1 s/.quant//g' result/22megahit_all/contig.*

### 2.2.5 Contig物种注释 kraken2

22kraken2_contig: 22megahit_all_salmon
	touch $@
	kraken2 --db ${kraken2_db} temp/22megahit_all/final.contigs.fa \
	--threads ${p} --use-names \
	--output result/22megahit_all/kraken2.tax
	# report 为 contig.count的物种注释


## 2.3. Gene annotation 基因注释


### 2.3.2 对合并组装的单个contig文件基因注释

23prodigal_all: 22megahit_all
	touch $@
	mkdir -p temp/23prodigal_all
	time prodigal -i temp/22megahit_all/final.contigs.fa -d temp/23prodigal_all/gene.fa  \
	  -o temp/23prodigal_all/gene.gff -p meta -f gff \
	  > temp/23prodigal_all/gene.log 2>&1 

23prodigal_all_sum: 23prodigal_all
	touch $@
	ls -lsh temp/23prodigal_all/gene.fa
	echo -ne 'metaProdigal all genes\t' > result/gene.log
	grep -c '>' temp/23prodigal_all/gene.fa>> result/gene.log
	cat result/gene.log
	# link to NR
	mkdir -p temp/23NRgene
	rm -rf temp/23NRgene/gene.fa
	ln temp/23prodigal_all/gene.fa temp/23NRgene/gene.fa

### 2.3.1 对单样品组装的每个contig文件基因注释(大数据可选)

	# ifseq, else, endif必须顶格写
23prokka_single: 22assemble_single
	touch $@
	echo -e "Assemble method: ${assemble_method}\nAssemble mode: assemble single\n"
ifeq (${assemble_method}, megahit)
	# 采用metahit单样品拼接
	time parallel --xapply -j ${j} \
		"prokka temp/22megahit/{1}/final.contigs.fa --outdir temp/23prokka/{1} \
		-prefix mg --metagenome --force --cpus ${p} \
		--kingdom Archaea,Bacteria,Mitochondria,Viruses" \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
	# DNA level
	cat temp/23prokka/*/mg.ffn > temp/23prokka/mg.ffn
else ifeq (${assemble_method}, metaspades)
	# metaspades单样品拼接
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in megahit or metaspades")
endif
	# 复制结果到新目录，单/多样本分析汇总
	mkdir -p temp/23NRgene
	cp temp/23prokka/mg.ffn temp/23NRgene/
	echo -ne 'Prokka single genes\t' > result/gene.log
	grep -c '>' temp/23NRgene/mg.ffn >> result/gene.log

### 2.3.2 对合并组装的单个contig文件基因注释

23prokka_all: 22megahit_all_quast
	touch $@
	time prokka temp/22megahit_all/final.contigs.fa --outdir temp/23prokka_all \
	-prefix mg --metagenome --force --cpus ${p1} \
	--kingdom Archaea,Bacteria,Mitochondria,Viruses
	mkdir -p temp/23NRgene
	ln temp/23prokka_all/mg.ffn temp/23NRgene/

23prokka_all_check: 
	# prokka默认步骤较多，但我们需要的基因注释可能已经完成
	# 检查文件时间，如小于当前已经完成，可手动转换并继续
	ll temp/23prokka_all/mg.ffn
	mkdir -p temp/23NRgene
	cp temp/23prokka_all/mg.ffn temp/23NRgene/
	echo -ne 'Prokka all genes\t' > result/gene.log
	grep -c '>' temp/23NRgene/mg.ffn >> result/gene.log
	cat result/gene.log


### 2.3.3 构建非冗余基集 Non-redundancy gene set(大数据可选)

23NRgeneSet: 
	touch $@
	# -c sequence identity threshold, default 0.9; -M  max available memory (Mbyte), default 400; -l  length of throw_away_sequences, default 10
	time cd-hit-est -i temp/23NRgene/gene.fa -o temp/23NRgene/NRgene.fa -aS ${cdhit_coverage} -c ${cdhit_similarity} -G 0 -g 0 -M ${cdhit_mem} -T ${p1} # -n 10 -d 0 -g 1
	# 统计基因数据，确定ID是否非冗余
	echo -ne 'NRgeneSet\t' >> result/gene.log
	grep -c '>' temp/23NRgene/NRgene.fa >> result/gene.log
	cat result/gene.log
	# 翻译核酸为对应蛋白序列， emboss
	transeq -sequence temp/23NRgene/NRgene.fa -outseq temp/23NRgene/protein.fa -trim Y 
	# 序列名自动添加了_1，为与核酸对应要去除
	sed -i 's/_1 / /' temp/23NRgene/protein.fa

### 2.3.4 基因定量 salmon genes

23salmon_gene: 23NRgeneSet
	touch $@
	mkdir -p temp/23salmon_gene
	# 1. 建索引
	# -t 转录本序列，--type 类型fmd/quasi，-k kmer长度默认31, -i 索引
	salmon index -t temp/23NRgene/NRgene.fa -p ${p} \
		-i temp/23salmon_gene/index # --type quasi -k ${salmon_kmer}
	# 2. 定量
	parallel -j ${j} \
		'salmon quant -i temp/23salmon_gene/index -l A -p ${p} --meta \
		-1 submit/{1}_1.fq.gz -2 submit/{1}_2.fq.gz \
		-o temp/23salmon_gene/{1}.quant' \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
	# 3. 合并
	mkdir -p result/23salmon_gene
	salmon quantmerge --quants temp/23salmon_gene/*.quant -o result/23salmon_gene/gene.TPM
	salmon quantmerge --quants temp/23salmon_gene/*.quant --column NumReads -o result/23salmon_gene/gene.count
	sed -i '1 s/.quant//g' result/23salmon_gene/gene.*

### 2.3.5 基因物种注释 kraken2 annotate gene(可选)

23kraken2_gene: 23salmon_gene
	touch $@
	kraken2 --db ${kraken2_db} temp/23NRgene/mg.ffn.nr \
	--threads ${p} --use-names \
	--output result/23salmon_gene/kraken2.tax


## 2.4 功能数据库注释

### 2.4.1 eggNOG

24eggnog: 23NRgeneSet
	touch $@
	mkdir -p temp/24eggnog
	# diamond比对基因至eggNOG数据库
	time emapper.py -m diamond --no_annot --no_file_comments \
	  --data_dir ${eggnog_db} --cpu ${p1} -i temp/23NRgene/protein.fa \
	  -o temp/24eggnog/protein --override
	# 比对结果功能注释
	time emapper.py --annotate_hits_table \
	  temp/24eggnog/protein.emapper.seed_orthologs --no_file_comments \
	  -o temp/24eggnog/output --cpu ${p1} --data_dir ${eggnog_db} --override


24eggnog_sum: 24eggnog

	touch $@
# 	mkdir -p result/24eggnog
# 	# 添加表头 Add header
# 	sed '1 i Name\tortholog\tevalue\tscore\ttaxonomic\tprotein\tGO\tEC\tKO\tPathway\tModule\tReaction\trclass\tBRITE\tTC\tCAZy\tBiGG\ttax_scope\tOG\tbestOG\tCOG\tdescription' temp/24eggnog/output.emapper.annotations > temp/24eggnog/output
#   # 合并
#   /conda/envs/humann3/bin/python /mnt/m1/yongxin/bin/summarizeAbundance.py \
#     -i result/salmon/gene.TPM \
#     -m temp/eggnog/output \
#     -c '9,16,21' -s ',+,+*' -n raw \
#     -o result/eggnog/eggnog
#   # eggnog.CAZy.raw.txt  eggnog.COG.raw.txt  eggnog.KO.raw.txt
#   # STAMP的spf格式，结合metadata.tsv进行COG/KO/Description差异比较
#   # KO
#   sed -i 's/^ko://' result/eggnog/eggnog.KO.raw.txt
#   awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
#     ${db}/eggnog/KO.anno \
#     result/eggnog/eggnog.KO.raw.txt| \
#     sed 's/^\t/Description\t/' > result/eggnog/eggnog.KO.TPM.spf
#   # CAZy
#   awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
#      ${db}/dbcan2/fam_description.txt result/eggnog/eggnog.CAZy.raw.txt | \
#     sed 's/^\t/Description\t/' > result/eggnog/eggnog.CAZy.TPM.spf
#   # COG
#   awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
#     ${db}/eggnog/COG.anno \
#     result/eggnog/eggnog.COG.raw.txt | sed '1 s/^/Level1\tLevel2/'> \
#     result/eggnog/eggnog.COG.TPM.spf

	# # 重点1序列名，7KO，12COG分类，13注释

	# # 1. KO注释表
	# # 提取基因KO表，基因1对多个KO时只提取第一个KO
	# cut -f 1,7 temp/24eggnog/output|cut -f 1 -d ','|grep -v -P '\t$$' > temp/24eggnog/1ko.list
	# # 基因丰度矩阵末尾添加对应KO编号，没注释的直接删除，可选注释为unclassified
	# awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print $$0,a[$$1]}' temp/24eggnog/1ko.list result/23salmon_gene/gene.count | \
	# 	sed '/\t$$/d' > temp/24eggnog/gene_ko.count
	# # 检查注释前后基因数量，2.6M，1M
	# wc -l result/23salmon_gene/gene.count 
	# # KO可注释基因数量
	# wc -l temp/24eggnog/gene_ko.count
	# # 合并基因表为KO表，输出count值和tpm值
	# time Rscript /home/meta/soft/Metagenome/denovo1/script/mat_gene2ko.R -i temp/24eggnog/gene_ko.count -o result/24eggnog/kotab -n ${unit}
	# # KO的数量
	# wc -l result/24eggnog/kotab.count
	# # 添加来自KEGG整理的KO对应的描述
	# # awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print $$0,a[$$1]}' /db/eggnog/KO.anno result/24eggnog/kotab.count > result/24eggnog/kotab.count.anno
	# # STAMP的spf格式，结果metadata.txt进行KO或Description差异比较
	# awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print a[$$1],$$0}' /db/eggnog/KO.anno result/24eggnog/kotab.count | sed 's/^\t/Undescription\t/' > result/24eggnog/kotab.count.spf
	# 
	# # 2. COG注释表
	# # 提取12列COG分类注释
	# cut -f 1,12 temp/24eggnog/output_file|cut -f 1 -d ','|grep -v -P '\t$$' > temp/24eggnog/1cog.list
	# # 检查注释前后基因数量，2.6M基因，1.6M基因有cog注释，远高于KEGG的1M
	# wc -l temp/24eggnog/1cog.list
	# # 基因丰度矩阵末尾添加对应cog编号，没注释的直接删除，可选注释为unclassified
	# awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print $$0,a[$$1]}' temp/24eggnog/1cog.list result/23salmon_gene/gene.count | \
	# 	sed '/\t$$/d' > temp/24eggnog/gene_cog.count
	# # 合并基因表为cog表，输出count值和tpm值
	# sed -i '1 s/\tCOG/\tKO/' temp/24eggnog/gene_cog.count
	# # 文件1.6M行，36列，按组合并2min
	# cat temp/24eggnog/gene_cog.count | datamash check
	# time Rscript /home/meta/soft/Metagenome/denovo1/script/mat_gene2ko.R -i temp/24eggnog/gene_cog.count -o result/24eggnog/cogtab -n ${unit}
	# # cog对应的描述，STAMP的spf格式，结合metadata.txt进行cog或Description差异比较
	# awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2"\t"$$3} NR>FNR{print a[$$1],$$0}' /db/eggnog/COG.anno result/24eggnog/cogtab.count | sed 's/^\t/Undescription\t/' > result/24eggnog/cogtab.count.spf
	# # 添加整理柱状图，和分组柱状图绘制
	# 
	# # 3. 基因功能描述
	# # 提取基因anno分类表，基因1对多个anno时只提取第一个anno
	# cut -f 1,13 temp/24eggnog/output_file|cut -f 1 -d ','|grep -v -P '\t$$' > temp/24eggnog/3anno.list
	# # 1.6M的功能注释，远高于KEGG的1M
	# wc -l temp/24eggnog/3anno.list
	# # anno对应的描述，STAMP的spf格式，没注释的基因删除，结合metadata.txt进行anno或Description差异比较
	# # 可选sed将末注释归类为Undescription sed 's/^\t/Undescription\t/'
	# awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print a[$$1],$$0}' temp/24eggnog/3anno.list result/23salmon_gene/gene.count | grep -v -P '^\t' > result/24eggnog/annotab.count.spf


### 2.4.2 KEGG(数据库无版权，请跳过)

24kegg: 23NRgeneSet
	touch $@
	mkdir -p temp/24kegg
	# 2.6 M genes, 8 threads, spend 6h
	time diamond blastp --db ${kegg_dmnd} --query temp/23NRgene/protein.fa \
		--outfmt 6 --threads ${p} --max-target-seqs 1 --quiet \
		--out temp/24kegg/gene_diamond.f6

24kegg_sum: 24kegg
	touch $@
	mkdir -p result/24kegg
	# 提取基因ID(Name)和KEGG基因ID(KgeneID)
	cut -f 1,2 temp/24kegg/gene_diamond.f6 | uniq | sed '1 i Name\tKgeneID' > temp/24kegg/gene_kegg.list
	# 追加KO编号(KO)和描述(Kdescription)
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2"\t"$$3} NR>FNR{print $$0,a[$$2]}' \
		/db/kegg/kegg_gene_ko_description.txt temp/24kegg/gene_kegg.list > temp/24kegg/gene_ko.list
	# 基因丰度矩阵末尾添加对应KO编号，没注释的直接删除，可选注释为unclassified
	# awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3} NR>FNR{print $0,a[$1]}' temp/24kegg/gene_ko.list result/23salmon_gene/gene.count | sed 's/\t$/\tunclassified/' > temp/24kegg/gene_ko.count
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$3} NR>FNR{print $$0,a[$$1]}' temp/24kegg/gene_ko.list result/23salmon_gene/gene.count | \
		sed '/\t$$/d' > temp/24kegg/gene_ko.count
	# 检查注释前后基因数量
	wc -l temp/24kegg/gene_ko.count
	wc -l result/23salmon_gene/gene.count 
	# 合并基因表为KO表，输出count值和tpm值
	Rscript ~/github/Metagenome/denovo1/script/mat_gene2ko.R -i temp/24kegg/gene_ko.count -o result/24kegg/kotab -n ${unit}
	# KO对应的描述
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print $$0,a[$$1]}' /db/kegg/ko_description.txt result/24kegg/kotab.count > result/24kegg/kotab.count.anno
 

### 2.4.3 CAZyome 碳水化合物数据库

24dbcan2: 23NRgeneSet
	touch $@
	mkdir -p temp/24dbcan2
	# 2.6M genes, 8 threads, 30min
	time diamond blastp --db ${dbcan2_dmnd} --query temp/23NRgene/protein.fa \
		--outfmt 6 --threads ${p} --max-target-seqs 1 --quiet \
		--out temp/24dbcan2/gene_diamond.f6

24dbcan2_sum: 24dbcan2
	mkdir -p result/24dbcan2
	# 提取基因对应基因家族，同一基因存在1对多，只取第一个
	cut -f 1,2 temp/24dbcan2/gene_diamond.f6 | uniq | sed 's/|/\t/g' | cut -f 1,3 | \
		cut -f 1,2 -d '_' |sed '1 i Name\tKO' > temp/24dbcan2/gene_fam.list
	# 基因丰度矩阵末尾添加对应FAM编号，没注释的直接删除
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print $$0,a[$$1]}' temp/24dbcan2/gene_fam.list result/23salmon_gene/gene.count | \
		sed '/\t$$/d' > temp/24dbcan2/gene_fam.count
	# 统计注释基因的比例
	wc -l result/23salmon_gene/gene.count
	wc -l temp/24dbcan2/gene_fam.count
	# 按基因家族合并
	Rscript ~/github/Metagenome/denovo1/script/mat_gene2ko.R -i temp/24dbcan2/gene_fam.count -o result/24dbcan2/cazytab
	# 结果中添加FAM注释，spf格式用于stamp分析
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print a[$$1],$$0}' ${dbcan2_anno} \
		result/24dbcan2/cazytab.count > result/24dbcan2/cazytab.count.spf


### 2.4.4 ResFams 抗生素抗性数据库
	
	# kraken2_gene
resfams: 
	touch $@
	mkdir -p temp/24resfams
	diamond blastp --db ${resfams_dmnd} --query temp/23prokka_all/mg.faa \
		--outfmt 6 --threads ${p} --max-target-seqs 1 --quiet \
		--out temp/24resfams/gene_diamond.f6

resfams_sum: resfams
	touch $@
	mkdir -p result/24resfams
	# 提取基因对应基因家族
	cut -f 1,2 temp/24resfams/gene_diamond.f6 | uniq |sed '1 i Name\tResGeneID' > temp/24resfams/gene_fam.list
	# 基因丰度矩阵末尾添加对应FAM编号，没注释的直接删除
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print a[$$1],$$0}' temp/24resfams/gene_fam.list result/23salmon_gene/gene.count | \
		sed '/^\t/d' > result/24resfams/resfam.count
	# 统计注释基因的比例
	wc -l result/23salmon_gene/gene.count
	wc -l result/24resfams/resfam.count # 172/7734=2.2%
	# 结果中添加FAM注释，spf格式用于stamp分析
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$4"\t"$$3"\t"$$2} NR>FNR{print a[$$1],$$0}' ${resfams_anno} \
		result/24resfams/resfam.count > result/24resfams/resfam.count.spf



# 3. 宏基因组分箱 Binning

## 3.1 混合分箱

binning:
	touch $@
	metawrap binning \
		-o temp/binning \
		-t ${p} \
		-a temp/megahit/final_assembly.fasta \
		--metabat2 --maxbin2 --concoct \
		seq/*.fastq


## 3.2 单样本分箱

binning_single: 
	touch $@
	metawrap binning \
		-o temp/binning \
		-t ${p} \
		-a temp/megahit/final_assembly.fasta \
		--metabat2 --maxbin2 --concoct \
		seq/*.fastq
	parallel --xapply -j ${j} \
		"cat temp/11qc/{1}*data_paired* > temp/12concat/{1}.fq" \
		::: `tail -n+2 result/metadata.txt | cut -f 1`

