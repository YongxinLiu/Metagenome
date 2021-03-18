#SHELL:=/bin/bash

	# 宏基因组分析流程第一版 
	# Metagenome pipeline version 1

# 1. 有参分析流程 Reference-based pipeline

	# 准备实验设计、测序数据和参数数据库
	# Prepare experiment design (result/metadata.txt), sequecing data (seq/*.fq.gz) and database

init:
	# Initial project
	# Clean zero-byte files of check points by touch (restart project)
	find . -name "*" -type f -size 0c | xargs -n 1 rm -f
	# Create directory
	mkdir -p seq temp result
	touch $@

## 1.1. 质控并移除宿主 Quality control & Remove host

### 1.1.1 (可选)质量评估 Quality access
	
qa: 
	touch $@
	# Quality access
	memusg -t fastqc -t ${p1} seq/*1.fq.gz
	memusg -t fastqc -t ${p1} seq/*2.fq.gz
	# Reports summary (multiqc call python3 from humann3 env)
	/conda/envs/humann3/bin/python /conda/envs/humann3/bin/multiqc -d seq/ -o result/qc/

### 1.1.2 质控并移除宿主 Quality control & Remove host

qc: 
	touch $@
	mkdir -p temp/qc result/qc
	# Quality control & Remove host in parallel
	# kneaddata call trimmomatic quality control and bowtie2 remove host
	# trimmomatic parameters ref：http://www.usadellab.org/cms/?page=trimmomatic
	# ILLUMINACLIP:去除接头文件TruSeq3-PE.fa:错配<=2:匹配>30:10 LEADING:切开头<3 TRAILING:切结尾<3 SLIDINGWINDOW:切4bp滑窗均值<15 MINLEN:72
	## 设置bt2参数示例：--bowtie2-options '--very-sensitive --dovetail'
	## --reorder排序，保证双端序列标签唯1时也可保持双端同序
	memusg -t parallel --xapply -j ${j} \
	  "kneaddata -i seq/{1}_1.fq.gz -i seq/{1}_2.fq.gz \
	  -o temp/qc -v -t ${p} --remove-intermediate-output \
	  --trimmomatic ${trimmomatic_path} \
	  --trimmomatic-options 'ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:72' \
	  --reorder -db ${host_bt2}" \
	  ::: `tail -n+2 result/metadata.txt | cut -f 1`
	# Summary Table
	kneaddata_read_count_table --input temp/qc --output temp/kneaddata.txt
	csvtk -t headers temp/kneaddata.txt
	cut -f 1,2,4,12 temp/kneaddata.txt | tail -n+2 | awk 'BEGIN{OFS=FS="\t"} {print $$0,$$3/$$2*100,$$4/$$3*100}' | \
	  sed 's/_1_kneaddata//' | \
	  sed '1 i SampleID\tRawReads\tTrimmedReads\tCleanReads\tQuality\tNonHost' \
	  > result/qc/kneaddata.txt
	cat result/qc/kneaddata.txt
	# 统计质控去宿主前后数据量(G)
	# csvtk -t summary -f RawReads:sum -i result/qc/kneaddata.txt|tail -n+2|awk '{print $0*150*2/1000000000}'
	tail -n+2 result/qc/kneaddata.txt|awk '{a=a+$$2}END{print a*150*2/1000000000}'
	tail -n+2 result/qc/kneaddata.txt|awk '{a=a+$$4}END{print a*150*2/1000000000}'
	# 检查输出双端文件是否成对(无报错行号信息即正确)
	seqkit seq -n -i temp/qc/${i}_1_kneaddata_paired_1.fastq|cut -f 1 -d '/' > temp/header_${i}_1
	seqkit seq -n -i temp/qc/${i}_1_kneaddata_paired_2.fastq|cut -f 1 -d '/' > temp/header_${i}_2
	cmp temp/header_${i}_?
	# 清理临时文件
	rm temp/qc/*contam*
	rm temp/qc/*unmatched* 

### 1.1.3 质控后再评估 Quality access after quality control

qa2: 
	touch $@
	mkdir -p submit
	# hard link to clean
	ln temp/qc/*data_paired* submit/
	# rename accroding to ID
	rename 's/_1_kneaddata_paired//;s/fastq/fq/' submit/*
	# compress for reducing space
	pigz submit/*
	ln result/kneaddata.txt submit/stat.txt
	ln result/metadata.txt submit/metadata.txt
	md5sum submit/*_1.fq.gz > submit/md5sum1.txt
	md5sum submit/*_2.fq.gz > submit/md5sum2.txt
	paste submit/md5sum1.txt submit/md5sum2.txt | \
	  awk '{print $$2"\t"$$1"\t"$$4"\t"$$3}' \
	  > submit/md5sum.txt
	sed -i 's/submit\///g' submit/md5sum.txt
	cat submit/md5sum.txt
	# Quality access
	memusg -t fastqc -t ${p1} submit/*1.fq.gz
	memusg -t fastqc -t ${p1} submit/*2.fq.gz
	# Summary quality control
	/conda/envs/humann3/bin/python /conda/envs/humann3/bin/multiqc -d submit/ -o result/qc/

### 1.1.4 (可选)Kraken2去宿主

kraken2_qc: 
	touch $@
	rm -rf temp/kraken2_qc
	mkdir -p temp/kraken2_qc
	memusg -t parallel -j 1 \
	  "echo {1} >> temp/kraken2_qc/log; \
	  kraken2 --db ${kraken2_db} --paired submit/{1}_*.${fq} \
	  --threads ${p} --use-names --report-zero-counts \
	  --report temp/kraken2_qc/{1}.report \
	  --output temp/kraken2_qc/{1}.output" >> temp/kraken2_qc/log 2>&1 \
	  ::: `tail -n+2 result/metadata.txt|cut -f1`

kraken2_qc_filter: 
	touch $@
	# 按Kraken2注释筛选
ifeq (${kraken2_qc_filter}, nohost)
	# 方案1. 提取非植物33090和动物(人)33208序列
	# 查看植物界Viridiplantae和动物界Metazoa含量
	grep -P 'Viridiplantae' temp/kraken2_qc/*.report
	grep -P 'Metazoa' temp/kraken2_qc/*.report
	# 370G,3jobs,24h; 330G,24x3p,22h;
	memusg -t parallel -j ${j} \
	  "extract_kraken_reads.py \
	  -k temp/kraken2_qc/{1}.output \
	  -r temp/kraken2_qc/{1}.report \
	  -1 submit/{1}_1.${fq}\
	  -2 submit/{1}_2.${fq} \
	  -t 33090 33208 --include-children --exclude \
	  --max ${read_max} --fastq-output \
	  -o temp/kraken2_qc/{1}_1.fq \
	  -o2 temp/kraken2_qc/{1}_2.fq" \
	  ::: `tail -n+2 result/metadata.txt|cut -f1`
else ifeq (${kraken2_qc_filter}, microbe)
	# 方案2. 提取细菌(2)、古菌(2157)、病毒(10239)和真菌(4751)，此数据库无真菌？
	# 查看植物界Viridiplantae和动物界Metazoa含量
	grep -P '\t    Bacteria$' temp/kraken2_qc/*.report
	grep -P '\t    Archaea$' temp/kraken2_qc/*.report
	grep -P '\t        Fungi$' temp/kraken2_qc/*.report # 无真菌吗？
	memusg -t parallel -j ${j} \
	  'memusg -t extract_kraken_reads.py \
	  -k temp/kraken2_qc/{1}.output \
	  -r temp/kraken2_qc/{1}.report \
	  -1 submit/{1}_1.fq.gz \
	  -2 submit/{1}_2.fq.gz \
	  -t 2 2157 10239 4751--include-children --fastq-output \
	  -o temp/kraken2_qc/{1}_1.fq \
	  -o2 temp/kraken2_qc/{1}_2.fq' \
	  ::: `tail -n+2 result/metadata.txt|cut -f1`
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of nohost or microbe")
endif
	echo -ne 'Kraken2 filter finished!!!\n' 

	# 链接/截取至qc区
	rm -rf temp/qc/*.fastq
	mkdir -p temp/qc/
ifeq (${subsample}, raw)
	# 链接质控后数据
	parallel -j 1 \
	  'ln temp/kraken2_qc/{1}_1.fq temp/qc/{1}_1_kneaddata_paired_1.fastq; \
	   ln temp/kraken2_qc/{1}_2.fq temp/qc/{1}_1_kneaddata_paired_2.fastq' \
	  ::: `tail -n+2 result/metadata.txt|cut -f1`
else ifeq (${subsample}, head)
	# 截取目标N行
	# 25s,1.5G,6m; 25s,6G,25m; 63s,6G,47m;2 23s,6G,55m; 
	memusg -t parallel -j 1 \
	  "head -n ${n} temp/kraken2_qc/{1}_1.fq > temp/qc/{1}_1_kneaddata_paired_1.fastq; \
	   head -n ${n} temp/kraken2_qc/{1}_2.fq > temp/qc/{1}_1_kneaddata_paired_2.fastq" \
	  ::: `tail -n+2 result/metadata.txt|cut -f1`
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of raw or head")
endif
	echo -ne 'subsample finished!!!\n' 
	pigz temp/kraken2_qc/*.fq

## 1.2. 物种和功能组成定量 humann2

### 1.2.1 humann2输入文件准备：双端文件cat合并

humann2_in: qc
	touch $@
	# 生成humann2输入要求的合并文件 cat pair-end for humann2
	# 截取最多6G数据分析
	mkdir -p temp/humann2_in
	parallel --xapply -j ${j} \
	  "head -n 80000000 temp/qc/{1}_1_kneaddata_paired_1.fastq > temp/humann2_in/{1}.fq" \
	  ::: `tail -n+2 result/metadata.txt | cut -f 1`
	# 查看样品数量和大小 show samples size
	ls -l temp/humann2_in/*.fq

### 1.2.2 humann2计算，包括metaphlan2

humann2: humann2_in
	touch $@
	mkdir -p temp/humann2
	memusg -t parallel -j ${j} \
	  'humann2 --input temp/humann2_in/{}.fq \
	  --threads ${p} \
	  --output temp/humann2/ ' \
	  ::: `tail -n+2 result/metadata.txt|cut -f1` 
	for i in `tail -n+2 result/metadata.txt|cut -f1`;do 
	  ln temp/humann2/$${i}_humann2_temp/$${i}_metaphlan_bugs_list.tsv temp/humann2/
	done
	# rm -rf temp/concat/* 
	# rm -rf temp/humann2/*_humann2_temp

### 1.2.3 功能组成整理 humann2_sum

humann2_sum: humann2
	touch $@
	mkdir -p result/humann2
	# 合并所有样品，通路丰度`pathabundance`包括各功能和具体的物种组成，还有基因家族`genefamilies`(太多)，通路覆盖度`pathcoverage`层面可以进分析
	humann2_join_tables \
	  --input temp/humann2/ \
	  --file_name pathabundance \
	  --output result/humann2/pathabundance.tsv
	sed -i 's/_Abundance//g' result/humann2/pathabundance.tsv
	# csvtk -t stat result/humann2/pathabundance.tsv
	# 标准化为相对丰度relab或百万分数cpm
	humann2_renorm_table \
	  --input result/humann2/pathabundance.tsv \
	  --units relab \
	  --output result/humann2/pathabundance_relab.tsv
	# 分层结果，结果stratified(每个菌的功能组成)和unstratified(功能组成)两个
	humann2_split_stratified_table \
	  --input result/humann2/pathabundance_relab.tsv \
	  --output result/humann2/
	sed -i 's/# Pathway/MetaCycPathway/' result/humann2/pathabundance_relab_*stratified.tsv


## 1.3. 整理物种组成表和基本绘图 Summary metaphlan2 and plot

### 1.3.1 整理物种组成表 Summary metaphlan2

metaphaln2_sum: 
	touch $@
	# metaphlan2功能组成
	mkdir -p result/metaphlan2
	# 合并单样品为表
	merge_metaphlan_tables.py temp/humann2/*_metaphlan_bugs_list.tsv | \
	  sed 's/_metaphlan_bugs_list//g' | sed '/^#/d' \
	  > result/metaphlan2/taxonomy.tsv
	# 转换为stamp的多级格式，株水不完整，不去掉列无法对齐
	metaphlan_to_stamp.pl result/metaphlan2/taxonomy.tsv \
	  > result/metaphlan2/taxonomy.spf
	# 绘制热图
	metaphlan_hclust_heatmap.py --in result/metaphlan2/taxonomy.tsv \
	  --out result/metaphlan2/tax_heat.pdf \
	  -c bbcry --top ${tax_top} --minv 0.1 -s log 
	# ValueError: The condensed distance matrix must contain only finite values.

### 1.3.2 GraPhlAn图

metaphaln2_graphlan: 
	touch $@
	# metaphlan2 to graphlan
	export2graphlan.py --skip_rows 1,2 \
	  -i result/metaphlan2/taxonomy.tsv \
	  --tree temp/metaphaln2_graphlan.tree \
	  --annotation temp/metaphaln2_graphlan.annot \
	  --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 \
	  --annotations 5,6 --external_annotations 7 --min_clade_size 1
	# graphlan annotation
	graphlan_annotate.py --annot temp/metaphaln2_graphlan.annot \
	  temp/metaphaln2_graphlan.tree \
	  temp/metaphaln2_graphlan.xml
	# output PDF figure, annoat and legend
	graphlan.py temp/metaphaln2_graphlan.xml \
	  result/metaphlan2/taxonomy_graphlan.pdf \
	  --external_legends --dpi 300 

### 1.3.3 物种组成LEfSe差异分析(可选)

	# 需要有完整的实验分组，且样本末尾数字为样本重复编号才可分析

metaphaln2_lefse: 
	touch $@
	# LEfSe差异分析和Cladogram
	# 修改样本品为组名，要求重复为结尾数值；下次改为metadata中通用代码
	sed '1 s/[0-9]*\t/\t/g' result/metaphlan2/taxonomy.tsv | grep -v '#' \
	  > result/metaphlan2/lefse.txt
	# 格式转换为lefse内部格式，需要安装于humann2中，或独立安装
	lefse-format_input.py result/metaphlan2/lefse.txt \
	  temp/lefse_input.in -c 1 -o 1000000
	# 运行lefse
	run_lefse.py temp/lefse_input.in \
	  temp/lefse_input.res
	# 绘制物种树注释差异
	lefse-plot_cladogram.py temp/lefse_input.res \
	  result/metaphlan2/lefse_cladogram.pdf --format pdf --dpi 600 
	# 绘制所有差异features柱状图
	lefse-plot_res.py temp/lefse_input.res \
	  result/metaphlan2/lefse_res.pdf --format pdf --dpi 600
	# 批量绘制所有差异features柱状图
	lefse-plot_features.py -f diff --archive none --format pdf \
	  temp/lefse_input.in temp/lefse_input.res \
	  result/metaphlan2/t_


## 1.4. kraken2物种组成

### 1.4.1 基于NCBI数据库的物种注释 Taxonomy assign by NCBI database

kraken2_read:
	touch $@
	rm -rf temp/kraken2_read
	mkdir -p temp/kraken2_read
	memusg -t parallel -j ${j} \
	  'echo {} >> temp/kraken2_read/log; \
	  kraken2 --db ${kraken2_db} --paired temp/qc/{1}_1_kneaddata_paired_*.fastq \
	  --threads ${p} --use-names --report-zero-counts \
	  --report temp/kraken2_read/{1}.report \
	  --output temp/kraken2_read/{1}.output' >> temp/kraken2_read/log 2>&1 \
	  ::: `tail -n+2 result/metadata.txt|cut -f1`
	mkdir -p result/kraken2_read
	grep -P 'Viridiplantae' temp/kraken2_read/*.report|cut -f 3 -d '/'|sed 's/.report:/\t/'|sed '1 i SampleID\tPercentage\tTotalReads\tRankReads\tRank\tTaxID\tTaxonomy' > result/kraken2_read/host.txt
	grep -P 'Metazoa' temp/kraken2_read/*.report|cut -f 3 -d '/'|sed 's/.report:/\t/' > result/kraken2_read/host.txt
	# format kraken2 report to metaphlan2 mpa
	for i in `tail -n+2 result/metadata.txt|cut -f1`;do
	  kreport2mpa.py -r temp/kraken2_read/$${i}.report --display-header \
	    -o temp/kraken2_read/$${i}.mpa; done
	# merge to table
	parallel -j ${j} \
	  'tail -n+2 temp/kraken2_read/{1}.mpa | sort | cut -f 2 | sed "1 s/^/{1}\n/" \
	  > temp/kraken2_read/{1}_count' \
	  ::: `tail -n+2 result/metadata.txt|cut -f1`
	## create colnames
	tail -n+2 temp/kraken2_read/${kraken2_header}.mpa | sort | cut -f 1 | \
	  sed "1 s/^/Taxonomy\n/" > temp/kraken2_read/0header_count
	# merge
	paste temp/kraken2_read/*count > result/kraken2_read/tax_count.mpa
	csvtk -t stat result/kraken2_read/tax_count.mpa

### 1.4.2 Bracken估计丰度

bracken2_read: 
	touch $@
	mkdir -p temp/bracken
	# estimate count in each rank
	for bracken2_tax in D P C O F G S;do
	for i in `tail -n+2 result/metadata.txt|cut -f1`;do
	bracken -d ${kraken2_db} \
	  -i temp/kraken2_read/$${i}.report \
	  -l $${bracken2_tax} -r ${bracken2_length} -t ${bracken2_t} \
	  -o temp/bracken/$${i}; done
	# Merge into table
	parallel -j 1 \
	  'tail -n+2 temp/bracken/{1}|sort|cut -f6| sed "1 s/^/{1}\n/" > temp/bracken/{1}.count ' \
	  ::: `tail -n+2 result/metadata.txt|cut -f1`
	tail -n+2 temp/bracken/${kraken2_header}|sort|cut -f1 | \
	  sed "1 s/^/Taxonomy\n/" > temp/bracken/0header.count
	# 检查文件数，为n+1
	ls temp/bracken/*count | wc
	# paste合并样本为表格，并删除非零行
	paste temp/bracken/*count > result/kraken2_read/bracken.$${bracken2_tax}.txt
	# 统计行列，默认去除表头
	csvtk -t stat result/kraken2_read/bracken.$${bracken2_tax}.txt
	done

## 1.5 Clean 清理有参过程临时文件

clean_ref: 
	touch $@
	# 1.1 质控qc质控后归档为submit，删除qc临时文件
	rm -rf temp/qc
	# 1.2 humman2筛选文件
	rm -rf temp/humann2_in
	# 备份临时文件中的metaphlan2重要结果
	mkdir temp/12metaphlan2
	cp -f temp/humann2/*_humann2_temp/*_metaphlan_bugs_list.tsv temp/12metaphlan2/
	# 删除humann2筛选文件
	rm -rf temp/humann2/*temp
	# 1.4 kraken2比对文件
	rm -r temp/kraken2_read/*output


# 2. 无参分析流程 De novo assemble pipeline

## 2.1. khmer质控(太耗时，结果对下游无明显改进，已删除)

## 2.2. Assemble 组装

### 2.2.1 序列混合组装(6-30个样品小项目)

megahit_all: 
	touch $@
	memusg -t megahit -t ${p1} \
	  -1 `tail -n+2 result/metadata.txt|cut -f1|sed 's/^/temp\/qc\//;s/$$/_1_kneaddata_paired_1.fastq/'|tr '\n' ','|sed 's/,$$//'` \
	  -2 `tail -n+2 result/metadata.txt|cut -f1|sed 's/^/temp\/qc\//;s/$$/_1_kneaddata_paired_2.fastq/'|tr '\n' ','|sed 's/,$$//'` \
	  -o temp/megahit_all --continue --min-contig-len ${min-contig-len} \
	  --k-min ${kmin} --k-max ${kmax} --k-step ${kstep} 
	seqkit stat temp/megahit_all/final.contigs.fa > temp/megahit_all/contigs.log
	cat temp/megahit_all/contigs.log

### 2.2.1 序列单样本组装(大项目 / 数据量>1/6内存)

assemble_single:
	touch $@
	echo -e "Assemble method: ${assemble_method}\nAssemble mode: assemble single\n"
	# ifseq, else, endif必须顶格写
ifeq (${assemble_method}, megahit)
	# 采用metahit单样品拼接
	# rm -rf temp/22megahit
	mkdir -p temp/22megahit
	memusg -t parallel -j ${j} \
		'megahit -t ${p1} --k-min ${kmin} --k-max ${kmax} --k-step ${kstep} \
		-1 temp/qc/{1}_1_kneaddata_paired_1.fastq \
		-2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
		-o temp/22megahit/{1} >> temp/22megahit/megahit.log 2>&1 ' \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
else ifeq (${assemble_method}, metaspades)
	# metaspades单样品拼接
	parallel --xapply -j ${j} \
		"metaspades.py -t ${p} -m 500 \
		-1 temp/qc/{1}_1_kneaddata_paired_1.fastq \
		-2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
		-o temp/22metaspades/{1}" \
		::: `tail -n+2 result/metadata.txt | cut -f 1`
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in megahit or metaspades")
endif
	echo -ne 'Assemble finished!!!\n' 

## 2.2.3 quast评估megahit混拼

megahit_all_quast: 
	touch $@
	# 拼接结果评估
	mkdir -p result/megahit_all
	memusg -t quast.py -o result/megahit_all/ \
	  -m ${quast_len} -t ${p} \
	  temp/megahit_all/final.contigs.fa
	# 备用重要中间文件至intermediate
	mkdir -p intermediate/megahit_all
	ln -f temp/megahit_all/final.contigs.fa intermediate/megahit_all/

### 2.2.5 (可选)kraken2重叠群物种注释 

kraken2_contig: 
	touch $@
	# 10M,9G,30m; 2G,4m; 7.7M,
	memusg -t kraken2 --db ${kraken2_db} \
	  temp/megahit_all/final.contigs.fa \
	  --threads ${p} --use-names \
	  --report temp/megahit_all/final.contigs.report \
	  --output temp/megahit_all/final.contigs.output
	grep -P 'Viridiplantae|Metazoa' temp/megahit_all/*.report
	# 从contig层面去宿主,10M,9G,18m,9.77/10.38M=94.1%; 2G,6m; 1.2G,1m,1775370/18.4M=96.2%
	memusg -t extract_kraken_reads.py \
	  -k temp/megahit_all/final.contigs.output \
	  -r temp/megahit_all/final.contigs.report \
	  -s temp/megahit_all/final.contigs.fa \
	  -t 33090 33208 --include-children --exclude \
	  -o temp/megahit_all/final.contigs.qc.fa 


## 2.3. Gene annotation 基因注释

### 2.3.1 对合并组装的单个contig文件基因注释

prodigal_all: 
	touch $@
	mkdir -p temp/prodigal_all
	# 2G,3h; 1.2G,1h; 6.7G,
	memusg -t prodigal -i temp/megahit_all/final.contigs.qc.fa \
	  -d temp/prodigal_all/gene.fa  \
	  -o temp/prodigal_all/gene.gff -p meta -f gff \
	  > temp/prodigal_all/gene.log 2>&1 
	ls -lsh temp/prodigal_all/gene.fa
	echo -ne 'metaProdigal all genes\t' > result/gene.log
	grep -c '>' temp/prodigal_all/gene.fa >> result/gene.log
	cat result/gene.log

prodigal_all_split: 
	touch $@
	# 文件太大，分块后再并行注释，然后合并
	seqkit split temp/megahit_all/final.contigs.qc.fa -s 1000000
	ls temp/megahit_all/final.contigs.qc.fa.split/final.contigs.qc.part_*.fa|cut -f 3 -d '_'|cut -f 1 -d '.' \
	  > temp/split.list
	mkdir -p temp/prodigal
	# 并行注释每个文件 1Gb/h
	memusg -t parallel -j 24 \
	  "prodigal -i temp/megahit_all/final.contigs.qc.fa.split/final.contigs.qc.part_{}.fa \
	  -d temp/prodigal/gene{}.fa  \
	  -o temp/prodigal/gene{}.gff -p meta -f gff \
	  > temp/prodigal/gene{}.log 2>&1 " \
	  ::: `cat temp/split.list`
	# 合并
	mkdir -p temp/prodigal_all
	cat temp/prodigal/gene*.fa > temp/prodigal_all/gene.fa
	cat temp/prodigal/gene*.gff > temp/prodigal_all/gene.gff

### 2.3.2 对单样品组装的每个contig文件基因注释(大数据可选)

prodigal_single: assemble_single
	touch $@
	echo -e "Assemble method: ${assemble_method}\nAssemble mode: assemble single\n"
	# ifseq, else, endif必须顶格写
ifeq (${assemble_method}, megahit)
	# 采用metahit单样品拼接
	memusg -t parallel --xapply -j ${j} \
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
	mkdir -p temp/NRgene
	cp temp/23prokka/mg.ffn temp/NRgene/
	echo -ne 'Prokka single genes\t' > result/gene.log
	grep -c '>' temp/NRgene/mg.ffn >> result/gene.log


### 2.3.3 构建非冗余基集 Non-redundancy gene set

NRgene: 
	touch $@
	mkdir -p temp/NRgene
	# -c sequence identity threshold, default 0.9; -M  max available memory (Mbyte), default 400; -l  length of throw_away_sequences, default 10
	# 4M,48p,2h; 11M,72p,7h
	memusg -t cd-hit-est \
	  -i temp/prodigal_all/gene.fa \
	  -o temp/NRgene/gene.fa \
	  -aS ${cdhit_coverage} -c ${cdhit_similarity} \
	  -G 0 -g 0 -M ${cdhit_mem} -T ${p1}
	# -n 10 -d 0 -g 1 -b ${band_width} 
	# 统计基因数据，确定ID是否非冗余
	echo -ne 'NRgene\t' >> result/gene.log
	grep -c '>' temp/NRgene/gene.fa >> result/gene.log
	cat result/gene.log
	# 翻译核酸为对应蛋白序列
	# 方法1. emboss, 14M,3m; 4M,1m
	transeq -sequence temp/NRgene/gene.fa \
	  -outseq temp/NRgene/protein.fa -trim Y 
	# 序列名自动添加了_1，为与核酸对应要去除, 14M,25s; 4M,36s
	sed -i 's/_1 / /' temp/NRgene/protein.fa
	# 方法2. seqkit，不加_1，结尾有*，14M,3m
	# memusg -t seqkit translate temp/NRgene/gene.fa > temp/NRgene/protein2.fa

### 2.3.4 基因定量 salmon genes

salmon_gene:
	touch $@
	mkdir -p temp/salmon_gene
	# 1. 建索引:-t转录本,-i索引; 14M,24p,83m; 4M,24p,15m
	memusg -t salmon index \
	  -t temp/NRgene/gene.fa \
	  -p ${p} \
	  -i temp/NRgene/salmon
	# 2. 定量: 38Gx4M,24px3,30m
	memusg -t parallel -j ${j} \
	  'salmon quant -i temp/NRgene/salmon -l A -p ${p} --meta \
	    -1 temp/qc/{1}_1_kneaddata_paired_1.fastq \
	    -2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
	    -o temp/salmon_gene/{1}.quant' \
	    ::: `tail -n+2 result/metadata.txt|cut -f1`
	# 3. 合并
	mkdir -p result/salmon_gene
	memusg -t salmon quantmerge --quants temp/salmon_gene/*.quant \
	  -o result/salmon_gene/gene.TPM
	memusg -t salmon quantmerge --quants temp/salmon_gene/*.quant \
	  --column NumReads \
	  -o result/salmon_gene/gene.count
	sed -i '1 s/.quant//g' result/salmon_gene/gene.*
	# 统计基因TPM， 11列   4,879,276行
	# csvtk -t stat result/salmon_gene/gene.TPM
	# 统计TPM和count值，TPM为1M，count值为比对的reads非整数
	# csvtk -t summary -f 2:sum result/salmon_gene/gene.TPM # 1M
	# csvtk -t summary -f 2:sum result/salmon_gene/gene.count # 9270775.02

### 2.3.5 (可选)基因物种注释 kraken2 annotate gene

kraken2_gene:
	touch $@
	# 9G,24p,3m;
	# 生成默认taxid结果
	memusg -t kraken2 --db ${kraken2_db} \
	  temp/NRgene/gene.fa \
	  --threads ${p} \
	  --report temp/NRgene/gene.report \
	  --output temp/NRgene/gene.output
	# 筛选基因-taxid列表
	grep '^C' temp/NRgene/gene.output|cut -f 2,3|sed '1 i Name\ttaxid' \
	  > temp/NRgene/gene.taxid
	# 基因添加7级注释
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$0} NR>FNR{print $$1,a[$$2]}' \
	  ${kraken2_db}/taxonomy.txt \
	  temp/NRgene/gene.taxid \
	  > temp/NRgene/gene.tax
	# 按7级合并; 4M,3m; 10M,7m
	mkdir -p result/kraken2_gene
	memusg -t ${py3}/python ${py3}/summarizeAbundance.py \
	  -i result/salmon_gene/gene.TPM \
	  -m temp/NRgene/gene.tax \
	  -c '2,3,4,5,6,7,8,9' -s ',+,+,+,+,+,+,+,' -n raw \
	  -o result/kraken2_gene/tax
	wc -l result/kraken2_gene/tax*|sort -n
	# 比较基因物种组成和读长物种组成
	# 从gene层面去宿主,9G,24m,13.5/13.92=96.9%; 
	#	memusg -t extract_kraken_reads.py \
	#	  -k temp/NRgene/gene.output \
	#	  -r temp/NRgene/gene.report \
	#	  -s temp/NRgene/gene.fa \
	#	  -t 33090 33208 --include-children --exclude \
	#	  -o temp/NRgene/gene.qc.fa 

## 2.4 功能数据库注释

### 2.4.1 eggNOG

eggnog: 
	touch $@
	mkdir -p temp/eggnog
	# Blast NRgene against eggNOG database by diamond
	memusg -t emapper.py -m diamond --no_annot --no_file_comments \
	  --data_dir ${eggnog_db} --cpu ${p1} -i temp/NRgene/protein.fa \
	  -o temp/eggnog/protein --override
	# Annotation blast result
	memusg -t emapper.py --annotate_hits_table \
	  temp/eggnog/protein.emapper.seed_orthologs --no_file_comments \
	  -o temp/eggnog/output --cpu ${p1} --data_dir ${eggnog_db} --override

eggnog_sum: 
	touch $@
	mkdir -p result/eggnog
	# Add header
	sed '1 i Name\tortholog\tevalue\tscore\ttaxonomic\tprotein\tGO\tEC\tKO\tPathway\tModule\tReaction\trclass\tBRITE\tTC\tCAZy\tBiGG\ttax_scope\tOG\tbestOG\tCOG\tdescription' temp/eggnog/output.emapper.annotations > temp/eggnog/output
	echo -n -e "Total protein\t " > result/eggnog/log
	grep -c '>' temp/NRgene/protein.fa >> result/eggnog/log
	echo -n -e "Annotated protein\t" >> result/eggnog/log
	wc -l temp/eggnog/output.emapper.annotations | cut -f1 -d ' ' >> result/eggnog/log
	echo -n -e "KO Annotated protein\t" >> result/eggnog/log
	grep -c 'ko:K' temp/eggnog/output  >> result/eggnog/log
	# group by annotation, generate sum.CAZy.raw.txt  sum.COG.raw.txt  sum.KO.raw.txt
	memusg -t ${py3}/python ${py3}/summarizeAbundance.py \
	  -i result/salmon_gene/gene.TPM \
	  -m temp/eggnog/output \
	  -c '9,10,11,16,21' -s ',+,+,+,+*' -n raw \
	  -o result/eggnog/sum
	wc -l result/eggnog/sum.*
	# KO summary
	sed -i 's/^ko://' result/eggnog/sum.KO.raw.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print a[$$1],$$0}' \
	  ${eggnog_db}/KO_description.txt \
	  result/eggnog/sum.KO.raw.txt| \
	  sed 's/^\t/Description\t/' > result/eggnog/sum.KO.TPM.spf
	# CAZy
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print a[$$1],$$0}' \
	   ${eggnog_db}/fam_description.txt result/eggnog/sum.CAZy.raw.txt | \
	  sed 's/^\t/Description\t/' > result/eggnog/sum.CAZy.TPM.spf
	# COG
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2"\t"$$3} NR>FNR{print a[$$1],$$0}' \
	  ${eggnog_db}/COG.anno \
	  result/eggnog/sum.COG.raw.txt | sed '1 s/^/Level1\tLevel2/'> \
	  result/eggnog/sum.COG.TPM.spf

ko_sum: 
	# 通路添加2，1级注释，有很多map和ko无法注释
	# 无法注释的ko在KEGG中也查不到，map号可查但无层级关系，暂时过滤后428个通路
	mv result/eggnog/sum.pathway.raw.txt result/eggnog/sum.Pathway.raw.txt
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$4]=$$1"\t"$$2"\t"$$3} NR>FNR{print a[$$1],$$0}' \
	 /db/kegg/brite/level1-3ko.txt result/eggnog/sum.Pathway.raw.txt  > result/eggnog/sum.Pathway.raw.spf
	# KO2pathway
	${py3}/python ${py3}/summarizeAbundance.py \
	  -i result/eggnog/sum.KO.raw.txt \
	  -m /db/kegg/brite/ko_path.list \
	  -c '2' -s '\t' -n raw \
	  -o result/eggnog/ko
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$3]=$$1"\t"$$2} NR>FNR{print a[$$1],$$0}' \
	 /db/kegg/brite/KO1-3.list result/eggnog/ko.Pathway.raw.txt > result/eggnog/ko.Pathway.raw.spf

### 2.4.2 KEGG (可选)

kegg: NRgene
	touch $@
	mkdir -p temp/kegg
	# 2.6M,8p,6h; 5M,72p,30m; 2.5M,72p,16m(--sensitive 108m,慢7倍)； 
	# KEGG基因注释
	memusg -t diamond blastp --db ${kegg_dmnd} \
	  --query temp/NRgene/protein.fa \
	  --outfmt 6 --threads ${p1} --quiet \
	  --evalue 1e-5 --max-target-seqs 1 --sensitive \
	  --out temp/kegg/gene_diamond.f6
	# 注释基因数量3853393/4879276=79%; 1919238/2550307=75%; 1413301/1780779=79%
	wc -l temp/kegg/gene_diamond.f6

kegg_sum: 
	touch $@
	mkdir -p result/kegg
	# 提取基因ID(Name)和KEGG基因ID(KgeneID)
	cut -f 1,2 temp/kegg/gene_diamond.f6 | uniq | sed '1 i Name\tKgeneID' \
	  > temp/kegg/gene_kegg.list
	# 添加基因ID对应ID，存在1对多
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print $$0,a[$$2]}' \
	  ${db}/kegg/genes/ko/gene2ko.txt \
	  temp/kegg/gene_kegg.list \
	  > temp/kegg/gene2ko.list
	# 合并基因-KO, 2.5M,1m; 
	memusg -t summarizeAbundance.py \
	  -i result/salmon_gene/gene.TPM \
	  -m temp/kegg/gene2ko.list \
	  -c 3 -s ',' -n raw \
	  -o result/kegg/TPM
	# csvtk -t summary -f 2:sum result/kegg/TPM.KO.raw.txt # 255747
	# 合并KO-Pathway
	memusg -t ${py3}/python ${py3}/summarizeAbundance.py \
	  -i result/kegg/TPM.KO.raw.txt \
	  -m /db/kegg/brite/KO1-4.txt \
	  -c 2,3,4 -s ',+,+,' -n raw \
	  -o result/kegg/TPM
	# KO表添加注释
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print a[$$1],$$0}' \
	  ${eggnog_db}/KO_description.txt \
	  result/kegg/TPM.KO.raw.txt \
	  > result/kegg/TPM.KO.raw.spf
	wc -l result/kegg/TPM*

### 2.4.3 CAZy 碳水化合物数据库

dbcan2: 
	touch $@
	mkdir -p temp/dbcan2
	# 3M,8p,0.5h; 14M,24p,1.5h; 2.5M,72p,
	# http://bcb.unl.edu/dbCAN2, diamond默认e值为1e-102
	memusg -t diamond blastp --db ${dbcan2_dmnd} \
	  --query temp/NRgene/protein.fa \
	  --outfmt 6 --threads ${p1} --max-target-seqs 1 \
	  --sensitive -e 1e-102 --quiet \
	  --out temp/dbcan2/gene_diamond.f6

dbcan2_sum: dbcan2
	touch $@
	# 整理比对数据为表格
	mkdir -p result/dbcan2
	# 提取基因与dbcan分类对应表，1.5M,12s
	perl /db/script/format_dbcan2list.pl \
	  -i temp/dbcan2/gene_diamond.f6 \
	  -o temp/dbcan2/gene.list
	# 按对应表累计丰度，14M,8m; 2.5M,39m
	memusg -t ${py3}/python ${py3}/summarizeAbundance.py \
	  -i result/salmon_gene/gene.TPM \
	  -m temp/dbcan2/gene.list \
	  -c 2 -s ',' -n raw \
	  -o result/dbcan2/TPM
	# 添加注释生成STAMP的spf格式，结合metadata.txt进行差异比较
	awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$$1]=$$2} NR>FNR{print a[$$1],$$0}' \
	  ${db}/dbcan2/fam_description.txt result/dbcan2/TPM.CAZy.raw.txt | \
	  sed 's/^\t/Unannotated\t/' > result/dbcan2/TPM.CAZy.raw.spf
	# 检查未注释数量，少量正常，大量则需要检查原因
	grep 'Unannotated' result/dbcan2/TPM.CAZy.raw.spf|wc -l

### 2.4.4 CARD 抗生素抗性数据库
	
card: 
	touch $@
	# 加载数据库
	rgi load -i ${db}/card/card.json \
	  --card_annotation ${db}/card/card_database_v3.1.0.fasta --local
	mkdir -p result/card
	# card无法自动截取ID，需手动截取为纯ID
	cut -f 1 -d ' ' temp/NRgene/protein.fa > temp/protein.fa
	# --include_loose 会显著增加上百倍结果，14M,24p,30m; 2.5M,24p,6m
	memusg -t rgi main \
	  --input_sequence temp/protein.fa \
	  --output_file result/card/rgi \
	  --input_type protein --clean \
	  --num_threads ${p1} --alignment_tool DIAMOND
	# 筛选结果和注释，json可直接在CARD在线可视化，txt可本地分析
	# 查看表头，1基因ID/9抗性分类ARO/15药物分类(19)/16抗性机制(6)/17基因家族(26)
	# csvtk -t headers result/card/rgi.txt
	# tail -n+2 result/card/rgi.txt|cut -f 17|sort|uniq|wc
	sed -i '1 s/ORF_ID/Name/;s/ //g' result/card/rgi.txt
	memusg -t ${py3}/python ${py3}/summarizeAbundance.py \
	  -i result/salmon_gene/gene.TPM \
	  -m result/card/rgi.txt \
	  -c '9,15,16,17' -s ';+;+;+;' -n raw \
	  -o result/card/ARG

# 3. 宏基因组分箱 Binning

## 3.1 混合分箱

	# 分箱中输入文件不配对，检查文件ID是否一致，见manual.md ### 检查双端配对
bin_all:
	touch $@
	mkdir -p temp/bin_all
	# 筛选序列用于分箱，200bp-1000bp,5.6G-2.3G
	seqkit stat temp/megahit_all/final.contigs.fa
	seqkit seq -m ${l} temp/megahit_all/final.contigs.fa > temp/bin_all/final.contigs.fa
	seqkit stat temp/bin_all/final.contigs.fa
	memusg -t metawrap binning \
	  -o temp/bin_all \
	  -t ${p1} \
	  -a temp/bin_all/final.contigs.fa \
	  --metabat2 --maxbin2 --concoct \
	  temp/qc/*.fastq > temp/bin_all/log
	tail temp/bin_all/log
	rm temp/bin_all/log

## 3.2 分组分箱

bin_batch:
	touch $@
	mkdir -p temp/bin
	for g in ${bin_list};do 
	# rm -r temp/bin/b_metahit_${g}
	memusg -t megahit -t ${p1} \
	  -1 `grep -P "\t${g}\t" result/metadata.txt|cut -f1|sed 's/^/temp\/qc\//;s/$$/_1_kneaddata_paired_1.fastq/'|tr '\n' ','|sed 's/,$$//'` \
	  -2 `grep -P "\t${g}\t" result/metadata.txt|cut -f1|sed 's/^/temp\/qc\//;s/$$/_1_kneaddata_paired_2.fastq/'|tr '\n' ','|sed 's/,$$//'` \
	  -o temp/bin/b_metahit_${g} --continue --min-contig-len ${l} \
	  --k-min ${kmin} --k-max ${kmax} --k-step ${kstep} 
	done

## 3.3 单样本分箱

binning_single: 
	touch $@
	metawrap binning \
		-o temp/binning \
		-t ${p} \
		-a temp/megahit/final_assembly.fasta \
		--metabat2 --maxbin2 --concoct \
		seq/*.fastq
	parallel --xapply -j ${j} \
		"cat temp/qc/{1}*data_paired* > temp/humann2_in/{1}.fq" \
		::: `tail -n+2 result/metadata.txt | cut -f 1`


help:
	# 帮助文档: 流程所需的软件、脚本及数据库版本
	# Help: Pipeline dependency version of softwares, scripts and databases
	
	# 软件 Softwares
	# 绿色软件直接默认安装或放在系统环境变量
	fastp -v # version 0.12.5 fastq文件质控、质量类型转换
	usearch10 --help # v10.0.240_i86linux64 扩增子分析软件
	clustalo --version # 1.2.1 多序列对齐
	biom --version # 2.1.5， OTU表格式转换

	# conda软件分布于meta、metaRef、metawrap、drep、gtdbtk等环境
	conda activate meta
		fastqc -v # 0.11.9 测序数据质量评估
		multiqc --version # 1.8 质量报告汇总
		kneaddata --version # v0.7.4 质控和去宿主
			kneaddata_read_count_table -h # 质控合并
		parallel --version # 20201122 并行任务管理
		run_lefse.py -h # 1.0 生物标志物鉴定
			lefse-format_input.py # 格式转换
			lefse-plot_cladogram.py # 绘制树状图
			lefse-plot_res.py # 所有差异柱状图
			lefse-plot_features.py # 绘制单特征图
		kraken2 --version # 2.0.9-beta
		megahit -h # v1.1.3
		quast.py # 5.0.0
		salmon -v # 非比对基因定量软件 0.11.2
		prokka --version # 基因组注释 1.13.3
		cd-hit # 去冗余 
			cd-hit-est # 核酸水平去冗余
		diamond help # 类blast比对工具 v0.8.22.84
	conda activate humann2
		humann2 --version # v2.8.1 功能组成定量
			humann2_join_tables # 样本功能丰度合并为特征表
			humann2_renorm_table # 特征表重新标准为relab(1)或cpm
			humann2_split_stratified_table # 生成功能和功能物种组成表
		metaphlan2.py --version # 2.7.5 物种组成定量
			merge_metaphlan_tables.py -h # 样本物种丰度合并为特征表
			metaphlan_hclust_heatmap.py -h # 样本特征表画热图
		export2graphlan.py -h # ver. 0.22 of 05 May 2020 mpa转换为graphlan输入
		graphlan.py -h # 1.1.3 (5 June 2018) # 绘制树图
			graphlan_annotate.py # 添加树注释

	# 脚本存放于/home/meta/soft/Metagenome/denovo1/script/ Scripts
	# microbiome_helper # 2017, 输助脚本 https://github.com/mlangill/microbiome_helper
		metaphlan_to_stamp.pl # 转换mpa为spf格式
	# Amplicon http://github.com/yongxinliu/Amplicon
		alpha_boxplot.sh -h # 1.0 基于usearch alpha_div绘制箱线图
	# Metagenome http://github.com/yongxinliu/Metagenome
		kegg_ko00001_htext2tsv.pl # 生成KEGG层级注释

	# 数据库 Databases


# 更新日志 Update log
	# 2018/09/12 Add kneaddata, humann2 for reference based pipeline
	# 2018/09/23 Add kraken2 for taxonomy in reads levels, salmon genes quaitity
	# 2019/06/06 Add de novo pipeline
	# 2020/12/30 Test on medicago root and gut metagenome
	# 2021/1/10 Update eggnog to emapper2, and add KEGG; remove number in pipeline
	# 2021/1/18 add bracken for kraken2 summay in each level
