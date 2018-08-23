SHELL:=/bin/bash

	# 16S扩增子分析流程第二版 16S Amplicon pipeline version 2

	# 帮助文档
	
	# 流程所需的软件、脚本及数据库版本
	# Help: Pipeline dependency version of softwares, scripts and databases

version: 
	# 2018/4/27 2.0 Standard 16S anlysis report
	# 
	# 软件 Softwares
	# fastqc -v # 0.11.5 测序数据质量评估
	# parallel # 并行/多线程任务管理
	# fastp -v # version 0.12.5 fastq文件质控、质量类型转换
	# usearch10 --help # v10.0.240_i86linux64 扩增子分析软件
	# clustalo --version # 1.2.1 多序列对齐
	# filter_alignment.py # from qiime1.9.1 筛选比对序列及区域
	# make_phylogeny.py # from qiime1.9.1 默认调用fastree建树
	# biom --version # 2.1.5， OTU表格式转换
	# 
	# 脚本 Scripts
	# alpha_boxplot.sh # 1.0 基于usearch alpha_div绘制箱线图
	#
	# 数据库 Databases
	# greengene13_5 # reference, and/or taxonomy database 
	# rdp train set 16 # taxonomy database
	# silva132 # chimera reference, and/or taxonomy database 


# 1. 标准流程 Standard pipeline

	# 建立程序必须目录 Create work directory

init:
	touch $@
	mkdir -p seq doc temp result script


## 1.1. 拆分下机数据为文库
	
	# Split lane into library
	# 需要将lane文件放入seq目录，对应的index和文库放在doc/library.txt
	# 注意：makefile中注释行的#顶格则不输出，如果缩进则输出
	# 参考Illumina index列表见doc/IlluminaIndex.fa，我们常用的为反向互补序列(revseq命令)为IlluminaIndexRC.fa

lane_split:
	touch $@
	# 质量原始数据 Quality control of lane files
	fastqc -t ${p} seq/${lane}_* &
	# 并行处理任务列表 parallel grep each index
	# parallel并行任务，--xapply参数平行而非相乘组合，-j进程数；zcat解压，grep -A匹配barocde，-v删除--行，输出lane名
	parallel --xapply -j ${p} "zcat seq/lane_1.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$$' > seq/{2}_1.fq" \
		::: `tail -n+2 doc/library.txt | cut -f 2` ::: `tail -n+2 doc/library.txt | cut -f 1`
	parallel --xapply -j ${p} "zcat seq/lane_2.fq.gz | grep -A 3 '#{1}'| grep -v -P '^--$$' > seq/{2}_2.fq" \
		::: `tail -n+2 doc/library.txt | cut -f 2` ::: `tail -n+2 doc/library.txt | cut -f 1`
	# 统计每个文库序列数量
	echo -e "libraryID\treads" > ${lib_log}
	for l in `tail -n+2 doc/library.txt | cut -f 1 | tr '\n' ' '`; do \
	echo -ne "$${l}\t" >> ${lib_log}; wc -l seq/$${l}_1.fq | awk '{print $$1/4}' >> ${lib_log}; done
	cat ${lib_log}

## 1.2 拆分文库为样品 Split library into sample

	# makefile中使用for循环，调用的列表要在同一行用空格分隔，换行分隔会报错；代码也要在同一行用分号断句，不可用回车，但可用\换行
	# 默认只拆分单左端barcode类型的样品:先匹配左端，再提取序列ID，再提取右端，最后改名，注意实验设计要严格规范无空格

library_split:
	touch $@
	# 文库按barcode拆分样品，保存至seq/sample目录中
	# s/^@//;去掉序列名标记；s/1$$/2/;修改/1端编号格式；s/ 1/ 2/; 修改 1端编号格式
	mkdir -p seq/sample
	for l in `tail -n+2 doc/library.txt | cut -f 1 | tr '\n' ' '`; do \
	echo $${l}; \
	parallel --xapply -j ${p} "grep -A 2 -B 1 -P "^{1}" seq/$${l}_1.fq | grep -v -P '^--$$' \
		> seq/sample/{2}_1.fq" ::: `tail -n+2 doc/$${l}.txt | cut -f 2` ::: `tail -n+2 doc/$${l}.txt | cut -f 1`; \
	parallel -j ${p} "awk 'NR%4==1' seq/sample/{1}_1.fq | sed 's/^@//;s/1$$/2/;s/ 1/ 2/' \
		> temp/id.{1}" ::: `tail -n+2 doc/$${l}.txt | cut -f 1`; \
	parallel -j ${p} "usearch10 -fastx_getseqs seq/$${l}_2.fq -labels temp/id.{1} \
		-fastqout seq/sample/{1}_2.fq" ::: `tail -n+2 doc/$${l}.txt | cut -f 1`; \
	done

library_split_stat: library_split
	touch $@
	rm -fr result/split
	mkdir -p result/split
	# 统计每个文库中样品测序数量
	# plot_bar_library.sh添加-A可修改组名列、-l可设置图片中显示正确的库编号
	for l in `tail -n+2 doc/library.txt | cut -f 1 | tr '\n' ' '|sed 's/$$ //'`; do echo $${l}; \
	for s in `tail -n+2 doc/$${l}.txt | cut -f 1 | tr '\n' ' '`; do \
	echo -ne "$${s}\t" >> result/split/$${l}.txt; wc -l seq/sample/$${s}_1.fq | awk '{print $$1/4}' >>  result/split/$${l}.txt; done; \
	plot_bar_library.sh -i result/split/$${l}.txt -d doc/design.txt -o result/split/$${l} -A ${g1} -l $${l};\
	done
	cat result/split/L*.txt > ${sample_split}
	# 实验设计样本数
	wc -l doc/design.txt
	# 显示总样品数
	ls seq/sample/* | wc -l


## 1.3 双端序列合并 Merge pair-end reads

sample_merge:
	touch $@
	# 双端序列合并
	mkdir -p seq/merge
	parallel -j ${p} "usearch10 -fastq_mergepairs seq/sample/{1}_1.fq -reverse seq/sample/{1}_2.fq \
		-fastqout seq/merge/{1}.fq -relabel {1}. -threads 1" ::: `tail -n+2 doc/design.txt | cut -f 1`
	# 合并所有双端合并的样品
	cat seq/merge/* > seq/all.fq

sample_merge_stat: sample_merge
	touch $@
	# 统计每个样品merge前后reads
	rm -f ${sample_merge}
	touch ${sample_merge}
	for s in `tail -n+2 doc/design.txt | cut -f 1 | tr '\n' ' '`; do \
	echo -ne "$${s}\t" >> ${sample_merge}; wc -l seq/merge/$${s}.fq | awk '{print $$1/4}' >> ${sample_merge}; done


## 1.4 切除引物 Cut primers and quality filter
# Cut barcode 10bp + V5 19bp in left and V7 18bp in right
fq_trim: sample_merge_stat
	touch $@
	usearch10 -fastx_truncate seq/all.fq \
		-stripleft ${stripleft} -stripright ${stripright} \
		-fastqout temp/stripped.fq 

## **1.5 质量控制fastq filter**
# 默认过滤错误率超1%的序列 Keep reads error rates less than 1%
fq_qc: fq_trim
	touch $@
	usearch10 -fastq_filter temp/stripped.fq \
		-fastq_maxee_rate ${fastq_maxee_rate} \
		-fastaout temp/filtered.fa -threads ${p}



## 1.6 序列去冗余 Remove redundancy
# miniuniqusize为8，去除低丰度，增加计算速度
fa_unqiue: 
	touch $@
	usearch10 -fastx_uniques temp/filtered.fa \
		-minuniquesize ${minuniquesize} -sizeout \
		-fastaout temp/uniques.fa -threads ${p}
	echo -ne 'Unique reads\t' > ${otu_log}
	grep -c '>' temp/uniques.fa >> ${otu_log}
	cat ${otu_log}


## 1.7 挑选OTU Pick OTUs

	# 可选97% cluster_otus，或100% unoise3，默认unoise3，不支持多线程
	# ifeq条件必须顶格，否则报错

otu_pick: fa_unqiue
	touch $@
	echo -e "OTU method\t${otu_method}" >> ${otu_log}
ifeq (${otu_method}, unoise3)
	# 类似100%聚类，只去除嵌合体和扩增及错误，保留所有高丰度序列
	usearch10 -unoise3 temp/uniques.fa -zotus temp/Zotus.fa -minsize ${minuniquesize} -threads ${p}
else ifeq (${otu_method}, cluster_otus)
	# cluster_otus无法修改聚类参数，想使用不同聚类相似度，使用cluster_smallmem命令
	usearch10 -cluster_otus temp/uniques.fa -otus temp/Zotus.fa
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in unoise3 or cluster_otus") 
endif
	awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' temp/Zotus.fa > temp/otus.fa
	echo -ne 'OTU number\t' >> ${otu_log}
	grep -c '>' temp/otus.fa >> ${otu_log}
	cat ${otu_log}


## 1.8 基于参考序列去嵌合 Remove chemira by silva

	# 可选，推荐使用最新SILVA大数据库 https://www.arb-silva.de/
	# 官方模式推荐sensitive错，推荐balanced, high_confidence

chimera_ref: otu_pick
	touch $@
	usearch10 -uchime2_ref temp/otus.fa \
		-db ${chimera_ref} -strand plus -mode ${chimera_mode} \
		-chimeras temp/otus_chimeras.fa -threads ${p}
	# 获得非嵌合体序列ID
	cat temp/otus.fa temp/otus_chimeras.fa| grep '>' | sort | uniq -u | sed 's/>//' > temp/no_chimeras.id
	# 筛选非嵌合体
	usearch10 -fastx_getseqs temp/otus.fa -labels temp/no_chimeras.id -fastaout temp/otus_no_chimeras.fa
	echo -ne 'no chimeras\t' >> ${otu_log}
	grep -c '>' temp/otus_no_chimeras.fa >> ${otu_log}
	cat ${otu_log}


## 1.9 去除宿主 remove host

host_rm: chimera_ref
	touch $@
ifeq (${host_method}, blast)
	# 方法1. 基于宿主基因组(含叶绿体/线粒体)比对
	blastn -query temp/otus_no_chimeras.fa -db ${host} -out temp/otus_no_chimeras.blastn \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' \
	-num_alignments 1 -evalue 1 -num_threads ${p}
	awk '$$3>${host_similarity} && $$13>${host_coverage}' temp/otus_no_chimeras.blastn | cut -f 1 | sort | uniq > temp/otus_host.id
	cat <(grep '>' temp/otus_no_chimeras.fa|sed 's/>//') temp/otus_host.id | sort | uniq -u > temp/otus_no_host.id
else ifeq (${host_method}, sintax_gg)
	# 方法2. 基于GG13_5的usearch注释结果筛选，排除线粒体、叶绿体和非16S
	usearch10 -sintax temp/otus_no_chimeras.fa \
		-db ${usearch_gg} -sintax_cutoff ${sintax_cutoff} -strand both \
		-tabbedout temp/otus_no_chimeras.tax -threads ${p}
	grep -P -v 'mitochondria|Chloroplast|\t$$' temp/otus_no_chimeras.tax | cut -f 1 > temp/otus_no_host.id
else ifeq (${host_method}, sintax_silva)
	# 方法3. 基于silva132的usearch注释结果筛选，排除线粒体、叶绿体、真核和非16S
	usearch10 -sintax temp/otus_no_chimeras.fa \
		-db ${usearch_silva} -sintax_cutoff ${sintax_cutoff} -strand both \
		-tabbedout temp/otus_no_chimeras.tax -threads ${p}
	grep -P -v 'Mitochondria|Chloroplast|Eukaryota|\t$$' temp/otus_no_chimeras.tax | cut -f 1 > temp/otus_no_host.id
else ifeq (${host_method}, sintax_silva_its)
	# 方法4. 基于silva132的usearch注释ITS，排除线粒体、叶绿体、细菌16S
	usearch10 -sintax temp/otus_no_chimeras.fa \
		-db ${usearch_silva} -sintax_cutoff ${sintax_cutoff} -strand both \
		-tabbedout temp/otus_no_chimeras.tax -threads ${p}
	grep -P -v 'Mitochondria|Chloroplast|Bacteria|\t$$' temp/otus_no_chimeras.tax | cut -f 1 > temp/otus_no_host.id
else ifeq (${host_method}, sintax_unite)
	# 方法5. 基于blast和usearch注释结果筛选，排除宿主和非真菌
	# blastn去除宿主
	blastn -query temp/otus_no_chimeras.fa -db ${host} -out temp/otus_no_chimeras.blastn \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' \
	-num_alignments 1 -evalue 1 -num_threads ${p}
	awk '$$3>${host_similarity} && $$13>${host_coverage}' temp/otus_no_chimeras.blastn | cut -f 1 | sort | uniq > temp/otus_host.id
	# usearch使用unite注释ITS排除宿主和非真菌
	usearch10 -sintax temp/otus_no_chimeras.fa \
		-db ${usearch_unite} -sintax_cutoff ${sintax_cutoff} -strand both \
		-tabbedout temp/otus_no_chimeras.tax -threads ${p}
	grep -P -v 'Fungi|\t$$' temp/otus_no_chimeras.tax | cut -f 1 > temp/otus_no_fungi.id
	cat <(grep '>' temp/otus_no_chimeras.fa|sed 's/>//') temp/otus_host.id temp/otus_no_fungi.id | sort | uniq -u > temp/otus_no_host.id
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in blast, usearch_gg or usearch_silva") 
endif
	# 最终筛选结果复制入result
	usearch10 -fastx_getseqs temp/otus_no_chimeras.fa -labels temp/otus_no_host.id -fastaout temp/otus_no_host.fa
	echo -ne 'no host\t' >> ${otu_log}
	grep -c '>' temp/otus_no_host.fa >> ${otu_log}
	cat ${otu_log}
	cp temp/otus_no_host.fa result/otu.fa


## 1.10 生成OTU表 Creat OTUs table

otutab_create: host_rm
	touch $@
ifeq (${map_method}, usearch10)
	# 方法1：usearch10, temp/stripped.fq建议新流程使用，与旧体系只有filtered.fa兼容
	usearch10 -otutab temp/filtered.fa -otus result/otu.fa -id ${map_identify} \
		-otutabout temp/otutab.txt -threads ${p}
else ifeq (${map_method}, vsearch)
	# 方法2：vsearch，比usearch10更快
	time vsearch --usearch_global temp/filtered.fa --db result/otu.fa --id ${map_identify} \
		--otutabout temp/otutab.txt --threads ${p}
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in usearch10 or vsearch") 
endif


## 1.11 OTU表筛选 Filter OTU table

otutab_filter: otutab_create
	touch $@
	# 统计OTU表的基本信息
	usearch10 -otutab_stats temp/otutab.txt -output temp/otutab.txt.stat
	echo -ne "\nOTU table summary\nType\tThreshold\tReads\tSamples\tOTUs\nTotal\t-\t" > ${log_otutable}
	head -n3 temp/otutab.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	# 按样本测序量筛选：通常低于5000的样本会删除
	usearch10 -otutab_trim temp/otutab.txt -min_sample_size ${min_sample_size} -output temp/otutab_trim1.txt
	usearch10 -otutab_stats temp/otutab_trim1.txt -output temp/otutab_trim1.txt.stat
	echo -ne "SampleSize\t${min_sample_size}\t" >> ${log_otutable}
	head -n3 temp/otutab_trim1.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	# 按OTU测序量筛选：通常低于8的OTU会删除
	usearch10 -otutab_trim temp/otutab_trim1.txt -min_otu_size ${min_otu_size} -output temp/otutab_trim2.txt
	usearch10 -otutab_stats temp/otutab_trim2.txt -output temp/otutab_trim2.txt.stat
	echo -ne "OtuSize\t${min_otu_size}\t" >> ${log_otutable}
	head -n3 temp/otutab_trim2.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	# 按OTU相对丰度筛选：通常低于1 RPM 的OTU会删除
	usearch10 -otutab_trim temp/otutab_trim2.txt -min_otu_freq ${min_otu_freq} -output temp/otutab_trim3.txt
	usearch10 -otutab_stats temp/otutab_trim3.txt -output temp/otutab_trim3.txt.stat
	echo -ne "OtuFreq\t${min_otu_freq}\t" >> ${log_otutable}
	head -n3 temp/otutab_trim3.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	# 复制最终版OTU表到结果目录
	cp temp/otutab_trim3.txt result/otutab.txt
	# 转换为biom格式
	biom convert -i result/otutab.txt -o result/otutab.biom --table-type="OTU table" --to-json
	# 统计OTU表
	biom summarize-table -i result/otutab.biom > result/otutab.biom.sum
	head -n 30 result/otutab.biom.sum

## OTU表抽样标准化

otutab_norm: otutab_filter
	touch $@
	# 依据最小样本量，设置OTU表标准化的阈值，如我们看到最小样品数据量为3.1万，可以抽样至3万
	usearch10 -otutab_norm result/otutab.txt -sample_size ${sample_size} -output result/otutab_norm.txt 
	usearch10 -otutab_stats result/otutab_norm.txt -output result/otutab_norm.txt.stat
	echo -ne "OtuNorm\t${sample_size}\t" >> ${log_otutable}
	head -n3 result/otutab_norm.txt.stat|awk '{print $$1}'|tr '\n' '\t'|sed 's/\t$$/\n/' >> ${log_otutable}
	cat ${log_otutable}
	# 转换为biom格式
	biom convert -i result/otutab_norm.txt -o result/otutab_norm.biom --table-type="OTU table" --to-json


## 1.12 物种注释 Assign taxonomy

tax_assign: otutab_norm
	touch $@
	usearch10 -sintax result/otu.fa \
		-db ${sintax_db} -sintax_cutoff ${sintax_cutoff} -strand both \
		-tabbedout temp/otu.fa.tax -threads ${p}


## 1.13 物种分类汇总 Taxonomy summary

tax_sum: tax_assign
	touch $@
	mkdir -p result/tax
	# 未分类的添加末注释标记，否则汇总时报错
	sed -i 's/\t$$/\td:Unassigned/' temp/otu.fa.tax
	# 按门、纲、目、科、属水平分类汇总
	for i in p c o f g;do \
		usearch10 -sintax_summary temp/otu.fa.tax -otutabin result/otutab_norm.txt -rank $${i} \
			-output result/tax/sum_$${i}.txt; \
	done
	# 删除Taxonomy中异常字符如() " - /
	sed -i 's/(//g;s/)//g;s/\"//g;s/\/Chloroplast//g;s/\-/_/g;s/\//_/' result/tax/sum_*.txt
	# 格式化物种注释：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
	cut -f 1,4 temp/otu.fa.tax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > result/taxonomy_2.txt
	# 生成物种表格：注意OTU中会有末知为空白，补齐分类未知新物种为Unassigned
	awk 'BEGIN{OFS=FS="\t"} {delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned"; split($$2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} print $$1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' result/taxonomy_2.txt | sed '1 i #OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > result/taxonomy_8.txt
	# 去除#号和空格，会引起读取表格分列错误
	sed -i 's/#//g;s/ //g' result/taxonomy_8.txt
	# 添加物种注释
	biom add-metadata -i result/otutab.biom --observation-metadata-fp result/taxonomy_2.txt -o result/otutab_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
	# 添加物种注释
	biom add-metadata -i result/otutab_norm.biom --observation-metadata-fp result/taxonomy_2.txt -o result/otutab_norm_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy
	# 制作门+变形菌纲混合比例文件
	cat <(grep -v 'Proteobacteria' result/tax/sum_p.txt) <(grep 'proteobacteria' result/tax/sum_c.txt) > result/tax/sum_pc.txt


# 1.14 多序列比对和进化树 Multiply alignment and make_phylogeny

tree_make: tax_sum
	touch $@
	# Reference based alignment
	#align_seqs.py -i result/otu.fa -t /mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta -o temp/aligned/
	# clustalo+qiime1.9.1
	clustalo -i result/otu.fa -o temp/otu_align.fa --seqtype=DNA --full --force --threads=${p}
	#此步用于pynast去除非目标区域，不适合clustalo，可能是导致unifrac结果不好的原因
	#filter_alignment.py -i temp/otu_align.fa  -o temp/
	#make_phylogeny.py -i temp/otu_align_pfiltered.fasta -o result/otu.tree
	make_phylogeny.py -i temp/otu_align.fa -o result/otu.tree

# 1.15 Alpha多样性指数计算 Calculate alpha diversity index

alpha_calc: tree_make
	touch $@
	mkdir -p result/alpha
	# 计算14种alpha多样性指数
	usearch10 -alpha_div result/otutab_norm.txt -output result/alpha/index.txt 
	# 稀释曲线：取1%-100%的序列中OTUs数量 Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)
	# method fast / with_replacement / without_replacement ref: https://drive5.com/usearch/manual/cmd_otutab_subsample.html
	usearch10 -alpha_div_rare result/otutab_norm.txt -output result/alpha/rare.txt -method ${rare_method}

# 1.16 Beta多样性进化树和距离矩阵计算 Beta diversity tree and distance matrix

beta_calc: alpha_calc
	touch $@
	# 计算距离矩阵，有多种方法结果有多个文件，需要目录
	mkdir -p result/beta/
	# 基于OTU构建进化树 Make OTU tree
ifeq (${tree_method}, usearch10)
	# usearch10 culster_agg建树+beta_div计算矩阵
	usearch10 -cluster_agg result/otu.fa -treeout result/otu_usearch.tree
	# ---Fatal error--- ../calcdistmxu.cpp(32) assert failed: QueryUniqueWordCount > 0 致信作者
	# 生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
	usearch10 -beta_div result/otutab.txt -tree result/otu_usearch.tree -filename_prefix result/beta/
	# ---Fatal error--- 1(91), expected ')', got '0.993' 它只依赖于cluster_agg的结果，fasttree结果不可用
else ifeq (${tree_method}, qiime)
	# clustalo + qiime1.9.1
	# 转换txt为biom才可以用qiime分析
	biom convert -i result/otutab_norm.txt -o result/otutab_norm.biom --table-type="OTU table" --to-json
	# 计算4种距离矩阵 http://qiime.org/scripts/beta_diversity.html -s显示矩阵列表有34种距离可选
	beta_diversity.py -i result/otutab_norm.biom -o result/beta/ -t result/otu.tree -m ${dis_method}
	# 删除文件名中多余字符，以方法.txt为文件名
	rename 's/_otutab_norm//' result/beta/*.txt
else ifeq (${tree_method}, mafft)
	# 方法3：mafft+fasttree
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in usearch10 or qiime, mafft") 
endif


## 1.17 有参比对Greengenes用于picurst, bugbase分析

otutab_gg: 
	touch $@
	# 根据方法选择usearch10/vsearch比对至gg13_5数据
ifeq (${map_method}, usearch10)
	usearch10 -otutab temp/filtered.fa -otus ${otutab_gg} -id ${map_identify} \
		-otutabout result/otutab_gg.txt -threads ${p}
else ifeq (${map_method}, vsearch)
	time vsearch --usearch_global temp/filtered.fa --db ${otutab_gg} --id ${map_identify} \
		--otutabout result/otutab_gg.txt --threads ${p}
else
	# 其它：没有提供正确的方法名称，报错提示
	$(error "Please select the right method: one of in usearch10 or vsearch") 
endif
	# 统计OTU表
	usearch10 -otutab_stats result/otutab_gg.txt -output result/otutab_gg.stat
	cat result/otutab_gg.stat


## 1.18 生成STAMP使用格式

stamp_input: tax_sum
	#touch $@
	# Data does not form a strick hierarchy. Child Unassigned has multiple parents. 不允许在各级别有重名，可直接使用sum_tax的各级别给STAMP使用
	#awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$$1]=$$0} NR>FNR {print a[$$1],$$0}' result/taxonomy_8.txt result/otutab_norm.txt | cut -f 2- | sed '1 s/#OTU ID/OTU/' | grep -v -P "Unassigned\tUnassigned" | grep -v '' > result/otutab_stamp.spf


## 1.19 test按丰度过滤OTU表计算beta

tree_make2: 
	touch $@
	# Reference based alignment
	align_seqs.py -i result/otu.fa -t /mnt/bai/public/ref/gg_13_8_otus/rep_set_aligned/97_otus.fasta -o temp/aligned/
	filter_alignment.py -i temp/aligned/otu_aligned.fasta -o temp/aligned/  # rep_seqs_align_pfiltered.fa, only very short conserved region saved
	make_phylogeny.py -i temp/aligned/otu_aligned_pfiltered.fasta -o temp/otu.tree # generate tree by FastTree

beta_calc2: tree_make2
	touch $@
	# OTU表筛选，根据具体组筛选
	rm -rf temp/beta/
	filter_otus_by_group_median.sh -i result/otutab_norm.txt \
		-d ${Dc_design} -A ${Dc_group_name} -B ${Dc_group_list} -t ${abundance_thre} \
		-o temp/otutab_norm.txt
	biom convert -i temp/otutab_norm.txt -o temp/otutab_norm.biom --table-type="OTU table" --to-json
	beta_diversity.py -i temp/otutab_norm.biom -o temp/beta/ -t temp/otu.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
	rename 's/_otutab_norm//' temp/beta/*.txt

beta_pcoa2: beta_calc2
	touch $@
	rm -f result/beta/*.p??
	beta_pcoa.sh -i temp/beta/ -m ${bp_method} \
		-d ${bp_design} -A ${bp_group_name} -B ${bp_group_list} -E ${bp_ellipse} \
		-c ${bp_compare} \
		-o ${bp_output} -h ${bp_height} -w ${bp_width}
# 过滤后整体PCoA样式不变，但解析率增加

# 2. 统计绘图 Statistics and plot

	# 绘图所需Shell脚本位于script目录中script目录下，会按参数生成R脚本于工作目录中的script下

plot_test:
	echo ${g1_list}

## 2.1 Alpha多样性指数箱线图 Alpha index in boxplot

alpha_boxplot: 
	touch $@
	rm -f result/alpha/*.p??
	alpha_boxplot.sh -i ${ab_input} -m ${ab_method} \
		-d ${ab_design} -A ${ab_group_name} -B ${ab_group_list} \
		-o ${ab_output} -h ${ab_height} -w ${ab_width}


## 2.2 Alpha丰富度稀释曲线 Alpha rarefracation curve

alpha_rare: alpha_boxplot
	touch $@
	alpha_rare.sh -i ${ar_input} \
		-d ${ar_design} -A ${ar_group_name} -B ${ar_group_list} \
		-o ${ar_output} -h ${ar_height} -w ${ar_width}


## 2.3 beta_pcoa 主坐标轴分析距离矩阵 PCoA of distance matrix

beta_pcoa: alpha_rare
	touch $@
	# 按组过滤OTUs再计算PCoA+CPCoA
	rm -rf temp/beta/
	filter_otus_by_group_median.sh -i result/otutab_norm.txt \
		-d ${Dc_design} -A ${Dc_group_name} -B ${Dc_group_list} -t ${abundance_thre} \
		-o temp/otutab_norm.txt
	biom convert -i temp/otutab_norm.txt -o temp/otutab_norm.biom --table-type="OTU table" --to-json
	beta_diversity.py -i temp/otutab_norm.biom -o temp/beta/ -t temp/otu.tree -m bray_curtis,weighted_unifrac,unweighted_unifrac
	rename 's/_otutab_norm//' temp/beta/*.txt
	rm -f result/beta/*.p??
	beta_pcoa.sh -i ${bp_input} -m ${bp_method} \
		-d ${bp_design} -A ${bp_group_name} -B ${bp_group_list} -E ${bp_ellipse} \
		-c ${bp_compare} \
		-o ${bp_output} -h ${bp_height} -w ${bp_width}


## 2.4 beta_cpcoa 限制性主坐标轴分析: OTU表基于bray距离和CCA  CCA of bray distance matrix

beta_cpcoa: beta_pcoa
	touch $@
	beta_cpcoa.sh -i ${bc_input} -m ${bc_method} \
		-d ${bc_design} -A ${bc_group_name} -B ${bc_group_list} -E ${bc_ellipse} \
		-o ${bc_output} -h ${bc_height} -w ${bc_width}


## 2.5 tax_stackplot 样品和组分类学各级别的堆叠柱状图 Stackplot showing taxonomy in each level

tax_stackplot: beta_cpcoa
	touch $@
	rm -f result/tax/*.p??
	tax_stackplot.sh -i ${ts_input} -m ${ts_level} -n ${ts_number} \
		-d ${ts_design} -A ${ts_group_name} -B ${ts_group_list} -O ${ts_order} \
		-o ${ts_output} -h ${ts_height} -w ${ts_width}


## 2.6 DA_compare 组间差异比较 edgeR or wilcox

DA_compare: tax_stackplot
	touch $@
	rm -fr ${Dc_output}
	mkdir -p ${Dc_output}
	compare.sh -i ${Dc_input} -c ${Dc_compare} -m ${Dc_method} \
		-p ${Dc_pvalue} -q ${Dc_FDR} -F ${Dc_FC} -t ${abundance_thre} \
		-d ${Dc_design} -A ${Dc_group_name} -B ${Dc_group_list} \
		-o ${Dc_output} -C ${Dc_group_name2}
### 计算门纲目科属水平秩合检验差异，绘图维恩图，并对维恩图中各部分进行重叠
DA_compare_tax: tax_stackplot
	touch $@
	for i in p pc c o f g;do \
	rm -fr result/compare_$${i}; \
	mkdir -p result/compare_$${i}; \
	compare.sh -i result/tax/sum_$${i}.txt -c ${Dct_compare} -m ${Dct_method} \
		-p ${Dct_pvalue} -q ${Dct_FDR} -F ${Dct_FC} -t ${Dct_thre} \
		-d ${Dct_design} -A ${Dct_group_name} -B ${Dct_group_list} \
		-o result/compare_$${i}/; \
	batch_venn.pl -i ${doc}/venn.txt -d result/compare_$${i}/diff.list; \
	batch2.pl -i result/compare_$${i}/diff.list.venn'*'.xls -d result/compare_$${i}/database.txt -o result/compare_$${i}/ -p vennNumAnno.pl; \
	done


# 2.7 plot_volcano (可选) 基于差异OTU表绘制火山图
plot_volcano: DA_compare
	touch $@
	# 指定文件绘制单组差异
	# plot_volcano.sh -i result/compare/ACT2KO-Col_all.txt -o result/compare/ACT2KO-Col
	# awk调用比较列表批量绘制文件，grep -v删空行
	awk 'BEGIN{OFS=FS="\t"}{system("plot_volcano.sh -i result/compare/"$$1"-"$$2"_all.txt -o result/compare/"$$1"-"$$2);}' \
		<(grep -v '^$$' ${Dc_compare})


# 2.8 差异OTU绘制热 图
plot_heatmap: plot_volcano
	touch $@
	# 指定文件绘制单个图
	# plot_heatmap.sh -i result/compare/V3703HnCp6-ZH11HnCp6_sig.txt -o result/compare/V3703HnCp6-ZH11HnCp6 -w 5 -h 7
	# awk调用批量绘制文件，grep -v删空行
	awk 'BEGIN{OFS=FS="\t"}{system("plot_heatmap.sh -i result/compare/"$$1"-"$$2"_sig.txt \
		-o result/compare/"$$1"-"$$2" -w ${ph_width} -h ${ph_height}");}' \
		<(grep -v '^$$' ${Dc_compare})

plot_heatmap_tax: DA_compare_tax
	touch $@
	awk 'BEGIN{OFS=FS="\t"}{system("plot_heatmap.sh -i result/compare_${ph_tax}/"$$1"-"$$2"_sig.txt \
		-o result/compare_${ph_tax}/"$$1"-"$$2" -w ${ph_width} -h ${ph_height}");}' \
		<(grep -v '^$$' ${Dc_compare})


# 2.9 差异OTU绘制曼哈顿图

plot_manhattan: plot_heatmap
	touch $@
	#f输入文件，x为X轴，y为Y轴，g为简化和物种，s为上下调，p为丰度
	#sp_manhattan2.sh -f result_k1-c/otu_INDvsTEJ -x otu -y PValue -g tax -s level -p A_mean
	#sp_manhattan2.sh -f result/compare/HTEJ-HIND_all.txt -x HTEJ_HIND -y PValue -g Phylum -s level -p logCPM
	# 陈同程序有待优化，先用我的脚本重画
	#	awk 'BEGIN{OFS=FS="\t"}{system("sp_manhattan2.sh -f result/compare/"$$1"-"$$2"_all.txt \
	#		-x "$$1"_"$$2" -y PValue -g Phylum -s level -p logCPM");}' \
	#		<(grep -v '^$$' ${Dc_compare})
	# plot_manhattan.sh -i result/compare/LTEJ-LIND_all.txt
	awk 'BEGIN{OFS=FS="\t"}{system("plot_manhattan.sh -i result/compare/"$$1"-"$$2"_all.txt");}' \
		<(grep -v '^$$' ${Dc_compare})


# 2.10 单个差异OTU绘制箱线图

plot_boxplot: 
	touch $@
	mkdir -p ${pb_output}
	plot_boxplot_ggpubr.sh -i ${pb_input} -d ${pb_design} -A ${g1} -B ${bp_group_list} \
		-m ${pb_list} -t ${pb_trans} -o ${pb_output} -n ${pb_norm}


# 2.11 维恩图

plot_venn: plot_manhattan
	touch $@
	# mkdir -p result/venn
	# 绘制维恩图ABCDE，获得比较列表AB
	batch_venn.pl -i ${venn} -d result/compare/diff.list
	# 注释列表为OTU对应描述，目前添加培养注释
	rm -f result/compare/*.xls.xls
	batch2.pl -i 'result/compare/diff.list.venn*.xls' -d ${venn_anno} -o result/venn/ -p vennNumAnno.pl


# 2.11 Upsetview图



# 3. 高级分析 Advanced analysis

## 3.1 lefse 多种差异物种特征分析

## 3.2 picurst GG宏基因组预测
picrust_calc: 
	touch $@
	mkdir -p picurst
	biom convert -i result/otutab_gg.txt -o picurst/otutab.biom --table-type="OTU table" --to-json
	biom summarize-table -i picurst/otutab.biom > picurst/otutab.stat
	cat picurst/otutab.stat
	# 校正拷贝数
	normalize_by_copy_number.py -i picurst/otutab.biom -o picurst/otutab_norm.biom -c /db/picrust/16S_13_5_precalculated.tab.gz
	# 预测宏基因组KO表
	predict_metagenomes.py -i picurst/otutab_norm.biom -o picurst/ko.biom -c /db/picrust/ko_13_5_precalculated.tab.gz
	predict_metagenomes.py -f -i picurst/otutab_norm.biom -o picurst/ko.txt  -c /db/picrust/ko_13_5_precalculated.tab.gz
	# 按功能级别分类汇总, -c指输出类型，有KEGG_Pathways, COG_Category, RFAM三种，-l是级别，分4级，初始KO为4级，可全并为1-3级
	categorize_by_function.py -f -i picurst/ko.biom -c KEGG_Pathways -l 3 -o picurst/ko3.txt
	categorize_by_function.py -f -i picurst/ko.biom -c KEGG_Pathways -l 2 -o picurst/ko2.txt
	categorize_by_function.py -f -i picurst/ko.biom -c KEGG_Pathways -l 1 -o picurst/ko1.txt

## 3.3 faprotax 元素循环
	
	# 需要有物种注释的biom文件作为输入
faprotax_calc:
	touch $@
	# 原始OTU表没有标准化，后期结果再标准化可能不合理。采用抽样标准化的OTU表
	mkdir -p result/faprotax
	/usr/bin/python2.7 /mnt/bai/yongxin/software/FAPROTAX_1.1/collapse_table.py -i result/otutab_norm_tax.biom -o result/faprotax/element_tab.txt -g /mnt/bai/yongxin/software/FAPROTAX_1.1/FAPROTAX.txt --collapse_by_metadata 'taxonomy' -v --force --out_report result/faprotax/report 

	# 绘制指定marker的箱线图
plot_fa_barplot: faprotax_calc
#	touch $@
	alpha_boxplot.sh -i result/faprotax/element_tab.txt -d ${Dc_design} -A ${Dc_group_name} -B ${Dc_group_list} \
		-m ${fapro_list} -t TRUE -o result/faprotax/ -n FALSE


## 3.4 bugbase GG表型预测

## 3.5 tax4fun Silva宏基因组预测

## 3.6 humman2 宏基因组代谢通路分析

## 3.7 network 网络分析

## 3.8 vegan 环境因子分析


## 3.9 culture 培养菌

	# 筛选每个OTUs在菌库中的相似度和覆盖度，挑选高丰度的绘制Graphlan图
culture: 
	mkdir -p result/39culture
	# 比对OTU至可培养菌blast数据库
	blastn -query result/otu.fa -db ${culture_db} -out temp/culture_otu.blastn -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore' -num_alignments 1 -evalue 1 -num_threads ${p}
	# 添加blastn结果表头，最主要前4列：OTUID，培养菌ID，相似度，覆盖度
	sed -i '1 i OTUID\tsseqid\tpident\tqcovs\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' temp/culture_otu.blastn
	# 统计可培养菌所占种类和丰度比例
	echo -ne "Total OTUs\t" > result/39culture/summary.txt
	grep '>' -c result/otu.fa >> result/39culture/summary.txt
	echo -ne "Cultured OTUs\t" >> result/39culture/summary.txt
	awk '$$3>=97 && $$4>=99' temp/culture_otu.blastn|wc -l >> result/39culture/summary.txt
	# 计算平均丰度
	otutab_mean.sh -i result/otutab.txt -o temp/otutab.mean
	# 添加丰度至culture
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$$1]=$$2} NR>FNR {print $$0,a[$$1]}' temp/otutab.mean temp/culture_otu.blastn | cut -f 1-4,14 > temp/temp
	# 添加物种注释
	awk 'BEGIN{OFS=FS="\t"} NR==FNR {a[$$1]=$$4} NR>FNR {print $$0,a[$$2]}' ${culture_db}.tax temp/temp | sed '1 s/$$/Taxonomy/' > result/39culture/otu.txt
	echo -ne "Cultured abundance\t" >> result/39culture/summary.txt
	awk '$$3>=97 && $$4>=99' result/39culture/otu.txt | awk '{a=a+$$5} END {print a}' >> result/39culture/summary.txt
	cat result/39culture/summary.txt



# 4. 输出结果报告 Write Rmarkdown HTML report

rmd: tax_stackplot
	report_16S.pl -g ${g1} -b ${version} -m ${compare_method} -d ${design} -c ${compare} -v ${venn} -a ${abundance_thre} -F ${FC} -p ${pvalue} -q ${FDR} # -D ${g2_list} -F ${g3_list} -l ${library} -S ${elite_report} -s ${summary}  -t ${tern}
	ln -sf ${wd}/${version}/ /var/www/html/report/16Sv2/${version}
	rm -f ${version}/${version}
