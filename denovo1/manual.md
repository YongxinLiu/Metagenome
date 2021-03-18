	# 宏基因组分析流程第一版 —— 操作脚本
	# Metagenome pipeline version 1 —— Manual script

	# Biocloud中重置环境变量(可选)
	PATH=/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
	# 调置Perl模块
	PERL5LIB=/conda/envs/kraken2/lib/site_perl/5.26.2/x86_64-linux-thread-multi:/conda/envs/kraken2/lib/site_perl/5.26.2:/conda/envs/kraken2/lib/5.26.2/x86_64-linux-thread-multi:~/miniconda2/envs/kraken2/lib/5.26.2

	# 加载目标环境变量
	source /conda/bin/activate
	# 启动宏基因组通用分析环境
	conda activate meta

# 1. 有参分析流程 Reference-based pipeline

	# 0. 准备工作 Preparation

	# 设置工作目录 Set work directory
	wd=medicago/lyr4/
	mkdir -p ~/${wd}
	cd ~/$wd
	
	# 准备流程 Prepare makefile
	ln -s /home/meta/soft/Metagenome/denovo1/parameter.md makefile
	cp /home/meta/soft/Metagenome/denovo1/manual.md manual.sh
	
	# 建立初始工作目录 Create initial working directory
	# 代码预览、执行、保存
	make -n -B init
	make init
	echo "#" `date` > result/pipeline.sh
	make -n -B init >> result/pipeline.sh

	# 准备原始数据 sequencing raw data (多样本合并和统计见附录1)
	# 链接数据至工作目录
	ln -s /mnt/m2/data/meta/$wd/*.gz seq/

	# 准备实验设计上传到result目录，至少有两列样本名和组名 Experiment design
	# 方法1. 数据来源处复制实验设计
	cp /mnt/m2/data/meta/$wd/metadata.txt result/metadata.txt
	# 方法2. 复制实验设计模板并手动填写
	cp /home/meta/soft/Metagenome/denovo1/result/design.txt result/metadata.txt
	# 方法3. 从样本名中提取，并手动补充
	ls seq/*_1.fq.gz|cut -f 2 -d '/'|cut -f 1 -d '_'|awk '{print $1"\t"$1}'|sed '1 i SampleID\tGroupID' > result/metadata.txt

	# 检查实验设计ID唯一性，有结果则为重复
	cut -f1 result/metadata.txt|sort|uniq -d

## 1.1. 质控并移除宿主 Quality control & Remove host

	### 1.1.1 (可选)质量评估原始数据
	# 依赖软件版本：质量评估、评估报告汇总
	fastqc -v # v0.11.9
	multiqc --version # 1.8
	# 预览代码
	make -n -B qa
	# 计算：240G,72p,30m
	time make qa
	# 查看result/multiqc_report.html，重点检查数据量分布，重复率，接头含量

	### 1.1.2 KneadData移除低质量和宿主
	echo -e "\n#" `date` >> result/pipeline.sh
	# 依赖软件版本：并行、移除低质量和宿主
	echo -n "kneaddata --version # " >> result/pipeline.sh
	kneaddata --version >> result/pipeline.sh 2>&1
	# 预览代码
	make -n -B qc
	# 计算：150G,24p,20h; 240G,72p,12h; ; 126G,24x3p,2.5h,99G; 异常中断处理见附录2. KneadData中断
	time make qc
	# 结果见 result/11kneaddata_stat.txt 高质量、非宿主比例
	# 保存运行代码
	make -n -B qc >> result/pipeline.sh

	### 1.1.3 提取上传的Clean数据
	# 预览代码
	make -n -B qa2
	# 	# 240G,72p,1h，按GSA标准整理上传数据，见submit目录
	make qa2
	make -n -B qa2 >> result/pipeline.sh

	### 1.1.4 kraken2质控
	conda activate kraken2 # bicloud
	kraken2 --version # 2.1.1
	# 采用kraken完整数据库注释, 370G,24p,10h; 330G,24x3p,1h;
	make -n -B kraken2_qc
	make kraken2_qc
	make -n -B qa2 >> result/pipeline.sh

## 1.2. 物种和功能组成定量 humman2

	## 1.2.1 humman2输入文件准备：双端文件cat连接
	make humann2_concat

	# 启动humann2环境
	conda activate humann2
	humann2 --version # v2.8.1

	# 检查数据库位置
	humann2_config --print
	
	## 1.2.2 humman2计算，包括metaphlan2
	make humann2
	# 4.7 Tb水稻数据，8X12线程运行2.5 Days
	#  X Tb Data, 3 x 24，运行 xx

	## 1.2.3 功能组成整理 humman2_sum
	make humann2_sum
	# 结果见 result/12humann2目录，有功能通路及物种组成表uniref.tsv、标准化表uniref_relab.tsv，以及拆分功能表unstratified和功能物种对应表stratified
	
	# humann2转为KEGG
	humann2_regroup_table -i temp/humann2/genefamilies.tsv \
	  -g uniref90_ko -o temp/humann2/ko.tsv

## 1.3. 整理物种组成表和基本绘图 Summary metaphlan2 and plot

	# 结果见result/13metaphlan2目录

	### 1.3.1 整理物种组成表 Summary metaphlan2
	make metaphaln2_sum
	# 结果taxonomy*文件，有多级物种表.tsv、株水平表.spf用于STAMP分析，以及聚类热图_heatmap_top.pdf观察分组情况

	### 1.3.2 GraPhlAn图
	make metaphaln2_graphlan
	# 结果taxonomy_graphlan*.pdf，包括主图、图例和注释文本3个文件

	### 1.3.3 物种组成LEfSe差异分析(可选)
	# 依赖实验设计，分组比较末确定时，可选跳过此步
	make metaphaln2_lefse


## 1.4. kraken2物种组成

	### 1.4.1 基于NCBI完整基因组数据库的k-mer物种注释
	# Taxonomy assign by k-mer and based on NCBI database
	kraken2 --version # Kraken version 2.0.9-beta
	make -n -B kraken2_read
	# 500G,3x24p,3h; 100G,3x24p,11m
	make kraken2_read
	# record pipeline
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B kraken2_read >> result/pipeline.sh

	### 1.4.2 合并为矩阵 merge into matrix
	make kraken2_read_sum



# 2. 无参分析流程 De novo assemble pipeline


## 2.1. Assemble 组装

	# 启动组装环境
	conda activate meta
	megahit -v # v1.2.9

	# 方法1. 大项目推荐
	### 2.2.1 基于qc质控序列单样本拼接(大项目)
	make assemble_single
	# 66个样168G压缩数据分别装为5h 计算过程日志见 temp/22megahit/megahit.log

	# 方法2. 小项目推荐
	### 2.2.2 基于qc质控后序列拼接
	# metahit拼接所有样本, 496G,72p,3d; 38G,48p,16h; 360G,72p,58h; 297G,72p,8h; 142G,48p,7d；1.5d,243G
	make -n -B megahit_all
	make megahit_all
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B megahit_all >> result/pipeline.sh

	# 问题. 数据量大无法拼接
	# 126 - Too many vertices in the unitig graph (8403694648 >= 4294967294), you may increase the kmer size to remove tons
	# 需要增加k-mer，如21增加为29

	### 2.2.3 megahit_all_quast评估
	# 3G,15m; 2G,7m
	make -n -B megahit_all_quast
	make megahit_all_quast

	### 2.2.5 (可选)kraken2对contig物种注释
	# 10M,9G,48m; 2G,
	conda activate kraken2
	make -n -B kraken2_contig
	make kraken2_contig
	conda deactivate
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B kraken2_contig >> result/pipeline.sh

## 2.3. Genome annotation 基因组注释

	### 2.3.1 对单样品组装的每个contig文件基因注释(大数据可选)
	make prodigal_single

	### 2.3.2 对合并组装的单个contig文件基因注释(小样本可选)
	# 3G,10h; 9G,20h；2G,
	prodigal -v # V2.6.3: February, 2016
	make -n -B prodigal_all
	make prodigal_all
	
	# 并行处理，加速基因预测
	make -n -B prodigal_all_split
	# 2.5h,473M
	make prodigal_all_split
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B prodigal_all_split >> result/pipeline.sh

	### 2.3.3 构建非冗余基集 Non-redundancy gene set(大数据可选)

	# 90%覆盖度，95%相似度下，拼接结果再聚类基本不变少，如7736减少为7729
	# 5,514,236聚类为4,879,276，72线程155m，8756m；
	# 15,889,664聚类为13,924,854，48线程1651m(1d)，71427m；
	cd-hit-est -v # 4.8.1 (built on Oct 26 2019)
	make -n -B NRgene
	make NRgene
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B NRgene >> result/pipeline.sh

	### 2.3.4 基因定量 salmon genes
	conda activate metawrap1.3 # 中salmon仅为0.13.1
	salmon -v # 1.3.0
	# 14M，
	make -n -B salmon_gene
	make salmon_gene
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B salmon_gene >> result/pipeline.sh

	### 2.3.5 基因物种注释 kraken2 annotate gene

	make kraken2_gene


## 2.4 功能数据库注释

### 2.4.2 eggNOG
	# biocloud上: conda activate biobakery
	emapper.py --version # 2.0.1
	diamond --version # 2.0.6
	# 预览代码
	make -n -B eggnog
	# eggnog注释
	# 15M 72p, 6h
	# The host system is detected to have 1082 GB of RAM. It is recommended to use this parameter for better performance: -c1
	# 默认diamond比对，采用：--more-sensitive  -e 0.001 --top 3 --query-cover 0 --subject-cover 0
	make eggnog
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B eggnog >> result/pipeline.sh

	# 9m，26G
	make -n -B eggnog_sum
	make eggnog_sum
	# 注释基因数量 3781759/4879276=77.5; 9036098/13,924,854=64.9; 3025768/3886015=77.8; 7189141/9729931=73.9
	# 注释KO基因数量 2335006/4879276=47.8; 5836548/13,924,854=41.9; 1955582/3886015=50.3; 4545144/9729931=46.7
	# 合并，15Mx76,75m; 4Mx25,3m; 7Mx25,11m
	echo -e "\n#" `date` >> result/pipeline.sh
	make -n -B eggnog_sum >> result/pipeline.sh

### 2.4.2 KEGG

	# 比对kegg 2020 数据库
	make kegg
	# 汇总为result/24kegg/kotab.count/tpm两个表，分别为原始count和rpm
	make kegg_sum

### 2.4.3 CAZy

	# blast比对到CAZy
	make -n -B dbcan2
	make -n -B dbcan2_sum

### 2.4.4 CARD
	
	# 启动rgi环境
	conda activate rgi
	
	make -n -B card

# 3 分箱

	conda activate metawrap1.3
	metawrap -v # 1.3.2
	megahit -v # v1.1.3
	metaspades.py -v # v3.13.0

## 3.1 混合分箱

## 3.2 分批次分箱

## 3.3 单样本分箱

  # 样本名
  i=R108B3R9
  # 线程数
  p=48
  # 任务数
  j=3
  # 定义完整度和污染率的阈值(50, 5; Finn NBT 2020;50, 10, Bowers NBT 2017)
  c=50
  x=10
  
  for i in `cut -f 1 result/metadata.txt |tail -n+8`; do
  # 解压改名，满足metawrap要求
  gunzip -c submit/${i}_1.fq.gz > submit/${i}_1.fastq
  gunzip -c submit/${i}_2.fq.gz > submit/${i}_2.fastq
  # 单样本分箱2h
  mkdir -p temp/binning
  metawrap binning \
      -o temp/binning/${i}\
      -t ${p} \
      -a temp/22megahit/${i}/final.contigs.fa \
      --metabat2 --maxbin2 --concoct \
      submit/${i}_*.fastq
  # 单样本分箱提纯5h
  mkdir -p temp/bin_refinement
  metawrap bin_refinement \
    -o temp/bin_refinement/${i} -t ${p} \
    -A temp/binning/${i}/metabat2_bins/ \
    -B temp/binning/${i}/maxbin2_bins/ \
    -C temp/binning/${i}/concoct_bins/ \
    -c ${c} -x ${x}
  # 统计中、高质量bin的数量
  mkdir -p result/32bin_single
  echo -e $i"\t\c" >> result/32bin_single/medium.txt
  tail -n+2 temp/bin_refinement/${i}/metawrap_50_10_bins.stats| awk '$2>=50 && $3<10'| wc -l >> result/32bin_single/medium.txt
  echo -e $i"\t\c" >> result/32bin_single/high.txt
  tail -n+2 temp/bin_refinement/${i}/metawrap_50_10_bins.stats| awk '$2>90 && $3<5'| wc -l >> result/32bin_single/high.txt    
  done

  # 显示结果
  cat result/32bin_single/*.txt

## 3.3 Drep去冗余


## 3.4 物种注释和进化



# 附录1. 数据库准备工作

## 1.1 KEGG数据库注释整理
cd ~/github/Metagenome/denovo1
mkdir script
mkdir kegg
# https://www.kegg.jp/kegg-bin/get_htext?ko00001.keg 下载htext于kegg目录
cp ~/bin/kegg_ko00001_htext2tsv.pl script/
kegg_ko00001_htext2tsv.pl -i kegg/ko00001.keg -o kegg/ko00001.tsv
# geneID-KO-Description
cd /db/kegg
# 提取基因-KO-描述列 shell命令存在Klebsiella的bug问题
#grep '>' /db/kegg/kegg_all_clean.fa |sed 's/>//' > kegg_all_clean.title
#sed 's/  /; K/;s/; K/\t/g' kegg_all_clean.title|cut -f 1,3|sed 's/\t/\tK/;s/ /\t/' > /db/kegg/kegg_gene_KO.list
format_kegg76_title.pl -i kegg_all_clean.fa -o kegg_gene_ko_description.txt
# 统计KO描述对应表
cut -f 2- kegg_gene_ko_description.txt|sort|uniq|sed '1 i ID\tKDescription'>ko_description.txt
# 统计基因和KO数量
wc -l kegg_gene_KO.list # 7,857,187 genes in KEGG
cut -f 2 kegg_gene_KO.list|sort|uniq|wc -l # 18,648 KO in KEGG

# KEGG比对加KO结果整理
# 提取基因ID(Name)和KEGG基因ID(KgeneID)
cut -f 1,2 temp/24kegg/gene_diamond.f6|uniq | sed '1 i Name\tKgeneID' > temp/24kegg/gene_kegg.list
# 添加KO编号(KO)和描述(Kdescription)
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print $0,a[$2]}' /db/kegg/kegg_gene_ko_description.txt temp/24kegg/gene_kegg.list > temp/24kegg/gene_ko.list
# 基因丰度矩阵末尾添加对应KO编号，没注释的直接删除，可选注释为unclassified
# awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3} NR>FNR{print $0,a[$1]}' temp/24kegg/gene_ko.list result/23salmon_gene/gene.count | sed 's/\t$/\tunclassified/' > temp/24kegg/gene_ko.count
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$3} NR>FNR{print $0,a[$1]}' temp/24kegg/gene_ko.list result/23salmon_gene/gene.count | sed '/\t$/d' > temp/24kegg/gene_ko.count
# 检查注释前后基因数量
wc -l temp/24kegg/gene_ko.count
wc -l result/23salmon_gene/gene.count 
# 合并基因表为KO表，输出count值和tpm值
mkdir -p result/24kegg
Rscript ~/github/Metagenome/denovo1/script/mat_gene2ko.R -i temp/24kegg/gene_ko.count -o result/24kegg/kotab
# 结果中添加KO注释, 个别KO没有注释
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' /db/kegg/ko_description.txt result/24kegg/kotab.count > result/24kegg/kotab.count.anno

# 查看文件尾有异常
tail result/24kegg/kotab* # Klebsiella和Unknown有小数问题？

# 整理多级合并，4-3-2-1



## 1.2 COG/EggNOG

http://eggnogdb.embl.de 4.5.1 Nov 2016

# 下载软件和数据库 Downloads
## 软件
cd ~/software
wget https://github.com/jhcepas/eggnog-mapper/archive/1.0.3.tar.gz
tar xvzf 1.0.3.tar.gz
cd eggnog-mapper-1.0.3


# trimmed蛋白库
wget http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.trimmed_algs.tar.gz
tar xvzf NOG.trimmed_algs.tar.gz

## COG描述
http://www.sbg.bio.ic.ac.uk/~phunkee/html/old/COG_classes.html 保存为 /mnt/bai/yongxin/data/db/eggnog/COG_one_letter_code_descriptions.txt
手动整理为/mnt/zhou/yongxin/db/eggnog/COG_one_letter_code_descriptions.tsv


## 1.3 CAZyme碳水化合物数据库

# 安装结构域预测工具
conda install hmmer

# 下载相关dbCAN软件和数据库
cd ~/data/db
mkdir dbCAN2 && cd dbCAN2
# 蛋白库、描述(400多类)和基因酶学编号
wget http://cys.bios.niu.edu/dbCAN2/download/Databases/CAZyDB.07312018.fa
wget http://cys.bios.niu.edu/dbCAN2/download/Databases/CAZyDB.07312018.fam-activities.txt
wget http://cys.bios.niu.edu/dbCAN2/download/Databases/CAZyDB.07312018.pr-with-ec.txt
# 建索引，不支持 --threads 9
diamond makedb --in CAZyDB.07312018.fa --db CAZyDB.07312018
## 蛋白与分类存在1对多问题，如>AWI89010.1|GT2|GT4，暂时只考虑其第一类
grep -v '#' CAZyDB.07312018.fam-activities.txt|sed 's/  //'|sed '1 i ID\tDescription' > fam_description.txt # cat -A|les

# diamond比对结果整理
diamond blastp --db /mnt/zhou/yongxin/db/dbCAN2/CAZyDB.07312018 --query temp/23prokka_all/mg.faa \
        --outfmt 6 --threads 9 --max-target-seqs 1 --quiet \
        --out temp/24dbcan2/gene_diamond.f6
mkdir -p result/24dbcan2
# 提取基因对应基因加族，同一基因存在1对多，只取第一个
cut -f 1,2 temp/24dbcan2/gene_diamond.f6 | uniq | sed 's/|/\t/g' | cut -f 1,3 | cut -f 1,2 -d '_' |sed '1 i Name\tKO' > temp/24dbcan2/gene_fam.list
# 基因丰度矩阵末尾添加对应KO编号，没注释的直接删除，可选注释为unclassified
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/24dbcan2/gene_fam.list result/23salmon_gene/gene.count | sed '/\t$/d' > temp/24dbcan2/gene_fam.count
wc -l result/23salmon_gene/gene.count
wc -l temp/24dbcan2/gene_fam.count
Rscript ~/github/Metagenome/denovo1/script/mat_gene2ko.R -i temp/24dbcan2/gene_fam.count -o result/24dbcan2/cazytab
# 结果中添加KO注释, 个别KO没有注释
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' /mnt/bai/yongxin/data/db/dbCAN2/fam_description.txt result/24dbcan2/cazytab.count > result/24dbcan2/cazytab.count.anno

# hmm注释，小数据更快，远同源更准，但不适合大数据库
cd ~/data/db/dbCAN2/Tools/run_dbcan
hmmsearch /mnt/bai/yongxin/data/db/dbCAN2/dbCAN-HMMdb-V7.txt temp/23prokka_all/mg.faa


## 1.4 ResFams 抗生素抗性数据库

cd ~/data/db/ResFams
wget -c http://dantaslab.wustl.edu/resfams/Resfams-proteins.tar.gz
tar xvzf Resfams-proteins.tar.gz
# 制作基因与resfam对应表
cd proteins/
for file in *.faa; do
  # 提取样品名
  name=${file%%.*}
  # echo $name
  # 将每个序列基因名后+文件名分类
  sed "s/ /\t$name\t/" $file > ${file}.mod
done
# ArmA.faa异常没有空格，需要单独匹配
awk '{if(/>/){print $0"\tArmA\t"}else{print $0}}' ArmA.faa > ArmA.faa.mod
cd ..
# 添加resfam注释
# 合并蛋白为单个文件，直接合并存在末换行的序列，需要添加换行，再去除空行
cat proteins/*.mod | sed 's/>/\n>/' | grep -v -P '^$'> Resfams-proteins.faa
# 提取基因与ResFam分类对应表
grep '>' Resfams-proteins.faa|cut -f 1,2|sed 's/>//'|sed '1 i ResGeneID\tResfam'|less>Resfams-proteins.id
# 添加ResFam按机制分类
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$2]=$3"\t"$8} NR>FNR{print $0,a[$2]}' resfam_metadata.tsv Resfams-proteins.id > Resfams-proteins_class.tsv




# 附录2. 常用问题及处理


## S2.1 多批数据合并基因集

	# 手动运行 23prokka_single
	cat temp/23prokka/*/mg.ffn temp/23prokka2/*/mg.ffn ../miniCore2/temp/23prokka/*/mg.ffn > temp/23prokka/mg.ffn # 19.7GB
	mkdir -p temp/23NRgene
	cp temp/23prokka/mg.ffn temp/23NRgene/
	echo -ne 'Prokka single genes\t' > result/gene.log
	grep -c '>' temp/23NRgene/mg.ffn >> result/gene.log




# 附录3. shell分析流程

# 一、宏基因组有参分析流程 metagenome reference-based pipeline (humann2)

# 系统要求 System: Linux Ubuntu 18.04 / CentOS 7.5
# 依赖软件 Sofware: KneadData、metaphlan2、humann2
# 运行前准备
# 1. 按7software, 8database目录中软件和数据库按课件说明安装，并添加环境变量
# 2. 学员U盘复制测序数据3metagenome目录至服务器~目录
# 3. Rstudio打开pipeline_ref_humann2.sh文件，Terminal中切换至工作目录

# 中文教程：https://mp.weixin.qq.com/s/XkfT5MAo96KgyyVaN_Fl7g
# 英文教程：https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Humann2)

# 上课演示Linux服务器：数据己经拷备到服务器/db/3metagenome目录，
# 用户登陆Rstudio网页版：192.168.1.107:8787(内网) 即可


# 1. 了解工作目录和文件

# 建立并进入工作目录
mkdir -p 3metagenome
cd 3metagenome
ln -s /db/3metagenome/* ./ # 准备工作流程

# 目录
tree # 显示文件结构
# 实验设计 design file
head -n3 doc/design.txt

# 文件说明
# pipeline_ref_humann2.sh 有参分析主流程

# seq/*.fq 原始测序数据，公司返回的测序结果，通常为一个样品一对fastq格式文件
# 如果测序数据是.gz结尾的压缩文件，使用gunzip解压，注意测序文件命名，结果可以是fastq/fq，左端可以_1/R1.fq
# gunzip seq/* # 如果压缩文件还需要解压
# 测序数据有12个样本，共24个文件。只取1%抽样用于测试。
head -n4 seq/p136C_1.fq
mkdir -p temp # 临时文件 temp directory for intermediate files



# 2. 序列质控和去宿主 Qaulity control and remove host contamination

# kneaddata是一套工作流，依赖trimmomatic进行质控和去接头；依赖bowtie2比对宿主，并筛选非宿主序列
# kneaddata -h # 显示帮助
# kneaddata_database # 查看可用数据库
# kneaddata_database --download human_genome bowtie2 ./ # 如下载人类基因组bowtie2索引至当前目录，并脚本

# 以p136C单样品合并为例
time kneaddata -i seq/p144C_1.fq -i seq/p144C_2.fq \
  -o temp/qc -v -t 4 --remove-intermediate-output \
  --trimmomatic /conda2/share/trimmomatic-0.36-3/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
  --bowtie2-options "--very-sensitive --dovetail" -db /db/bowtie2/Homo_sapiens   
# 单样品，质控17s
    
# 现实中是有一大堆样品，你可以逐个修改样品名运行，如果服务器性能允许可以并行加速分析
# parallel --citation # 打will cite以后不再提醒
time parallel -j 3 --xapply \
  'kneaddata -i {1} -i {2} \
  -o temp/qc -v -t 3 --remove-intermediate-output \
  --trimmomatic /conda2/share/trimmomatic-0.36-3/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
  --bowtie2-options "--very-sensitive --dovetail" -db /db/bowtie2/Homo_sapiens' \
 ::: seq/*_1.fq ::: seq/*_2.fq
# 32s，系统用时2m38s

	# 单样本分析代码
	time kneaddata -i 02seq/HnZH11R3_1.fq.gz -i 02seq/HnZH11R3_2.fq.gz \
		-o temp/11qc -v -t 8 --remove-intermediate-output \
		--trimmomatic /conda2/share/trimmomatic-0.36-3/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
		--bowtie2-options "--very-sensitive --dovetail" -db /db/rice/IndJap

	# 多样本并行
	parallel --xapply -j 6 \
		"kneaddata -i 02seq/{1}_1.fq.gz -i 02seq/{1}_2.fq.gz \
		-o temp/11qc -v -t 8 --remove-intermediate-output \
		--trimmomatic /conda2/share/trimmomatic-0.36-3/ --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:50' \
		--bowtie2-options '--very-sensitive --dovetail' -db /db/rice/IndJap" \
		::: `tail -n+2 01doc/design.txt | cut -f 1`
	
	# 结果汇总
	kneaddata_read_count_table --input 03temp/11qc --output 04result/11kneaddata_stat.txt
	cat 04result/11kneaddata_stat.txt

	# 批量质控和汇总

# 质控结果统计
kneaddata_read_count_table --input temp/qc --output temp/kneaddata_read_counts.txt
cat temp/kneaddata_read_counts.txt

# 合并质控后样品文件：有参宏基因组不考虑双端，将单双端和双端不完全序列共4个文件合并
mkdir -p temp/concat
for i in `tail -n+2 doc/design.txt | cut -f 1`;do \
  cat temp/qc/${i}*data_paired* > temp/concat/${i}.fq; done
ll temp/concat/*.fq # 查看样品数量和大小



# 3. 计算物种组成和代谢通路

# 物种组成调用metaphlan2, bowtie2比对至核酸序列；功能为humann2调用diamond比对至蛋白库11Gb

humann2_config # 查看参数和数据库位置是否正确
# metaphlan2数据库默认位于程序所在目录的db_v20和databases下各一份

# 单个文件运行示例
time humann2 --input temp/concat/p136C.fq  \
  --output temp/ \
  --threads 8 &
# 大约20 min

# 并行处理所有样品
time parallel -j 3 \
  'humann2 --input {}  \
  --output temp/ ' \
  ::: temp/concat/*.fq > log &
  
# 62m，累计时间1167m 约等于 20小时



# 4. 物种组成分析

mkdir -p metaphlan2

# 样品结果合并
merge_metaphlan_tables.py temp/*_humann2_temp/*_metaphlan_bugs_list.tsv | \
  sed 's/_metaphlan_bugs_list//g' > metaphlan2/taxonomy.tsv

# 如果上面运行没有结果，可使用提前过多成的文件继续分析
# cp /db/bak/taxonomy.tsv metaphlan2/

# 转换为spf格式方便stamp分析
metaphlan_to_stamp.pl metaphlan2/taxonomy.tsv > metaphlan2/taxonomy.spf

# 下载design.txt和metaphlan2.spf使用stamp分析



# 5. 物种组成分析和可视化进阶

# 绘制热图
metaphlan_hclust_heatmap.py --in metaphlan2/taxonomy.tsv \
  --out metaphlan2/heatmap.pdf \
  -c bbcry --top 25 --minv 0.1 -s log 
# c设置颜色方案，top设置物种数量，minv最小相对丰度，s标准化方法，log为取10为底对数，文件名结尾可先pdf/png/svg三种图片格式。更多说明详见 metaphlan_hclust_heatmap.py -h

# GraPhlAn图
# metaphlan2 to graphlan
export2graphlan.py --skip_rows 1,2 -i metaphlan2/taxonomy.tsv\
  --tree temp/merged_abundance.tree.txt \
  --annotation temp/merged_abundance.annot.txt \
  --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 \
  --annotations 5,6 --external_annotations 7 --min_clade_size 1
# graphlan annotation
graphlan_annotate.py --annot temp/merged_abundance.annot.txt temp/merged_abundance.tree.txt \
  temp/merged_abundance.xml
# output PDF figure, annoat and legend
graphlan.py temp/merged_abundance.xml metaphlan2/graphlan.pdf --external_legends --dpi 300 

# LEfSe差异分析和Cladogram
# 修改样本品为组名
sed '1 s/p[0-9]*//g' metaphlan2/taxonomy.tsv | grep -v '#' > metaphlan2/lefse.txt
# 格式转换为lefse内部格式
lefse-format_input.py metaphlan2/lefse.txt temp/input.in -c 1 -o 1000000
# 运行lefse
run_lefse.py temp/input.in temp/input.res
# 绘制物种树注释差异
lefse-plot_cladogram.py temp/input.res metaphlan2/lefse_cladogram.pdf --format pdf --dpi 600 
# 绘制所有差异features柱状图
lefse-plot_res.py temp/input.res metaphlan2/lefse_res.pdf --format pdf --dpi 600
# 绘制单个features柱状图(同STAMP中barplot)
# sort -k3,3n temp/input.res |less -S # 查看差异features列表
lefse-plot_features.py -f one --feature_name "k__Bacteria.p__Firmicutes" --format pdf \
  temp/input.in temp/input.res metaphlan2/Firmicutes.pdf 
# 批量绘制所有差异features柱状图
lefse-plot_features.py -f diff --archive none --format pdf \
  temp/input.in temp/input.res metaphlan2/features



# 6. 功能组成分析

mkdir -p humann2

# 合并所有样品，通路包括各功能和具体的物种组成，还有基因家族(太多)，通路覆盖度层面可以进分析
humann2_join_tables --input temp/ --file_name pathabundance --output humann2/pathabundance.tsv
sed -i 's/_Abundance//g' humann2/pathabundance.tsv

# 如果上面运行没有结果，可使用提前过多成的文件继续分析
# cp /db/bak/pathabundance.tsv metaphlan2/

# 标准化为相对丰度relab或百万分数cpm
humann2_renorm_table --input humann2/pathabundance.tsv --units relab \
  --output humann2/pathabundance_relab.tsv

# 分层结果
humann2_split_stratified_table --input humann2/pathabundance_relab.tsv \
  --output humann2/
# 结果stratified(每个菌的功能组成)和unstratified(功能组成)两个

# 筛选某个通路结果用STAMP展示
head -n 1 humann2/pathabundance_relab_stratified.tsv | \
  sed 's/# Pathway/MetaCyc_pathway/' \
  > humann2/pathabundance_relab_stratified_LACTOSECAT-PWY.spf
grep "LACTOSECAT-PWY" humann2/pathabundance_relab_unstratified.tsv \
  >> humann2/pathabundance_relab_stratified_LACTOSECAT-PWY.spf

# 可以使用stamp用分层，或末分层的结果进行统计分析


## 1.7 kraken2 reads on NCBI

# 定义输出目录、数据库位置和输入文件名
mkdir 03temp/17kraken2_reads

DBNAME=/mnt/zhou/yongxin/db/kraken2

i=HnZH11R1
# 双端序列输出metaphlan2格式报告，双端序列注释比例提高，48%-57.5%
kraken2 --db $DBNAME --paired 03temp/11qc/${i}_1_kneaddata_paired*.fastq \
	--threads 8 --use-names --use-mpa-style --report-zero-counts \
	--report 03temp/17kraken2_reads/${i}_report \
	--output 03temp/17kraken2_reads/${i}

# shell 合并文件
for file in *counts; do
  # 提取样品名
  name=${file%%.*}
  # 将每个文件中的count列改为样品名
  sed -i "1 s/count/$name/g" $file
  cut -f 2 $file > $file.count
done
# 生成基因丰度表
i=`tail -n 1 ../doc/design.txt | cut -f 1`
paste <(cut -f 1 $i.quant.counts) *count > gene_count.txt
head gene_count.txt



# 三、宏基因组无参分析流程 metagenome de novo pipeline (megahit)

# 系统要求 System: Linux Ubuntu 18.04 / CentOS 7.5
# 依赖软件 Sofware: KneadData、
# 运行前准备
# 1. 按附录说明安装软件和数据库
# 2. Rstudio打开pipeline_meta_denovo.sh文件，Terminal中切换至工作目录



# 实战流程

# 数据来自2016年mBio的文章，https://www.ncbi.nlm.nih.gov/bioproject/PRJNA278302/  

# 创建目录和准备文件
cd
mkdir -p 3meta_denovo
cd ~/3meta_denovo
mkdir -p seq
ln -s /db/3meta_denovo/seq/* ./seq/
cp -r /db/3meta_denovo/doc/ ./
cp -r /db/3meta_denovo/*.sh ./
tree
# seq下为原始数据抽样1M
# doc下为实验设计，仅2个样品用于演示
mkdir -p temp
mkdir -p result



# 1. 质量控制

# 1.1 FastQC质量评估
# 进入子目录，使输入文件没有目录更简洁
cd seq
# fastqc每个文件一个线程，2个双端样本4个文件设置4线程
fastqc *.gz -t 4

# 1.2 Trimmomatic去接头和低质量
for i in `tail -n+2 ../doc/design.txt|cut -f 1`; do
trimmomatic PE -threads 8 \
  ${i}_1.fastq.gz ${i}_2.fastq.gz \
  ${i}_1.qc.fq.gz ${i}_s1 ${i}_2.qc.fq.gz ${i}_s2 \
  ILLUMINACLIP:/db/TruSeq2-PE.fa:2:40:15 LEADING:3 \
  TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75
done

# 1.3 FastQC再评估 (可选)
fastqc *.qc.fq.gz -t 4

# 1.4 生成多样品报告比较
multiqc .
# 查看右侧seq目录中multiqc_report.html，可交互式报告

# 1.5 khmer质控
# 只对高覆盖度中的低丰度kmer剪切(更可能是测序错误)；低覆盖度保留
for i in `tail -n+2 ../doc/design.txt|cut -f 1`; do
# i=SRR1976948
# 合并双端序列，去除低频K-mer
interleave-reads.py temp/11qc/${i}_1.qc.fq.gz ${i}_2.qc.fq.gz | \
  trim-low-abund.py -V -M 8G -C 3 -Z 10 - -o ${i}.trim.fq
# 质控后拆分类单端和双端
split-paired-reads.py -f -0 ${i}_0 -1 ${i}_1.kh.fq -2 ${i}_2.kh.fq ${i}.trim.fq &
done

# 1.6 khmer过滤后评估(可选)
fastqc *.kh.fq -t 4
# 比较质控前后的fastqc
multiqc . 
# multiqc_report_1.html报告中接头进一步降低

# 1.7 比较质控前后kmer数量(可选)
i=SRR1976948
unique-kmers.py ${i}_1.qc.fq.gz ${i}_2.qc.fq.gz
# 32-mers in SRR1976948_1.qc.fq.gz: 50379132
# 32-mers in SRR1976948_2.qc.fq.gz: 42821212
# Total estimated number of unique 32-mers: 63192978
unique-kmers.py ${i}_1.kh.fq ${i}_2.kh.fq
# 32-mers in SRR1976948_1.kh.fq: 49429713
# 32-mers in SRR1976948_2.kh.fq: 42206294
# Total estimated number of unique 32-mers: 61989488
# 大数据时为下游分析速度和质量帮助较大


# 1.8 kraken物种注释 https://ccb.jhu.edu/software/kraken/
kraken # 显示帮助
i=SRR1976948
kraken -db /db/kraken -threads 24 --fastq-input --paired \
  --classified-out ${i}_csseq.fq \
  --output ${i}.krk ${i}_1.kh.fq ${i}_2.kh.fq
# 结果文件 C分类的，ID，0末分类，序列长，LCA比对结果ID
head ${i}.krk
# 转换TaxonomyID为物种描述
kraken-translate --db /db/kraken ${i}.krk > ${i}.tax
# 序列物种注释
head ${i}.tax

## 数据库下载, kraken小数据库
cd /mnt/zhou/yongxin/db/kraken
nohup wget -c http://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz &bg # 2.7%
nohup wget -c http://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz &bg # 5%

# kraken2

https://github.com/DerrickWood/kraken2

帮助 

https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.html

# 2 拼接 Assembly
cd ~/3meta_denovo

	# 单样品组装示例代码 
	mkdir -p temp/22megahit
	i=HnNrtR1
	megahit -t 12 \
		-1 03temp/11qc/${i}_1_kneaddata_paired_1.fastq \
		-2 03temp/11qc/${i}_1_kneaddata_paired_2.fastq \
		-o 03temp/22megahit/${i}
	

# 2.1 Megahit拼接
time megahit -t 24 \
  -1 seq/SRR1976948_1.kh.fq,seq/SRR1977249_1.kh.fq \
  -2 seq/SRR1976948_2.kh.fq,seq/SRR1977249_2.kh.fq \
  -o megahit &
# 默认192线程，3分54秒，实际时间445分钟，其中系统占用115m
# 设定24线程，5m, 54m, 等待1m；虽然刚才快了1min，但浪费了5倍资源
# 查看拼接结果
head megahit/final.contigs.fa

# 2.2 metaSPAdes拼接(可选)
metaspades.py -h # 详细参数说明
time metaspades.py -t 24 -m 500 \
  -1 seq/SRR1976948_1.kh.fq -1 seq/SRR1977249_1.kh.fq \
  -2 seq/SRR1976948_2.kh.fq -2 seq/SRR1977249_2.kh.fq \
  -o metaspades
# real 14m, user 214m，sys 7m

# 2.3 quast评估
quast.py -h # 显示帮助
quast.py megahit/final.contigs.fa -o megahit/
quast.py metaspades/contigs.fasta -o metaspades/
# 生成report文本tsv/txt、网页html、PDF等格式报告

# 2.4 基因注释 Prokka
ll megahit/final.contigs.fa # 47Mb
time prokka 03temp/22megahit/HnNrtR1/final.contigs.fa --outdir 03temp/24prokka/HnNrtR1 \
  --prefix mg --metagenome --kingdom Archaea,Bacteria,Bacteria,Bacteria,Mitochondria,Viruses --force --cpus 24
# 以mg开头，注释宏基因组，细菌类型，强制覆盖输出
# 默认8线程，21m, 65m；24线程，11m, 53m

# 2.4.1 建立非冗余基因集
# 合并所有基因
mkdir -p 03temp/24NRgene
# DNA level
cat 03temp/23prokka/*/mg.ffn > 03temp/24NRgene/mg.ffn
grep -c '>' 03temp/24NRgene/mg.ffn
# -c  sequence identity threshold, default 0.9; -M  max available memory (Mbyte), default 400; -l  length of throw_away_sequences, default 10
cd-hit-est -i 03temp/24NRgene/mg.ffn -o 03temp/24NRgene/mg.ffn.nr -c 0.9 -M 90000 -l 70 # -T 8
# 统计基因数据，确定ID是否非冗余
grep -c '>' 03temp/24NRgene/mg.ffn.nr
grep '>' 03temp/24NRgene/mg.ffn.nr|cut -f 1 -d ' '|sort|uniq|wc -l
# protein level
cat 03temp/23prokka/*/mg.faa > 03temp/24NRgene/mg.faa
grep -c '>' 03temp/24NRgene/mg.faa
cd-hit -i 03temp/24NRgene/mg.faa -o 03temp/24NRgene/mg.faa.nr -c 0.9 -M 90000 -l 70 # -T 8
grep -c '>' 03temp/24NRgene/mg.faa.nr



# 2.5 基因定量
# 2.5.1 转录本建索引
salmon -h # 查看帮助
salmon index -h # 索引帮助
ll prokka/mg.ffn # 注释基因大小33M
mkdir -p 04temp/25salmon
salmon index -t 03temp/24NRgene/mg.ffn.nr -p 8 \
  -i 04temp/25salmon/mRNA_index #  --type quasi -k 31
# -t 转录本序列，--type 类型fmd/quasi，-k kmer长度默认31, -i 索引
# 2.5.2 salmon定量
salmon quant -h
salmon quant --help-reads

for i in `tail -n+2 01doc/design.txt | cut -f 1` ; do
# i=HnZH11R1
# 关于-l详细 https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype
time salmon quant -i 04temp/25salmon/mRNA_index -l A -p 8 --meta \
    -1 03temp/11qc/${i}_1_kneaddata_paired_1.fastq -2 03temp/11qc/${i}_1_kneaddata_paired_2.fastq \
    -o 04temp/25salmon/${i}.quant; done

# 定量合并
salmon quantmerge -h
# 合并TPM列
salmon quantmerge --quants 04temp/25salmon/*.quant \
	-o 04temp/25salmon/gene.TPM
# 合Numreads列
salmon quantmerge --quants 04temp/25salmon/*.quant \
	--column NumReads \
	-o 04temp/25salmon/gene.NumReads
# 删除冗余标记
sed -i '1 s/.quant//g' 04temp/25salmon/gene.*



## 7 S， total 38s，一定要控制线程，默认用所有资源反而会慢10倍。把任务拆分了200份分发和回收比计算时间都长
## 结果位于quant.sf文件中
#head salmon/SRR1977249.quant/quant.sf
## 2.5.3 合并样品表
#cd salmon
#gather-counts.py # 搜索所有quant.sf文件并提取count值，于当前目录
## 查看定量结果
#head SRR1976948.quant.counts
## 更正列样品名
#for file in *counts; do
#  # 提取样品名
#  name=${file%%.*}
#  # 将每个文件中的count列改为样品名
#  sed -i "1 s/count/$name/g" $file
#  cut -f 2 $file > $file.count
#done
## 生成基因丰度表
#i=`tail -n 1 ../doc/design.txt | cut -f 1`
#paste <(cut -f 1 $i.quant.counts) *count > gene_count.txt
#head gene_count.txt
#cd ..


# 2.6 物种组成

DBNAME=/mnt/zhou/yongxin/db/kraken2
i=HnZH11R1
kraken2 --db $DBNAME --threads 8 --paired 03temp/11qc/${i}_1_kneaddata_paired_#.fastq --output --report

for i in `tail -n+2 01doc/design.txt | cut -f 1` ; do
# 
# 关于-l详细 https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype
time salmon quant -i 04temp/25salmon/mRNA_index -l A -p 8 --meta \
    -1 03temp/11qc/${i}_1_kneaddata_paired_1.fastq -2 03temp/11qc/${i}_1_kneaddata_paired_2.fastq \
    -o 04temp/25salmon/${i}.quant; done


kraken2 --db $DBNAME --threads 8 seqs.fa

# 3. Bining分箱


## 3.1 bwa比对进行contig定量

# 3.1.1 建索引
bwa index megahit/final.contigs.fa

# 3.1.2 比对
for i in `tail -n+2 doc/design.txt | cut -f 1` ; do
# i=SRR1976948
time bwa mem -t 8 megahit/final.contigs.fa \
  seq/${i}_1.kh.fq seq/${i}_2.kh.fq > temp/${i}.sam; done
    
# 3.1.3 samtools基因组索引、压缩、排序、建索引
samtools faidx megahit/final.contigs.fa
for i in temp/*.sam; do
 samtools import megahit/final.contigs.fa $i $i.bam
 samtools sort $i.bam -o $i.sorted.bam
 samtools index $i.sorted.bam
done &

# 3.1.4 查看某条contig上reads分布(可选)
# 按contig的reads数量排序，找高丰度的查看，获取reads数第9的contig
i=`grep -v ^@ temp/SRR1976948.sam | cut -f 3 | sort | uniq -c | sort -n |awk '{print $2"\t"$1}'| tail |head -n1|cut -f 1`
# 查看 i 序列400bp开始，每次拼接的结果和编号会不同
samtools tview temp/SRR1976948.sam.sorted.bam megahit/final.contigs.fa -p $i:400
# 方向可以上下左右移动查看，q退出，最好同时查看两个文件比较
samtools tview temp/SRR1977249.sam.sorted.bam megahit/final.contigs.fa -p $i:400

# 3.1.5 bedtools中的genomeCoverageBed定量contig
for i in temp/*sorted.bam; do
    genomeCoverageBed -ibam $i > ${i/.pe*/}.histogram.tab &
done
# 查看结果格式(可选)
head temp/SRR1976948.sam.sorted.bam.histogram.tab
# 1. Contig name
# 2. Depth of coverage 覆盖深度
# 3. Number of bases on contig depth equal to column 2
# 4. Size of contig (or entire genome) in base pairs
# 5. Fraction of bases on contig (or entire genome) with depth equal to column 2
# To get an esimate of mean coverage for a contig we sum (Depth of coverage) * 
# (Number of bases on contig) / (Length of the contig). We have a quick script that will do this calculation.

# 3.1.6 计算平均覆盖度，添加.coverage.tab
for hist in temp/*histogram.tab; do
    calculate-contig-coverage.py $hist
done
# 简化文件名
ll temp/*his*
rename 's/sam.sorted.bam.histogram.tab.coverage.tab/cov/' temp/*.tab
# 查看平均深度
head temp/SRR1976948.cov


# 3.2 maxbin分箱

# 建立输出目录，显示脚本帮助
mkdir -p maxbin
run_MaxBin.pl -h

# 注意多个丰度文件可制作列表，改为-abund_list
time run_MaxBin.pl \
  -contig megahit/final.contigs.fa \
  -abund temp/SRR1976948.cov -abund2 temp/SRR1977249.cov \
  -out maxbin/mb -max_iteration 50 -thread 24

# 查看结果摘要，自己评估了完整度
cat maxbin/mb.summary



# 3.3 MetaBAT分箱(可选)

mkdir -p metabat

# 统计contig覆盖度
jgi_summarize_bam_contig_depths --outputDepth metabat/depth_var.txt temp/*.sam.sorted.bam

# 运行MetaBAT
metabat -h # 显示帮助
time metabat -i megahit/final.contigs.fa -a metabat/depth_var.txt \
  --verysensitive -o metabat/mb -v --seed 315 -t 24
# 每次结果不同，设置seed才可保证重复,v输出计算过程，t线程数
# 6s, 1m36m



# 3.4 可视化VizBin

# 合并Maxbin序列和编号用于VizBin展示
# 将所有的bin文件合并，并将序列名后面添加bin编号
rm -rf maxbin/maxbin.fa
for file in maxbin/mb.*.fasta; do
    num=${file//[!0-9]/}
    sed -e "/^>/ s/$/ ${num}/" maxbin/mb.$num.fasta >> maxbin/maxbin.fa
done
# 生成一个用于bin中序列的列表
echo label > maxbin/maxbin.anno
grep ">" maxbin/maxbin.fa |cut -f2 -d ' '>> maxbin/maxbin.anno

# 合并metabat序列和编号用于VizBin展示
# 将所有的bin文件合并，并将文件名作为序列名
rm -rf metabat/metabat.fa
for file in metabat/mb.*.fa; do
    num=${file//[!0-9]/}
    sed -e "/^>/ s/$/ ${num}/" metabat/mb.$num.fa >> metabat/metabat.fa
done
# 可视化的列表
echo label > metabat/metabat.anno
grep ">" metabat/metabat.fa | cut -f2 -d ' '>> metabat/metabat.anno

# 打开maxbin/metabat中合并的文件 metabin/maxbin.fa和.anno即可展示结果



# 3.5 Bin评估CheckM https://github.com/Ecogenomics/CheckM

# 3.5.1 在物种树中鉴定位置
time checkm tree -t 8 -x fasta maxbin/ checkm
# t: 8个线程; x: fa调置文件扩展名；输入目录，输出目录,5m45s

# 3.5.2 展示统计每个bin中的标记基因信息和物种
checkm tree_qa checkm

# 3.5.3 获得标记基因集
checkm lineage_set checkm/ checkm/marker_file

# 3.5.4 分析
time checkm analyze -x fasta checkm/marker_file maxbin/ checkm
# 12m8s

# 3.5.5 统计
checkm qa checkm/marker_file checkm/ > checkm/qa.txt



# 附录1：软件安装，测试环境为Ubuntu 18.04 LTS

# 更新软件库列表
sudo apt-get update
# blast比对
sudo apt-get -y install python ncbi-blast+
# multiqc多样品质量比较
conda install multiqc
# sourmash K-mer质控
conda install sourmash
# megahit 快速组装
conda install megahit
# spades 高质量拼接
conda install spades
# idba De Bruijn Graph De Novo Assembler 拼接工具
conda install idba
# QUEST 组装评估
conda install quast
# prokka 细菌基因组注释
conda install prokka
# 升级prokka依赖的krona数据库
ktUpdateTaxonomy.sh /conda2/opt/krona/taxonomy
# sourmash K-mer分析比较工具
conda install sourmash
# 定量工具salmon
conda install salmon
# 合并定量结果脚本
wget https://raw.githubusercontent.com/ngs-docs/2016-metagenomics-sio/master/gather-counts.py
chmod +x gather-counts.py
sudo mv gather-counts.py /usr/local/bin
# 比对工具bwa
conda install bwa
# 计算contig平均覆盖度脚本，用于bin
wget https://raw.githubusercontent.com/ngs-docs/2017-cicese-metagenomics/master/files/calculate-contig-coverage.py
sudo mv calculate-contig-coverage.py /usr/local/bin
# 分箱工具maxbin2
conda install maxbin2
wget https://downloads.jbei.org/data/microbial_communities/MaxBin/getfile.php?MaxBin-2.2.5.tar.gz
mv getfile.php\?MaxBin-2.2.5.tar.gz MaxBin-2.2.5.tar.gz
tar xvzf MaxBin-2.2.5.tar.gz
cd ./MaxBin-2.2.5/src/
make
cd ../..
cp -r MaxBin-2.2.5 /conda2/ 
# 分箱工具metabat
conda config --add channels ursky
conda install -c ursky metabat2
# 分箱评估
conda install checkm-genome
# 分箱评估数据库下载
mkdir /db/checkm
checkm data setRoot # 输入 /db/checkm
cd /db/checkm
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz
# Bin可视化VizBin
sudo apt-get install libatlas3-base libopenblas-base default-jre
curl -L https://github.com/claczny/VizBin/blob/master/VizBin-dist.jar?raw=true > VizBin-dist.jar
mv VizBin-dist.jar /usr/local/bin # 或~/bin
# 物种分类karken
apt install kraken
db=/db/kraken
mkdir -p $db
kraken-build --standard --threads 24 --db $db # > 450GB space and > 250GB RAM, 24cores > 5hours
kraken-build --db $db --clean

# Circos——绘制基因组圈图
conda install circos

# Anvi'o工具箱——组装assembly和分箱bin结果可视化
conda install anvio
# 分析采用网页Pythone服务器8080端口，

# 物种注释和分箱流程 https://github.com/bxlab/metaWRAP
conda create -n metawrap python=2.7
source activate metawrap
conda install -c ursky metawrap-mg 



## 附录

### 附录1. 宏基因组样本多文件合并和统计

	# 水稻测试数据6个样品，实验和对照各3个，数据量109-236M，PE100, 21.8-47.2G，共198.8GB
	# 拟南芥测试数据36个样本，4个实验组，共850GB数据
	# 样有多个文件，合并各样品文件()
	# merge_sample ~ 5h, 合并单个样品, 20GB 10min;
	cd /mnt/m2/data/meta/ath/3T/
	# 按实验设计按文件夹批量合并再改名，需要输入文件每个样本一个目录
	p=3
	for i in `tail -n+2 design.txt|cut -f3`; do
		zcat `find 01.filter/${i}/ -name *.gz | grep '_1.fq'` | pigz -p ${p} > seq/${i}_1.fq.gz &
		zcat `find 01.filter/${i}/ -name *.gz | grep '_2.fq'` | pigz -p ${p} > seq/${i}_2.fq.gz &
	done
	awk 'BEGIN{OFS=FS="\t"}{system("mv seq/"$3"_1.fq.gz seq/"$1"_1.fq.gz ");system("mv seq/"$3"_2.fq.gz seq/"$1"_2.fq.gz ");}' <(tail -n+2 design.txt)
	# 统计样本md5值
	cd seq
	md5sum *_1.fq.gz > md5sum.txt
	md5sum *_2.fq.gz >> md5sum.txt
	cat md5sum.txt

### 附录2. KneadData中断的手动继续
	# qc中kneaddata多任务中断后仅重跑不完整样本的手动处理
	# 检查log中报告完成的样本，并提取样本名
	grep 'Final output files created' temp/11qc/*.log|cut -f 3 -d '/'|cut -f 1 -d '_' > temp/qc.finished
	# 筛选末完成样本
	cat temp/qc.finished <(tail -n+2 result/design.txt) | cut -f 1 | sort | uniq -u > temp/qc.unfinished
	# 手动运行末完成样本
time parallel --xapply -j 8 \
  "kneaddata -i seq/{1}_1.fq.gz -i seq/{1}_2.fq.gz \
  -o temp/11qc -v -t 12 --remove-intermediate-output \
  --trimmomatic /conda/share/trimmomatic-0.38-1/ --trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:50' \
  --bowtie2-options '--very-sensitive --dovetail' -db /db/host/human/Homo_sapiens" \
  ::: `cat temp/qc.unfinished`
  kneaddata_read_count_table --input temp/11qc --output temp/11kneaddata_stat.txt
  cut -f 1,2,4,12 temp/11kneaddata_stat.txt | awk 'BEGIN{OFS=FS="\t"} {print $0,$3/$2*100,$4/$3*100}' | sed 's/_1_kneaddata//' | sed '1 s/-nan/Hi-Q%/;s/-nan/rm_host%/' > result/11kneaddata_stat.txt
  cat result/11kneaddata_stat.txt
	# 常见错误
	# Critical Error: Unable to gunzip input file: seq/QBD7BPX1T_1.fq.gz
	# zlib.error: Error -3 while decompressing: invalid distance too far back
	# fastq文件可能传输中出错，重复找测序公司重传；这些数据一般用fastqc评估也无法通过

### 检查双端配对
	
	# med中25个样本，lyk9开始第二批测序中的15个全有错误
	for i in `tail -n+2 result/metadata.txt|cut -f1`; do
	# i=lyk9b3r10
	echo $i
	seqkit seq -n -i temp/qc/${i}_1_kneaddata_paired_1.fastq|cut -f 1 -d '/' > temp/header_${i}_1
	seqkit seq -n -i temp/qc/${i}_1_kneaddata_paired_2.fastq|cut -f 1 -d '/' > temp/header_${i}_2
	#md5sum temp/header_${i}_1 > temp/md5sum_${i}_1
	#md5sum temp/header_${i}_2 > temp/md5sum_${i}_2
	#cat temp/md5sum_${i}* | cut -f 1 -d ' ' | uniq -u
	cmp temp/header_${i}_1 temp/header_${i}_2 | cut -f 2 -d '_' >> temp/diff.list
	done

	# 错误的修复，无配对序列
	i=lyk9b3r10
	mkdir -p temp/pair/
	time seqkit pair -1 temp/qc/${i}_1_kneaddata_paired_1.fastq -2 temp/qc/${i}_1_kneaddata_paired_2.fastq -O temp/pair/ -u

	# 检查错误来源，以小样本为例：来源去宿主后即不配对
	# kraken2_qc过滤后，因为只是head提前，肯定文件本来有问题
	# 去宿主后，存在ID不一致
	i=lyk9b3r2
	seqkit seq -n -i submit/${i}_1.fq.gz|cut -f 1 -d '/' > submit/header_${i}_1
	seqkit seq -n -i submit/${i}_2.fq.gz|cut -f 1 -d '/' > submit/header_${i}_2
	ls -l submit/header_${i}_?
	cmp submit/header_${i}_? | cut -f 2 -d '_' >> submit/diff.list
	cat submit/diff.list
	
	# 去宿主前(ID一致)
	i=`tail -n+12 result/metadata.txt | head -n 1 | cut -f 1`
	# 18G,9m,29M；多线程序列一样，时间一样。效果没有明显提交
	seqkit stat seq/${i}_1.fq.gz
	memusg -t seqkit seq -n -i -j 10 seq/${i}_1.fq.gz |cut -f 1 -d '/' > seq/header_${i}_1
	memusg -t seqkit seq -n -i seq/${i}_2.fq.gz|cut -f 1 -d '/' > seq/header_${i}_2
	ls -l seq/header_${i}_?
	cmp seq/header_${i}_? | cut -f 2 -d '_' >> seq/diff.list
	cat seq/diff.list

	# 去宿主后
	# ID不致
	memusg -t seqkit seq -n -i temp/11qc/${i}_1_kneaddata_paired_1.fastq|cut -f 1 -d '/' | head > temp/header_${i}_1
	memusg -t seqkit seq -n -i temp/11qc/${i}_1_kneaddata_paired_2.fastq|cut -f 1 -d '/' | head > temp/header_${i}_2
	# 但行数一致
	wc -l temp/11qc/${i}_1_kneaddata_paired_1.fastq temp/11qc/${i}_1_kneaddata_paired_2.fastq
	
	# 质控去宿主部分有错，重新运行10个样
	conda activate meta
	# 以lyk9b3r1为例测试结果是否成对
	# 36G,2h,43G; 0.3G,1m,1G
	i=`tail -n+12 result/metadata.txt | head -n 1 | cut -f 1`
	zcat seq/${i}_1.fq.gz|head -n 4000000 > seq/${i}_1.fq
	zcat seq/${i}_2.fq.gz|head -n 4000000 > seq/${i}_2.fq
	memusg -t parallel --xapply -j 3 \
	  "kneaddata -i seq/{1}_1.fq -i seq/{1}_2.fq \
	  -o temp/qc -v -t 48 --remove-intermediate-output \
	  --trimmomatic /conda/envs/meta/share/trimmomatic/ --trimmomatic-options 'ILLUMINACLIP:/conda/envs/meta/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:100' \
	  --reorder -db /db/kneaddata/medicago/bt2" \
	  ::: `tail -n+12 result/metadata.txt | head -n 1 | cut -f 1`
	# 去宿主后
	memusg -t seqkit seq -n -i temp/qc/${i}_1_kneaddata_paired_1.fastq|cut -f 1 -d '/' | head > temp/header_${i}_1
	memusg -t seqkit seq -n -i temp/qc/${i}_1_kneaddata_paired_2.fastq|cut -f 1 -d '/' | head > temp/header_${i}_2
	# 行数一致，但序列ID不一致
	wc -l temp/11qc/${i}_1_kneaddata_paired_1.fastq temp/11qc/${i}_1_kneaddata_paired_2.fastq
	cmp temp/header_${i}_?


	# 依赖软件版本：改名、压缩、链接、md5值计算、
	parallel --version # 20201122
	csvtk version # v0.22.0
	rename --version # 0.20
	pigz --version # 2.4
	ln --version # (GNU coreutils) 8.28, 包括Linux常用命令md5sum, cat
	sed --version # 4.8
	md5sum --version # 8.28
