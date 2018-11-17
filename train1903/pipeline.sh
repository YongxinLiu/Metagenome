# 一、分析准备工作

# 系统要求 System: Linux Ubuntu 18.04 / CentOS 7.5
# 依赖软件 Sofware: KneadData、MetaPhlAn2、HUMAnN2
# 运行前准备
# 1. 第一次使用，请参见`soft_db.sh`脚本，安装分析所需软件和数据库(3-5天)
# 2. 学员U盘和服务器家目录(~)有项目文件夹meta, 包含测序数据和实验设计
# 3. 使用谷歌浏览器访问服务器，具体IP地址课上通知
# 4. 网页版Rstudio打开~/meta/pipeline.sh流程文件


# 新用户无文件时上传或复制准备的文件(己有文件请跳过)
mkdir -p ~/meta && cd ~/meta
cp -r /db/meta/seq ./
cp /db/meta/*.sh ./
mkdir -p temp result
cp /db/meta/result/design.txt result/

## 1.1 了解工作目录和文件

# Terminal中切换至工作目录
cd ~/meta
# 显示文件结构
tree 
# .
# ├── pipeline.sh
# ├── result
# │   └── design.txt
# ├── seq
# │   ├── p136C_1.fq.gz
# │   ├── p136C_2.fq.gz
# │   ├── p136N_1.fq.gz
# │   ├── p136N_2.fq.gz
# │   ├── p144C_1.fq.gz
# │   ├── p144C_2.fq.gz
# │   ├── p144N_1.fq.gz
# │   ├── p144N_2.fq.gz
# │   ├── p153C_1.fq.gz
# │   ├── p153C_2.fq.gz
# │   ├── p153N_1.fq.gz
# │   └── p153N_2.fq.gz
# ├── soft_db.sh
# └── temp
# pipeline.sh是分析流程代码；
# seq目录中有6个样本双端测序，共12个序列文件；
# temp是临时文件夹，存储分析中间文件，结束可全部删除节约空间
# result是重要节点文件和整理化的分析结果图表，实验设计design.txt也在此


## 1.2 FastQC质量评估(可选)

# 进入子目录，使输入文件没有目录更简洁
# fastqc每个文件一个线程，6个双端样本12个文件，设置6线程
fastqc seq/*.gz -t 6 # 9s
# 结果见seq目录，解读见[数据的质量控制软件——fastQC](https://mp.weixin.qq.com/s/MMfierO-8H2MVRkJKGCOVQ)

# 生成多样品报告比较
multiqc -d seq/ -o result/qc
# 查看右侧result/qc目录中multiqc_report.html，可交互式报告
# GC都不同，组间差异太明显了


## 1.3 序列质控和去宿主 Qaulity control and remove host contamination

# kneaddata是一套工作流，依赖trimmomatic进行质控和去接头；
# 依赖bowtie2比对宿主，并筛选非宿主序列
# kneaddata -h # 显示帮助

# 以p144C单样品质控为例，readl 9s, user 18s.
# time kneaddata -i seq/p144C_1.fq.gz -i seq/p144C_2.fq.gz \
#   -o temp/qc -v -t 3 --remove-intermediate-output \
#   --trimmomatic /conda2/share/trimmomatic-0.36-3/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
#   --bowtie2-options "--very-sensitive --dovetail" -db /db/bowtie2/Homo_sapiens   

# 现实中是有一大堆样品，你可以逐个修改样品名运行，如果服务器性能允许可以并行加速分析
# parallel --citation # 打will cite以后不再提醒
time parallel -j 3 --xapply \
  'kneaddata -i {1} -i {2} \
  -o temp/qc -v -t 3 --remove-intermediate-output \
  --trimmomatic /conda2/share/trimmomatic-0.36-3/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
  --bowtie2-options "--very-sensitive --dovetail" -db /db/bowtie2/Homo_sapiens' \
 ::: seq/*_1.fq.gz ::: seq/*_2.fq.gz
# real 16s, user 1m28s

# 质控结果统计
kneaddata_read_count_table --input temp/qc --output result/01kneaddata_sum.txt
cat result/01kneaddata_sum.txt


# 1.4 质控后质量再评估 (可选)
fastqc temp/qc/*_1_kneaddata_paired_* -t 6
multiqc -d temp/qc/ -o result/qc/
# 整理bowtie2, trimmomatic, fastqc报告



# 二、宏基因组有参分析流程 Reference-based Metagenomic Pipeline (HUMAnN2)

# 中文教程：https://mp.weixin.qq.com/s/XkfT5MAo96KgyyVaN_Fl7g
# 英文教程：https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Humann2)


## 2.1 合并质控后文件

# 有参宏基因组不考虑双端，HUMAnN2是求将双端序列共合并后作为输入
mkdir -p temp/concat
for i in `tail -n+2 result/design.txt | cut -f 1`;do \
  cat temp/qc/${i}*_1_kneaddata_paired_?.fastq > temp/concat/${i}.fq; done
ls -l temp/concat/*.fq # 查看样品数量和大小


## 2.2 HUMAnN2计算物种和功能组成

# 物种组成调用MetaPhlAn2, bowtie2比对至核酸序列；功能为humann2调用diamond比对至蛋白库11Gb

# 单个文件运行示例，大约20 min
# time humann2 --input temp/concat/p136C.fq  \
#   --output temp/ \
#   --threads 8 &

# 并行HUMAnN2处理所有样品
parallel --citation # 输入提示文件
mkdir -p temp/humann2
time parallel -j 3 \
  'humann2 --input {}  \
  --output temp/humann2/ ' \
  ::: temp/concat/*.fq > temp/log &
cat temp/log
# 3 jobs X 8p; real 31m, user 583m, sys 29m
# 6 jobs X 8p; real 17m, user 613m, sys 27m
# 核心步骤，测序数据3X8=24线程，用时30min，真实数据可能要几小时至几天


## 2.3 功能组成分析

mkdir -p result/humann2

# 合并所有样品，通路包括各功能和具体的物种组成，还有基因家族(太多)，通路覆盖度层面可以进分析
humann2_join_tables --input temp/humann2 --file_name pathabundance \
  --output result/humann2/pathabundance.tsv
sed -i 's/_Abundance//g' result/humann2/pathabundance.tsv

# 标准化为相对丰度relab或百万分数cpm
humann2_renorm_table --input result/humann2/pathabundance.tsv --units relab \
  --output result/humann2/pathabundance_relab.tsv

# 分层结果
humann2_split_stratified_table --input result/humann2/pathabundance_relab.tsv \
  --output result/humann2/
# 结果stratified(每个菌的功能组成)和unstratified(功能组成)两个
# 可以使用stamp用分层，或末分层的结果进行统计分析


## 2.4 物种组成分析

mkdir -p result/metaphlan2

### 2.4.1 样品结果合并
merge_metaphlan_tables.py temp/humann2/*_humann2_temp/*_metaphlan_bugs_list.tsv | \
  sed 's/_metaphlan_bugs_list//g' > result/metaphlan2/taxonomy.tsv

### 2.4.2 转换为spf格式方便stamp分析
metaphlan_to_stamp.pl result/metaphlan2/taxonomy.tsv > result/metaphlan2/taxonomy.spf
# 下载design.txt和taxonomy.spf使用stamp分析

### 2.4.3 绘制热图
metaphlan_hclust_heatmap.py --in result/metaphlan2/taxonomy.tsv \
  --out result/metaphlan2/heatmap.pdf \
  -c bbcry --top 50 --minv 0.1 -s log 
# c设置颜色方案，top设置物种数量，minv最小相对丰度，s标准化方法，log为取10为底对数，文件名结尾可先pdf/png/svg三种图片格式。更多说明详见 metaphlan_hclust_heatmap.py -h

### 2.4.4 GraPhlAn图
# metaphlan2 to graphlan
export2graphlan.py --skip_rows 1,2 -i result/metaphlan2/taxonomy.tsv \
  --tree temp/merged_abundance.tree.txt \
  --annotation temp/merged_abundance.annot.txt \
  --most_abundant 1000 --abundance_threshold 20 --least_biomarkers 10 \
  --annotations 3,4 --external_annotations 7
# graphlan annotation
graphlan_annotate.py --annot temp/merged_abundance.annot.txt \
  temp/merged_abundance.tree.txt  temp/merged_abundance.xml
# output PDF figure, annoat and legend
graphlan.py temp/merged_abundance.xml result/metaphlan2/graphlan.pdf \
  --external_legends 

# 2.4.5 LEfSe差异分析和Cladogram
# 修改样本品为组名
sed '1 s/p[0-9]*//g' result/metaphlan2/taxonomy.tsv | grep -v '#' > result/metaphlan2/lefse.txt
# 格式转换为lefse内部格式
lefse-format_input.py result/metaphlan2/lefse.txt temp/input.in -c 1 -o 1000000
# 运行lefse
run_lefse.py temp/input.in temp/input.res
# 绘制物种树注释差异
lefse-plot_cladogram.py temp/input.res result/metaphlan2/lefse_cladogram.pdf --format pdf --dpi 600 
# 绘制所有差异features柱状图
lefse-plot_res.py temp/input.res result/metaphlan2/lefse_res.pdf --format pdf --dpi 600
# 绘制单个features柱状图(同STAMP中barplot)
grep -v '-' temp/input.res | sort -k3,3n  # 查看显著差异features，按丰度排序
lefse-plot_features.py -f one --feature_name "k__Bacteria.p__Firmicutes" --format pdf \
  temp/input.in temp/input.res result/metaphlan2/Firmicutes.pdf 
# 批量绘制所有差异features柱状图
lefse-plot_features.py -f diff --archive none --format pdf \
  temp/input.in temp/input.res result/metaphlan2/features



# 三、宏基因组无参分析流程 De novo Metagenomic Pipeline (MEGAHIT)

## 3.1 kraken2基于NCBI数据库注释reads层面
# 还可以进行contig、gene、bin层面的序列物种注释

### 3.1.1 多样本并行物种注释
mkdir -p temp/kraken2
time parallel -j 3 \
  'kraken2 --db /db/kraken2 --paired temp/qc/{1}_1_kneaddata_paired*.fastq \
  --threads 3 --use-names --use-mpa-style --report-zero-counts \
  --report temp/kraken2/{1}_report \
  --output temp/kraken2/{1}_output' \
  ::: `tail -n+2 result/design.txt | cut -f 1`
# 屏幕会输出各样品注释比例，和运行时间 real 39s, user 11s, sys 1m54s

### 3.1.2 汇总样品物种组成表
mkdir -p result/kraken2
parallel -j 6 \
  'cut -f 2 temp/kraken2/{1}_report | sed "1 s/^/{1}\n/" > temp/kraken2/{1}_count ' \
  ::: `tail -n+2 result/design.txt | cut -f 1`
header=`tail -n 1 result/design.txt | cut -f 1`
cut -f 1 temp/kraken2/${header}_report | sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
paste temp/kraken2/*count > result/kraken2/taxonomy_count.txt


## 3.2 拼接 Assembly

### 3.2.1 Megahit拼接
# 删除输去旧文件夹，否则不重新计算
rm -rf temp/megahit
# 组装核心程序，计算和内存需求密集，几百G或上T数据需几天至几周
time megahit -t 9 \
  -1 `ls temp/qc/*_1_kneaddata_paired_1.fastq|tr '\n' ','|sed 's/,$//'` \
  -2 `ls temp/qc/*_1_kneaddata_paired_2.fastq|tr '\n' ','|sed 's/,$//'` \
  -o temp/megahit # --k-min 27 --k-max 191 --k-step 20
# `ls...`用于自动获得文件列表，无需手动添写
# 默认使用所有线程，如本机192线程，real 43s, user 96m, sys 1m
# 设定9线程，real 1m, user 8m, sys 1m；虽才快了17s，但浪费了10倍资源
# 查看拼接结果
head temp/megahit/final.contigs.fa

### 3.2.2 quast评估
mkdir -p result/megahit/
ln temp/megahit/final.contigs.fa result/megahit/
quast.py result/megahit/final.contigs.fa -o result/megahit/
# 生成report文本tsv/txt、网页html、PDF等格式报告


## 3.3 基因注释和定量 Gene annotation & quantitfy

### 3.3.1 prokka基因注释
# 查看文件大小，预估时间
ll temp/megahit/final.contigs.fa # 47Mb
time prokka temp/megahit/final.contigs.fa --outdir temp/prokka \
  --prefix mg --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
  --force --cpus 8
# 以mg开头，注释宏基因组，多界，强制覆盖输出
# 8线程，耗时30s, 1m32s, 29s

### 3.3.2 cd-hit构建非冗余基因集
mkdir -p temp/NR
# 输入文件可由多个样本、组、批次序列合并文件，方便整合分析
# aS覆盖度，c相似度，G局部比对，M内存0不限制，T多线程，g最优解
time cd-hit-est -i temp/prokka/mg.ffn -o temp/NR/mg.ffn.nr \
  -aS 0.9 -c 0.95 -G 0 -M 0 -T 9 -g 1
# real 0.3s,  user 1.8s
# 挑选重要结果保存到结果目录
mkdir -p result/NR
ln temp/NR/mg.ffn.nr result/NR/nucleotide.fa
# 翻译核酸为对应蛋白序列
transeq -sequence result/NR/nucleotide.fa -outseq result/NR/protein.fa -trim Y 
# 序列名自动添加了_1，为与核酸对应要去除
sed -i 's/_1 / /' result/NR/protein.fa

### 3.3.3 基因定量salmon
mkdir -p temp/salmon
# 建索引,-t转录本，--type类型fmd/quasi，-k kmer长度默认31, -i 索引
salmon index -t result/NR/nucleotide.fa -p 9 -k 31 \
  -i temp/salmon/index # --type quasi 
# 定量，l文库类型自动选择，p线程 ，--meta宏基因组模式, 3s, 8s
time parallel -j 3 \
  'salmon quant -i temp/salmon/index -l A -p 3 --meta \
  -1 temp/qc/{1}_1_kneaddata_paired_1.fastq -2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
  -o temp/salmon/{1}.quant' \
  ::: `tail -n+2 result/design.txt | cut -f 1`
# 合并
mkdir -p result/salmon
salmon quantmerge --quants temp/salmon/*.quant \
  -o result/salmon/gene.TPM
salmon quantmerge --quants temp/salmon/*.quant \
  --column NumReads -o result/salmon/gene.count
sed -i '1 s/.quant//g' result/salmon/gene.*




## 3.4 功能基因注释


### 3.4.1 基因注释eggNOG/COG/KEGG/GO

# diamond比对基因至eggNOG数据库, real 13m, user 122m
mkdir -p temp/eggnog
time emapper.py -m diamond --no_annot --no_file_comments --data_dir /db/eggnog \
  --cpu 9 -i result/NR/protein.fa -o temp/eggnog/protein --override
# 比对结果功能注释, real 6s, user 44s 
time emapper.py --annotate_hits_table temp/eggnog/protein.emapper.seed_orthologs --no_file_comments \
		-o temp/eggnog/output --cpu 9 --data_dir /db/eggnog --override

# 结果注释表头, 重点1序列名，7KO，12COG分类，13注释
mkdir -p result/eggnog
sed '1 i Name\teggNOG\tEvalue\tScore\tGeneName\tGO\tKO\tBiGG\tTax\tOG\tBestOG\tCOG\tAnnotation' \
  temp/eggnog/output.emapper.annotations > temp/eggnog/output

# 1.整理COG表
# 提取12列COG分类
cut -f 1,12 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$'>temp/eggnog/1cog.list
# 基因丰度矩阵末尾添加对应cog编号，没注释的直接删除
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/eggnog/1cog.list result/salmon/gene.count | \
	sed '/\t$/d' | sed '1 s/COG/KO/' > temp/eggnog/gene_cog.count
# 按COG类型合并表格，输出count和RPM值，n设置标准化单位，默认1M，可选100/1
mat_gene2ko.R -i temp/eggnog/gene_cog.count -o result/eggnog/cogtab -n 1000000
# STAMP的spf格式，结果design.txt进行KO或Description差异比较
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
  /db/eggnog/COG.anno result/eggnog/cogtab.count > result/eggnog/cogtab.count.spf

# 2. 整理KO表
# 提取基因KO表，基因1对多个KO时只提取第一个KO
cut -f 1,7 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' > temp/eggnog/2ko.list
# 基因丰度矩阵末尾添加对应KO编号，没注释的直接删除
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/eggnog/2ko.list \
  result/salmon/gene.count | sed '/\t$/d' > temp/eggnog/gene_ko.count
# 合并基因表为KO表，输出count值和tpm值
mat_gene2ko.R -i temp/eggnog/gene_ko.count -o result/eggnog/kotab -n 1000000
# STAMP的spf格式，结果design.txt进行KO或Description差异比较
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' /db/eggnog/KO.anno \
  result/eggnog/kotab.count | sed 's/^\t/Undescription\t/' > result/eggnog/kotab.count.spf



### 3.4.2 碳水化合物dbCAN2
# 比对CAZy数据库, 35s, 5m20s
mkdir -p temp/dbcan2
time diamond blastp --db /db/dbcan2/CAZyDB.07312018 --query result/NR/protein.fa \
	--outfmt 6 --threads 9 --max-target-seqs 1 --quiet \
	--out temp/dbcan2/gene_diamond.f6
# 整理比对数据为表格
mkdir -p result/dbcan2
# 提取基因对应基因家族，同一基因存在1对多，只取第一个
cut -f 1,2 temp/dbcan2/gene_diamond.f6 | uniq | sed 's/|/\t/g' | cut -f 1,3 | \
	cut -f 1,2 -d '_' |sed '1 i Name\tKO' > temp/dbcan2/gene_fam.list
# 基因丰度矩阵末尾添加对应FAM编号，没注释的直接删除
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/dbcan2/gene_fam.list \
  result/salmon/gene.count | sed '/\t$/d' > temp/dbcan2/gene_fam.count
# 按基因家族合并
mat_gene2ko.R -i temp/dbcan2/gene_fam.count -o result/dbcan2/cazytab
# 结果中添加FAM注释，spf格式用于stamp分析
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' /db/dbcan2/fam_description.txt \
	result/dbcan2/cazytab.count > result/dbcan2/cazytab.count.spf



### 3.4.3 抗生素抗性ResFam

mkdir -p temp/resfam
# 比对至抗生素数据库 1s, 8s
time diamond blastp --db /db/resfam/Resfams-proteins.dmnd --query result/NR/protein.fa \
	--outfmt 6 --threads 9 --max-target-seqs 1 --quiet \
	--out temp/resfam/gene_diamond.f6

mkdir -p result/resfam
# 提取基因对应基因家族
cut -f 1,2 temp/resfam/gene_diamond.f6 | uniq | \
  sed '1 i Name\tResGeneID' > temp/resfam/gene_fam.list
# 基因丰度矩阵末尾添加对应FAM编号，没注释的直接删除
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
  temp/resfam/gene_fam.list result/salmon/gene.count | \
	sed '/^\t/d' > result/resfam/resfam.count
# 统计注释基因的比例
wc -l result/23salmon_gene/gene.count
wc -l result/resfam/resfam.count # 172/7734=2.2%
# 结果中添加FAM注释，spf格式用于stamp分析
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4"\t"$3"\t"$2} NR>FNR{print a[$1],$0}' \
  /db/resfam/Resfams-proteins_class.tsv  result/resfam/resfam.count \
  > result/resfam/resfam.count.spf




# 四、挖掘单菌基因组/分箱 Binning (MetaWRAP)

# 挖掘单菌基因组，需要研究对象复杂度越低、测序浓度越大，结果越好。推荐数据量30GB+，至少3个样本起。
# 上面的测序数据6个样一共有70MB，适合演示分析流程但无法挖掘到单菌基因组，这里使用官方测序数据演示讲解
# 软件和数据库布置需2-3天，演示数据分析过程超10h，标准30G样也需3-30天，由服务器性能决定。

# 准备原始数据从头分析，详见公众号或官网教程
# 这里我们从质控后数据和拼接结果开始
mkdir -p binning && cd binning
mkdir -p temp
# 7G质控数据
cp -r /db/meta/binning/temp/qc ./temp
# 拼接结果
cp -r /db/meta/binning/temp/megahit ./temp
# 这里基于质控clean数据和拼接好的contigs


# 加载运行环境
cd ~/meta/binning
source activate metawrap

## 4.1.1 物种组成预测(kraken+krona，可选)

# kraken物种注释, -o输出目录, -t线程数, -s抽样1M减少计算量, 质控数据, contig文件
# 7G演示数据，8线程耗时3m17s
metawrap kraken -o temp/kraken -t 8 -s 1000000 \
  temp/qc/ERR*.fastq temp/megahit/final.contigs.fa
# 结果文件夹中有注释结果文件*.kraken(每个reads的注释结果)、*.krona(注释结果分类汇总)
# 可视化的Krona网页图表kronagram.html，链接结果到result方便查看
mkdir -p result
ln temp/kraken/kronagram.html result/

### 4.1.2 运行三种bin软件

# 8线程耗时34m
metawrap binning -o temp/binning -t 8 -a temp/megahit/final.contigs.fa \
  --metabat2 --maxbin2 --concoct temp/qc/ERR*.fastq

### 4.1.3 Bin提纯

# 调用三大主流binning程序cococt, maxbin2, metabat2
# 一般要求c完整度70，x污染率5，这里数据少降低阈值保证有结果作演示
# 8线程耗时67m
metawrap bin_refinement -o temp/bin_refinement -t 8 \
  -A temp/binning/metabat2_bins/ \
  -B temp/binning/maxbin2_bins/ \
  -C temp/binning/concoct_bins/ \
  -c 50 -x 10
# 结果改进程度见temp/bin_refinement/figures/目录
# ImportError: /usr/lib/x86_64-linux-gnu/libk5crypto.so.3: symbol k5_os_mutex_destroy version krb5support_0_MIT not defined in file libkrb5support.so.0 with link time reference

### 4.1.4 Bin可视化

# 计算每个contig的GC含量和在每个样本中的丰度
# 8线程耗时31m
metawrap blobology -a temp/megahit/final.contigs.fa -t 8 \
  -o temp/bloblogy --bins temp/bin_refinement/metaWRAP_bins \
  temp/qc/ERR*.fastq
# 结果为final.contigs.binned.blobplot，方便使用ggplot2可视化
# 参考脚本/conda2/envs/metawrap/bin/metawrap-scripts/blobology/makeblobplot_with_bins.R

### 4.1.5 Bin定量

# 使用salmon计算每个bin在样本中相对丰度
# 耗时1m37m，系统用时43m，此处可设置线程，但salmon仍尽可能多调用资源
metawrap quant_bins -b temp/bin_refinement/metaWRAP_bins -t 8 \
  -o temp/quant_bins -a temp/megahit/final.contigs.fa temp/qc/ERR*.fastq

# 结果包括bin丰度热图`temp/quant_bins/genome_abundance_heatmap.png`。
# 如果想自己画图，原始数据位于`temp/quant_bins/abundance_table.tab`

### 4.1.7 重组装

# 需合并所有样本作为此步输入,6s
cat temp/qc/ERR*_1.fastq > temp/qc/all_1.fq
cat temp/qc/ERR*_2.fastq > temp/qc/all_2.fq
# 提纯的bin还可以通过再组装进一步改善结果，8核，100G内存，用时43m
metawrap reassemble_bins -o temp/reassemble_bins \
  -1 temp/qc/all_1.fq -2 temp/qc/all_2.fq -t 8 -m 100 \
  -c 50 -x 10 -b temp/bin_refinement/metaWRAP_bins
# 结果统计见`temp/reassemble_bins/reassembled_bins.stats`，
# `temp/reassemble_bins/reassembly_results.png`，比对重组装前后的变化，N50、这完整度和污染率均有改进。
# `temp/reassemble_bins/reassembled_bins.png`展示CheckM对bin评估结果的可视化。


## 4.1.8 bin物种注释
# Taxator-tk 进行每条contig物种注释，再估计bin整体的物种
# 11m
metawrap classify_bins -b temp/reassemble_bins/reassembled_bins \
  -o temp/classify_bins -t 8
# 注释结果见`temp/classify_bins/bin_taxonomy.tab`


## 4.1.9 Bin功能注释
# 基于prokka基因注释，4m
metaWRAP annotate_bins -o temp/annotate_bins -t 8 \
  -b temp/reassemble_bins/reassembled_bins/
# 每个bin基因注释的gff文件见 `FUNCT_ANNOT/bin_funct_annotations`

# 分析结束，退出虚拟环境
source deactivate

# 其它单菌基因组分析可结合进化、COG、KEGG、GO、碳水化合物、
# 抗生素抗性和基因簇分析(如antismash)等