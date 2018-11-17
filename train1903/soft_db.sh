# 附录1：软件和数据库安装，测试环境为Ubuntu 18.04 LTS


# 软件管理器miniconda2
mkdir -p soft && cd soft
wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
# 正常默认即可，为方便多人使用目录更短，我安装在/conda目录
# 安装说明详见：[Nature Method：Bioconda解决生物软件安装的烦恼](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)
# 添加通道
conda config --add channels conda-forge
conda config --add channels bioconda
# 添加清华镜像加速下载
site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
conda config --add channels ${site}/pkgs/free/ 
conda config --add channels ${site}/pkgs/main/
conda config --add channels ${site}/cloud/conda-forge/
conda config --add channels ${site}/pkgs/r/
conda config --add channels ${site}/cloud/bioconda/
conda config --add channels ${site}/cloud/msys2/
conda config --add channels ${site}/cloud/menpo/
conda config --add channels ${site}/cloud/pytorch/


# 数据库安装位置
# 多人使用，我们统一安装至/db目录
db=~/db
mkdir -p $db

## 1.1 质控软件
# 质量评估软件fastqc
conda install fastqc # 0.11.8
# 多样品评估报告汇总multiqc
conda install multiqc # v1.6

# 质量控制流程kneaddata
conda install kneaddata # v0.7.0
# 查看可用数据库，如宏基因组去宿主，宏转录组去核糖体
kneaddata_database
mkdir -p $db/kneaddata
# 如下载人类基因组bowtie2索引，3.44G，推荐晚上过夜下载
kneaddata_database --download human_genome bowtie2 $db/kneaddata/
# 如宏转录组用SILVA核糖体数据库，3.61G
kneaddata_database --download ribosomal_RNA bowtie2 $db/kneaddata/


## 1.2 有参分析流程MetaPhlAn2、HUMAnN2

# 安装MetaPhlAn2、HUMAnN2和所有依赖关系
conda install humann2 # 0.11.1
# 测试流程是否可用
humann2_test

humann2_databases # 显示可用数据库
wd=~/db/humann2
mkdir -p $wd # 建立下载目录
# 微生物物种核心基因 5.37G
humann2_databases --download chocophlan full $wd 
#  功能基因diamond索引 10.3G
humann2_databases --download uniref uniref90_diamond $wd
# 设置数据库位置
# 显示参数
humann2_config --print
# 如修改线程数
humann2_config --update run_modes threads 8
humann2_config --update database_folders protein $wd/uniref
humann2_config --update database_folders nucleotide $wd/chocophlan
# metaphlan2数据库默认位于程序所在目录的db_v20和databases下各一份

# 物种注释
# 基于LCA算法的物种注释kraken2  https://ccb.jhu.edu/software/kraken/
conda install kraken2 # 2.0.7_beta
# 下载数据库
kraken2-build --standard --threads 24 --db /db/kraken2
# 此步下载数据>50GB，下载时间由网速决定，索引时间4小时33分，多线程最快35min完成
# 标准模式下只下载5种数据库：古菌archaea、细菌bacteria、人类human、载体UniVec_Core、病毒viral

## 1.3 基因组拼接、注释和定量

# megahit 快速组装
conda install megahit # v1.1.3
# QUEST 组装评估
conda install quast # 5.0.1 

# prokka 细菌基因组注释
conda install prokka # 1.13.3
# cd-hit 非冗余基因集
conda install cd-hit # 4.6.8
# emboss transeq工具翻译核酸为蛋白
conda install emboss # 6.6.0
# 定量工具salmon
conda install salmon # 


## 1.4 基因功能注释

### COG/eggNOG http://eggnogdb.embl.de
# 安装eggnog比对工具
conda install eggnog-mapper
# 下载常用数据库，注意设置下载位置
mkdir -p ~/db/eggnog && cd ~/db/eggnog
download_eggnog_data.py --data_dir ./ -y -f euk bact arch viruses
# 如果内存够大，复制eggNOG至内存加速比对
cp eggnog.db /dev/shm
# 手工整理COG分类注释
wget -c wget http://bailab.genetics.ac.cn/share/COG.anno
# 手工整理KO注释
wget -c wget http://bailab.genetics.ac.cn/share/KO.anno


### 碳水化合物数据库dbCAN2 http://cys.bios.niu.edu/dbCAN2/
mkdir -p ~/db/dbCAN2 && cd ~/db/dbCAN2
wget -c http://cys.bios.niu.edu/dbCAN2/download/CAZyDB.07312018.fa
wget -c http://cys.bios.niu.edu/dbCAN2/download/CAZyDB.07312018.fam-activities.txt
diamond makedb --in CAZyDB.07312018.fa --db CAZyDB.07312018
# 提取fam对应注释
grep -v '#' CAZyDB.07312018.fam-activities.txt|sed 's/  //'| \
  sed '1 i ID\tDescription' > fam_description.txt

### 抗生素抗性基因Resfam http://dantaslab.wustl.edu/resfams
mkdir -p resfam && cd ~/db/resfam
# 官网的数据格式非常混乱, 推荐下载我手动整理的索引和注释
wget http://bailab.genetics.ac.cn/share/Resfams-proteins.dmnd
wget http://bailab.genetics.ac.cn/share/Resfams-proteins_class.tsv


## 1.5 分箱工具

# 物种注释和分箱流程 https://github.com/bxlab/metaWRAP
conda create -n metawrap python=2.7
source activate metawrap
conda config --add channels ursky
conda install -c ursky metawrap-mg 

# 相关数据库，大小近300GB
# 这里我们安装数据库到`~/db`目录，保证你有权限，但要保证至少有500GB的空间。请根据你的情况修改为自己有权限且空间足够的位置。
db=~/db
mkdir -p $db
# 多人使用，建议管理员统一安装节省空间

## CheckM用于Bin完整和污染估计和物种注释
cd $db
mkdir checkm && cd checkm
# 下载文件276MB，解压后1.4GB
wget -c https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz
# 设置数据库位置
checkm data setRoot
# 按提示输出你数据下载的路径

## KRAKEN物种注释数据库
# 下载建索引需要 > 300GB以上空间，完成后占用192GB空间
cd $db
mkdir kraken
kraken-build --standard --threads 24 --db kraken
kraken-build --db kraken --clean

## NCBI_nt核酸序列用于bin物种注释
# 41GB，我下载大约12h；解压后99GB
cd $db
mkdir NCBI_nt && cd NCBI_nt
wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
for a in nt.*.tar.gz; do tar xzf $a; done

## NCBI物种信息
# 压缩文件45M，解压后351M
cd $db
mkdir NCBI_tax
cd NCBI_tax
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz

## 人类基因组去宿主
cd $db
mkdir BMTAGGER_INDEX
cd BMTAGGER_INDEX
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
gunzip *fa.gz
cat *fa > hg38.fa
rm chr*.fa
bmtool -d hg38.fa -o hg38.bitmask
srprism mkindex -i hg38.fa -o hg38.srprism -M 100000

## 数据库位置设置
which config-metawrap
# 查使用vi/vim/gedit等文本编辑器来修改数据库的位置吧

# 退出环境
source deactivate

# 比对结果整理samtools
conda install samtools


# 其它输助脚本

# 按注释类型合并脚本
cd ~/bin
wget http://bailab.genetics.ac.cn/share/mat_gene2ko.R
chmod +x mat_gene2ko.R

# Bin可视化VizBin
sudo apt-get install libatlas3-base libopenblas-base default-jre
curl -L https://github.com/claczny/VizBin/blob/master/VizBin-dist.jar?raw=true > VizBin-dist.jar
mv VizBin-dist.jar /usr/local/bin # 或~/bin


