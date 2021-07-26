[TOC]

# 附录1：软件和数据库安装

    # 测试环境为Linux Ubuntu 18.04 / CentOS 7

    # 数据库安装位置，默认~/db目录(无需管理权限)，管理员可安装至/db，安装和运行必备环境变量
    db=~/db
    mkdir -p ${db} && cd ${db}
    # Conda软件安装目录，如管理员可能为/conda2，而自己下载可能为~/miniconda2或~/anaconda2
    soft=~/miniconda2


    # 软件管理器miniconda2
    wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    bash Miniconda2-latest-Linux-x86_64.sh
    # 安装时许可协议可打yes，默认目录为~/miniconda2，默认不运行conda直接回车
    conda -V # 查看版本号 4.6.14
    python --version # 2.7.16
    # 安装说明详见：[Nature Method：Bioconda解决生物软件安装的烦恼](https://mp.weixin.qq.com/s/SzJswztVB9rHVh3Ak7jpfA)
    # 添加生物频道，才能找到生物学软件
    conda config --add channels bioconda
    conda config --add channels conda-forge
    # conda默认配置文件为 ~/.condarc ，如不存在使用 conda config --show-sources 查看配置文件位置
    cat ~/.condarc

    # 软件管理器Anaconda(可选)
    # wget -c https://repo.anaconda.com/archive/Anaconda2-2019.03-Linux-x86_64.sh
    # bash Anaconda2-2019.03-Linux-x86_64.sh


    # 添将conda为环境变量，请根据实际情况修改目录位置
    export PATH="${soft}/bin:$PATH"
    
    # 可选创建虚拟环境，防止污染环境变量
    conda create -n metagenome_env python=2.7 # 22.2MB
    # 新版本时替换source为conda
    source activate metagenome_env


    # 并行计算
    sudo apt install parallel
    parallel --version
    
    
## 1.1 质控软件

    # 质量评估软件fastqc
    conda install fastqc # 180 MB
    fastqc -v # FastQC v0.11.8
    
    # 多样品评估报告汇总multiqc
    conda install multiqc # 111 MB
    multiqc --version # multiqc, version 1.7
    
    # 质量控制流程kneaddata
    # conda install kneaddata # kneaddata 最新版 v0.7.2，去宿主后数据极少且双端数据不对应
    # 指定安装0.6.1
    conda install kneaddata=0.6.1 # 175 MB
    kneaddata --version # 0.6.1
    trimmomatic -version # 0.39
    bowtie2 --version # 2.3.5

    # 查看可用数据库
    kneaddata_database
    # 包括人类、小鼠基因组、人类转录组和核糖体RNA
    # human_genome : bmtagger = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_BMTagger_v0.1.tar.gz
    # human_genome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_Bowtie2_v0.1.tar.gz
    # mouse_C57BL : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz
    # human_transcriptome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz
    # ribosomal_RNA : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.1.tar.gz
    
    # 下载人类基因组bowtie2索引3.4 GB至指定数据目录，当然也可选bmtagger索引
    cd ${db}
    mkdir -p ${db}/kneaddata/human_genome
    kneaddata_database --download human_genome bowtie2 ${db}/kneaddata/human_genome
    # # 如下载小鼠基因组用于小鼠宏基因组研究去宿主 2.8 GB
    # mkdir -p ${db}/kneaddata/mouse_genome
    # kneaddata_database --download mouse_C57BL bowtie2 ${db}/kneaddata/mouse_genome
    # # 人类转录组用于人类相关宏转录组去宿主 211 MB
    # mkdir -p ${db}/kneaddata/human_transcriptome
    # kneaddata_database --download human_transcriptome bowtie2 ${db}/kneaddata/human_transcriptome
    # # SILVA用于宏转录组研究去除核糖体 3.6 GB
    # mkdir -p ${db}/kneaddata/ribosomal_RNA
    # kneaddata_database --download ribosomal_RNA bowtie2 ${db}/kneaddata/ribosomal_RNA

    # (可选) 数据库下载慢，可选上方链接直接下载、百度云链接(公众号后台回复：数据库) 或 国内备份链接
    # mkdir -p ${db}/kneaddata/human_genome && cd ${db}/kneaddata/human_genome
    # wget -c http://210.75.224.110/share/meta/Homo_sapiens_Bowtie2_v0.1.tar.gz
    # tar -xvzf Homo_sapiens_Bowtie2_v0.1.tar.gz


## 1.2 有参分析流程MetaPhlAn2、HUMAnN2

    # 安装MetaPhlAn2、HUMAnN2和所有依赖关系
    conda install humann2 # 141M
    humann2 --version # humann2 v0.11.2
    metaphlan2.py -v # MetaPhlAn version 2.7.5 (6 February 2018)
    diamond help #  v0.8.22.84
    
    # 测试流程是否可用
    humann2_test
    
    # 下载数据库
    humann2_databases # 显示可用数据库
    # utility_mapping : full = http://huttenhower.sph.harvard.edu/humann2_data/full_mapping_1_1.tar.gz
    # chocophlan : full = http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan_plus_viral.v0.1.1.tar.gz
    # uniref : uniref90_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_annotated_1_1.tar.gz
    # uniref : uniref50_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref50_annotated_1_1.tar.gz
    cd ${db}
    mkdir -p ${db}/humann2 # 建立下载目录
    # 微生物泛基因组 5.37 GB
    humann2_databases --download chocophlan full ${db}/humann2
    # 功能基因diamond索引 10.3 GB
    humann2_databases --download uniref uniref90_diamond ${db}/humann2

    # 可选：数据库下载慢，可选上方链接直接下载、百度云链接(公众号后台回复：数据库) 或 国内备份链接
    mkdir -p ${db}/humann2/chocophlan && cd ${db}/humann2/chocophlan
    wget -c http://210.75.224.110/share/meta/full_chocophlan_plus_viral.v0.1.1.tar.gz
    tar xvzf full_chocophlan_plus_viral.v0.1.1.tar.gz
    # uniref90和50任选其1，推荐uniref90更全 5.9 GB
    cd ${db}/humann2
    wget -c http://210.75.224.110/share/meta/uniref90_annotated_1_1.tar.gz
    tar xvzf uniref90_annotated_1_1.tar.gz
    # 服务器小于32G内存选uniref50，更省内存更快 2.5 GB
    #wget -c http://210.75.224.110/share/meta/uniref50_annotated_1_1.tar.gz
    #tar xvzf uniref50_annotated_1_1.tar.gz
    # 不要同一文件中有两个文件，会先比90，再比50出现混乱
    
    # 可选：Metaphlan2数据库没成功的手动配置
    # 下载并构建索引
    mkdir -p ${db}/metaphlan2 && cd ${db}/metaphlan2
    wget -c http://210.75.224.110/share/meta/metaphlan2/mpa_v20_m200.tar
    tar xvf mpa_v20_m200.tar
    bzip2 -d mpa_v20_m200.fna.bz2
    bowtie2-build mpa_v20_m200.fna mpa_v20_m200
    # 链接到软件安装目录
    db=~/db
    soft=~/miniconda2
    mkdir -p ${soft}/bin/db_v20
    ln -s ${db}/metaphlan2/* ${soft}/bin/db_v20/
    mkdir -p ${soft}/bin/databases
    ln -s ${db}/metaphlan2/* ${soft}/bin/databases/

    # 设置数据库位置
    # 显示参数
    humann2_config --print
    # 如修改线程数，推荐3-8，根据实际情况调整
    humann2_config --update run_modes threads 3
    humann2_config --update database_folders nucleotide ${db}/humann2/chocophlan
    humann2_config --update database_folders protein ${db}/humann2
    # metaphlan2数据库默认位于程序所在目录的db_v20和databases下各一份
    humann2_config --print


### 1.2.1 物种组成美化GraPhlAn

    # GraPhlAn核心程序包
    conda install graphlan # 22 KB
    graphlan.py --version # GraPhlAn version 1.1.3 (5 June 2018)
    # GraPhlAn输入文件制作程序，如转换LEfSe、Metaphlan2结果格式为GraPhlAn用于绘图
    conda install export2graphlan # 38 KB
    export2graphlan.py -h # ver. 0.20 of 29th May 2017
    
### 1.2.2 物种差异比较和绘制

    conda install lefse # 57.5 MB, 1.0.8.post1

### 1.2.3 物种注释Kraken2

    # 物种注释
    # 基于LCA算法的物种注释kraken2  https://ccb.jhu.edu/software/kraken/
    conda install kraken2 # 2.8 MB
    kraken2 --version # 2.0.8-beta 2019

    # 下载数据库
    cd ${db}
    # --standard标准模式下只下载5种数据库：古菌archaea、细菌bacteria、人类human、载体UniVec_Core、病毒viral
    # 此步下载数据 > 50GB，占用100 GB空间，下载时间由网速决定，索引时间4小时33分，多线程最快35min完成
    kraken2-build --standard --threads 24 --db ${db}/kraken2

    # 个性化定制数据库(可选，详见 https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)


## 1.3 基因组拼接、注释和定量

    # megahit 快速组装
    conda install megahit # 6.4 MB
    megahit -v # v1.1.3
    
    # metaSPAdes拼接，只是spades系列中的一个定制脚本
    conda install spades # 13.7 MB
    metaspades.py -v # SPAdes v3.13.1 [metaSPAdes mode]
    
    # QUEST 组装评估
    conda install quast # 87.2 MB
    metaquast.py -v # QUAST v5.0.2 (MetaQUAST mode)
    
    # prokka 细菌基因组注释
    conda install prokka # 352.8 MB
    prokka -v # 1.13.3
    
    # cd-hit 非冗余基因集
    conda install cd-hit # 790 KB
    cd-hit -v # 4.8.1 (built on May 14 2019)
    
    # emboss transeq工具
    conda install emboss # 94.5 MB
    embossversion # 6.6.0.0
    
    # 定量工具salmon
    conda install salmon # 15.1 MB
    salmon -v # 0.13.1


## 1.4 基因功能注释

    ### COG/eggNOG http://eggnogdb.embl.de
    # 安装eggnog比对工具
    conda install eggnog-mapper # 0.13.1
    emapper.py --version # 1.0.3
    
    # 下载常用数据库，注意设置下载位置
    mkdir -p ${db}/eggnog && cd ${db}/eggnog
    download_eggnog_data.py --data_dir ./ -y -f euk bact arch viruses
    # 如果内存够大，复制eggNOG至内存加速比对
    # cp eggnog.db /dev/shm
    # 手工整理COG分类注释
    wget -c wget http://210.75.224.110/share/COG.anno
    # 手工整理KO注释
    wget -c wget http://210.75.224.110/share/KO.anno
    
    # eggNOG mapper v2 https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2#Installation
    git clone https://github.com/jhcepas/eggnog-mapper.git
    python eggnog-mapper/emapper.py --version # 1.0.3，同上

    
    
    ### 碳水化合物数据库dbCAN2 http://cys.bios.niu.edu/dbCAN2/
    mkdir -p ${db}/dbCAN2 && cd ${db}/dbCAN2
    # 最近此数据库无法访问，改用国内备份链接
    # wget -c http://cys.bios.niu.edu/dbCAN2/download/CAZyDB.07312018.fa
    # wget -c http://cys.bios.niu.edu/dbCAN2/download/CAZyDB.07312018.fam-activities.txt
    wget -c http://210.75.224.110/share/meta/dbcan2/CAZyDB.07312018.fa # 497 MB
    wget -c http://210.75.224.110/share/meta/dbcan2/CAZyDB.07312018.fam-activities.txt # 58 KB
    time diamond makedb --in CAZyDB.07312018.fa --db CAZyDB.07312018 # 28s
    # 提取fam对应注释
    grep -v '#' CAZyDB.07312018.fam-activities.txt|sed 's/  //'| \
      sed '1 i ID\tDescription' > fam_description.txt
    
    ### 抗生素抗性基因Resfam http://dantaslab.wustl.edu/resfams
    mkdir -p ${db}/resfam && cd ${db}/resfam
    # 官网的数据格式非常混乱, 推荐下载我手动整理的索引和注释
    wget http://210.75.224.110/share/Resfams-proteins.dmnd # 1.5 MB
    wget http://210.75.224.110/share/Resfams-proteins_class.tsv # 304 KB


## 1.5 分箱工具

    # 物种注释和分箱流程 https://github.com/bxlab/metaWRAP
    conda create -n metawrap python=2.7 # 22.2MB
    conda activate metawrap
    conda config --add channels ursky
    conda install -c ursky metawrap-mg # 1.14 GB, v1.2
    
    conda list # 显示软件列表
    conda env list # 显示环境列表
    kraken -v # 查看软件版本 1.1.1
    # conda deactivate # 退出环境
    
    # 相关数据库，大小近300GB
    # 这里我们安装数据库到`~/db`目录，保证你有权限，但要保证至少有500GB的空间。请根据你的情况修改为自己有权限且空间足够的位置。
    # 多人使用，建议管理员统一安装节省空间
    cd ${db}

    ## CheckM用于Bin完整和污染估计和物种注释
    mkdir -p checkm && cd checkm
    # 下载文件275 MB，解压后1.4 GB
    wget -c https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    # 国内备用链接
    # wget -c http://210.75.224.110/share/meta/checkm/checkm_data_2015_01_16.tar.gz
    tar -xvf *.tar.gz
    # rm *.gz
    # 设置数据库位置
    checkm data setRoot
    # 按提示输出你数据下载的路径或直接回车默认为当前位置
    
    ## NCBI_nt核酸序列用于bin物种注释
    # 41GB，我下载大约12h；解压后99GB
    cd ${db}
    mkdir -p NCBI_nt && cd NCBI_nt
    wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
    # 备用下载链接，或百度云下载
    # wget -c http://210.75.224.110/share/meta/NCBI_nt/filelist.txt
    # for a in `cat filelist.txt`; do wget -c http://210.75.224.110/share/meta/NCBI_nt/$a; done
    for a in nt.*.tar.gz; do tar xzf $a; done &
    
    ## NCBI物种信息
    # 压缩文件45M，解压后351M
    cd ${db}
    mkdir NCBI_tax
    cd NCBI_tax
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xvf taxdump.tar.gz
    
    ## 人类基因组去宿主
    cd ${db}
    mkdir -p metawrap/BMTAGGER && cd metawrap/BMTAGGER
    wget -c ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
    gunzip *fa.gz
    cat *fa > hg38.fa
    rm chr*.fa
    # 上方下载太慢，使用国内备份链接手动下载
    wget -c http://210.75.224.110/share/meta/metawrap/BMTAGGER/hg38.fa
    bmtool -d hg38.fa -o hg38.bitmask
    srprism mkindex -i hg38.fa -o hg38.srprism -M 100000
    
    ## KRAKEN物种注释数据库
    # 下载建索引需要 > 300GB以上空间，完成后占用192GB空间
    cd ${db}
    mkdir -p kraken
    kraken-build --standard --threads 24 --db kraken > log &
    kraken-build --db kraken --clean
    # 手动下载
    cd kraken
    wget -c http://210.75.224.110/share/meta/kraken/database.kdb
    wget -c http://210.75.224.110/share/meta/kraken/database.idx
    mkdir -p taxonomy && cd taxonomy
    wget -c http://210.75.224.110/share/meta/kraken/taxonomy/nodes.dmp
    wget -c http://210.75.224.110/share/meta/kraken/taxonomy/names.dmp
    # 从其它位置复制
    # cp -r /db/kraken/* ./
    
    ## 数据库位置设置
    which config-metawrap
    # 配置文件通常为~/miniconda2/envs/metawrap/bin/config-metawrap
    # 使用Rstudio/vim等文本编辑器来修改数据库的位置
    
    # 软件安装结果时，建议下载Krona物种数据库，个性化数据库位置的方法
    # 创建目标默认目录为软链，则可搬移数据库位置
    # rm -rf ~/miniconda2/envs/metawrap/opt/krona/taxonomy
    # mkdir -p ~/db/krona/taxonomy
    # ln -s ~/db/krona/taxonomy ~/miniconda2/envs/metawrap/opt/krona/taxonomy
    # ktUpdateTaxonomy.sh

    # QUAST评估默认没下载数据库，下载可增加评估报告的信息
    # - The default QUAST package does not include:
    # * GRIDSS (needed for structural variants detection)
    # * SILVA 16S rRNA database (needed for reference genome detection in metagenomic datasets)
    # * BUSCO tools and databases (needed for searching BUSCO genes) -- works in Linux only!
    # 下载QUAST依赖数据库
    # quast-download-gridss # 38 MB，~/miniconda2/envs/metawrap/lib/python2.7/site-packages/quast_libs/gridss/gridss-1.4.1.jar
    # quast-download-silva # 177 MB
    # quast-download-busco

  
    
## 1.6 其它软件和输助脚本
    
    # 按注释类型合并脚本，链接 https://github.com/YongxinLiu/Metagenome
    cd ${db}
    mkdir script && cd script
    # 矩阵按最后列注释合并的脚本，如Gene-KO-Module-Pathway的合并
    wget http://210.75.224.110/share/meta/script/mat_gene2ko.R
    # 基于Kraken2的结果计算alpha多样性
    wget http://210.75.224.110/share/meta/script/kraken2alpha.R
    # 基于alpha多样性和分组信息绘制箱线图
    wget http://210.75.224.110/share/meta/script/alpha_boxplot.R
    chmod +x *.R

    
    # Bin可视化VizBin
    # sudo apt-get install libatlas3-base libopenblas-base default-jre
    # curl -L https://github.com/claczny/VizBin/blob/master/VizBin-dist.jar?raw=true > VizBin-dist.jar
    # mv VizBin-dist.jar /usr/local/bin # 或~/bin
    
    # 比对结果整理samtools
    conda install samtools

    ### CARD(选学) https://card.mcmaster.ca/download 
    # 方法1. 直接conda安装
    conda install --channel bioconda rgi
    # 如果报错，尝试方法2。CondaMultiError: CondaFileIOError: '/home/liuyongxin/miniconda2/pkgs/prokka-1.11-0.tar.bz2'. contains unsafe path: db/cm/READM
    # 方法2. 虚拟环境安装
    conda activate metawrap
    conda install --channel bioconda rgi
    rgi main -v # 4.0.3
    # rgi教程 https://github.com/arpcard/rgi
    
  

## 附录

