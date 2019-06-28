# 宏基因组分析流程

**Metagenome analysis pipeline**

**目录说明**

- `denovo1`  # 宏基因组分析流程无参第一版 makefile格式
- `train1903` # 易生信2019年3月培训教程和代码
- `train1906` # 易生信2019年6月培训教程和代码

## 一、分析准备工作

**Data pre-treatment**

    # 系统要求 System: Linux Ubuntu 18.04 / CentOS 7
    # 依赖软件 Sofware: Rstudio server 1.2、KneadData v0.6.1、MetaPhlAn2 v2.7.5、HUMAnN2 v0.11.2 ......
    # Windows/Linux访问Rstudio服务器，推荐使用Chrome浏览器，可选Edge，IE有兼容问题

- fastqc 0.11.5 对所有原始序列进行质量评估
- multiqc 1.4 汇总多样品质量评估报告

### 1.1 了解工作目录和文件

### 1.2 FastQC质量评估(可选)

### 1.3 序列质控和去宿主 Qaulity control and remove host contamination

### 1.4 质控后质量再评估 (可选)

## 二、宏基因组有参分析流程 

**Reference-based Metagenomic Pipeline (HUMAnN2)**

### 2.1 合并质控文件为HUMAnN2输入

### 2.2 HUMAnN2计算物种和功能组成

### 2.3 物种组成表

### 2.4 功能组成分析

### 2.5 GraPhlAn图

### 2.6 LEfSe差异分析和Cladogram

### 2.7 kraken2基于NCBI数据库注释reads层面


## 三、宏基因组无参分析流程

**De novo Metagenomic Pipeline (MEGAHIT)**

### 3.1 拼接 Assembly

### 3.2 基因注释和定量 Gene annotation & quantitfy

### 3.3 功能基因注释

## 四、挖掘单菌基因组/分箱 

**Binning (MetaWRAP)**

* 官方测试数据

### 4.1 运行三种bin软件

### 4.2 Bin提纯

### 4.3 Bin定量

### 4.4 Bin注释
