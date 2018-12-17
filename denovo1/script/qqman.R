#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录
rm(list = ls())



# 1.1 程序功能描述和主要步骤

# 程序功能：基于gemma结果绘制manhattan plot
# Functions: GWAS, gemma result plot
# Main steps: 
# - Reads gemma pvaule
# - Draw manhattanplot and save in output.pdf



# 1.2 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("optparse", "reshape2","ggplot2","devtools","dplyr","qqman") # 
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 1.3 解析命令行

# 解析参数-h显示帮助信息
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="tassel/alpha_richness.txt",
                help="Input table file to read; Gemma结果pvalu文件 [default %default]"),
    make_option(c("-t", "--type"), type="numeric", default="4",
                help="column of pvalue; P值所在列 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="doc/design.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=16,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=4.5,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 有统计表txt和矢量图pdf [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$output==""){
    opts$output=paste(opts$input, ".pdf", sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("Column of pvalue is ", opts$type,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


# 2. 统计与绘图
gwas <- read.table(opts$input, header=T, check.names=F)

sub_gwas <- subset(gwas, gwas$P <=0.001,select = c(1,2,3,opts$type))

# colnames(sub_table) = c("CHR", "SNP", "BP", "P")

pdf(file=opts$output, width=opts$width, height=opts$height) # output to PDF or screen

manhattan(sub_gwas, main = "Manhattan Plot", ylim = c(3, 10), cex = 0.6, 
          cex.axis = 0.9, col = c("blue4", "orange3"))
dev.off()


