#!/usr/bin/env Rscript
# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo & Yang Bai. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 41, 1-16, doi:10.1007/s13238-020-00724-8 (2020).
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 1. 分析前准备：帮助、参数、依赖包和读取文件

# 命令行运行为当前目录；Rstudio手动运行脚本前需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录


# 1.1 程序功能描述和主要步骤

# 程序功能：Alpha多样性箱线图绘制及组间统计
# Functions: Boxplot show alpha richness among groups
# Main steps: 
# - Reads data table input.txt
# - Calculate pvalue and save in output.txt
# - Draw boxplot and save in output.pdf

# 程序使用示例 USAGE
# Default
# # 显示脚本帮助
# Rscript ./script/alpha_boxplot.R -h
# # 默认画richness
# Rscript ./script/alpha_boxplot.R
# # 完整参数，输出文件名默认为alpha指数类型 -t 指数列名; -d 实验设计文件;-n 分组列名; -w  图片宽;-e 图片高; -o 输出文件前缀, 有统计表txt和矢量图pdf
# Rscript ./script/alpha_boxplot.R -i result/kraken2/tax_alpha.txt -t richness \
# -d result/metadata.txt -n group \
# -o result/kraken2/alpha_richness \
# -w 4 -e 2.5


# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/kraken2/tax_alpha.txt",
                help="Input alpha table [default %default]"),
    make_option(c("-t", "--type"), type="character", default="richness",
                help="type of alpha index [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="design file [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="name of group type [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=4,
                help="Width of figure [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=2.5,
                help="Height of figure [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix, include pdf and txt [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  prefix = gsub(".txt$", "", opts$input, perl = T)
  if (opts$output==""){
    opts$output=paste(prefix, "_", opts$type, sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("Type of alpha index is ", opts$type,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}


# 2.1 安装CRAN来源常用包
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
# c("reshape2","ggplot2","devtools","bindrcpp","ggthemes","agricolae","dplyr","multcompView")
package_list = c("dplyr","ggplot2","agricolae") 
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 将Tukey检验结果P值转换为显著字母分组
generate_label_df = function(TUKEY, variable){
  library(multcompView)
  # 转换P值为字母分组
  ## 提取图基检验中分组子表的第4列P adjust值
  Tukey.levels = TUKEY[[variable]][,4]
  ## multcompLetters函数将两两p值转换为字母，data.frame并生成列名为Letters的数据框
  Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # 按分组名字母顺序
  ## 提取字母分组行名为group组名
  Tukey.labels$group = rownames(Tukey.labels)
  # 按组名的字线顺序排列，默认的Levels
  Tukey.labels=Tukey.labels[order(Tukey.labels$group) , ]
  return(Tukey.labels)
}



# 3. 读取输入文件

# 读取usearch alpha文件
alpha = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="") 

# 提取样品组信息,默认为group可指定
sampFile = as.data.frame(design[,opts$group],row.names = row.names(design))
colnames(sampFile)[1] = "group"

## richness index
# add design to alpha
#index = cbind(alpha[rownames(design),]$richness, sampFile) 
index = cbind(alpha[rownames(sampFile),][[opts$type]], sampFile) 
colnames(index) = c(opts$type,"group") # add richness colname is value



# 4. 统计与绘图

# 统计各组间差异

#model = aov(richness ~ group, data=index)
model = aov(index[[opts$type]] ~ group, data=index)

# 计算Tukey显著性差异检验
Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
# 提取比较结果
Tukey_HSD_table = as.data.frame(Tukey_HSD$group) 

# LSD检验，添加差异组字母
out = LSD.test(model,"group", p.adj="none")
stat = out$groups
# 分组结果添入Index
index$stat=stat[as.character(index$group),]$groups

# 设置分组位置为各组y最大值+高的3%
max=max(index[,c(opts$type)])
min=min(index[,opts$type])
x = index[,c("group",opts$type)]

# 下方的richness如何替换为变量
y = x %>% group_by(group) %>% summarise_(Max=paste('max(',opts$type,')',sep=""))

y=as.data.frame(y)
rownames(y)=y$group
index$y=y[as.character(index$group),]$Max + (max-min)*0.05

p = ggplot(index, aes(x=group, y=index[[opts$type]], color=group)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste(opts$type, "index")) + theme_classic() +
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)
# p



# 5. 保存图表

# 保存pdf和png格式方便查看和编辑
suppressWarnings(ggsave(paste(opts$output, ".pdf", sep=""), p, width = opts$width, height = opts$height))
# ggsave(paste(opts$output, ".png", sep=""), p, width = opts$width, height = opts$height)

# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
suppressWarnings(write.table(Tukey_HSD_table, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

# 提示工作完成
print(paste("Output in ", opts$output, ".txt/pdf finished.", sep = ""))
