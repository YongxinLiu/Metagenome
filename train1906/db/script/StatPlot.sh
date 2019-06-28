# 统计绘图

cd /c/meta

## 2.3 HUMAnN2

### 2.3.1 物种组成热图

    # 方法1. metaphlan_hclust_heatmap.py 脚本服务器绘制热图，依赖关系和环境变量复杂
    # 方法2. Excel筛选 metaphlan2/taxonomy.tsv 并在线绘制热图
    # 方法3. R数据筛选taxonomy.spf并用pheatmap绘制热图
    
    # 显示脚本帮助 help
    Rscript db/script/metaphlan_hclust_heatmap.R -h
    # 默认读取result/metaphlan2/taxonomy.spf，按指定列合并、排序并取Top25种绘制热图，输出至输入目录
    Rscript db/script/metaphlan_hclust_heatmap.R
    # 完整参数：-i输入MetaPhlAn2文件；
    # -t 分类级别，可选Kingdom/Phylum/Class/Order/Family/Genus/Species/Strain，界门纲目科属种株，推荐门，目，属
    # -n 输出物种数量，默认为25，最大为合并后的数量
    # -o输出图表前缀，默认根据输入文件、物种级别和数量自动生成；
    Rscript db/script/metaphlan_hclust_heatmap.R \
      -i result/metaphlan2/taxonomy.spf \
      -t Order \
      -n 30 \
      -o result/metaphlan2/heatmap_Order

## 2.7 kraken2基于NCBI数据库注释reads层面


### 2.7.3 物种多样性分析

    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # -d指定最小样本量，不指定为最小值，抽平文件为taxonomy_count.norm.txt，alpha多样性/taxonomy_count.alpha.txt
    # Rscript db/script/kraken2alpha.R -h # 查看参数
    Rscript db/script/kraken2alpha.R -i result/kraken2/taxonomy_count.txt
    
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson/invsimpson
    # Rscript db/script/alpha_boxplot.R -h # 查看参数
    Rscript db/script/alpha_boxplot.R -i result/kraken2/taxonomy_count.alpha.txt -t shannon \
      -d result/metadata.txt -n group -w 4 -e 2.5
    # 批量计算6种指数的箱线图+统计
    for i in richness chao1 ACE shannon simpson invsimpson;do
    Rscript db/script/alpha_boxplot.R -i result/kraken2/taxonomy_count.alpha.txt -t ${i} \
      -d result/metadata.txt -n group -w 4 -e 2.5
    done
    # 解读：richness和shannon无显著差异，chao1有显著差异，证明测序深度不足(测试数据仅有0.25%数据)


### 2.7.4 物种组成

    # 物种组成热图
    # 转换为metaphalan2 spf格式，但ncbi注释不完整
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["d"]="unclassified";a["p"]="unclassified";a["c"]="unclassified";a["o"]="unclassified";a["f"]="unclassified";a["g"]="unclassified";a["s"]="unclassified";a["S"]="unclassified"; \
      split($1,x,"|");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print a["d"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"],a["S"],$0;}' \
      result/kraken2/taxonomy_count.norm.txt > result/kraken2/temp.txt
      cut -f 1-8,10- result/kraken2/temp.txt > result/kraken2/taxonomy_count.norm.spf
      sed -i '1 s/unclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified\tunclassified/Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain/' result/kraken2/taxonomy_count.norm.spf
    # 绘制热图
    Rscript db/script/metaphlan_hclust_heatmap.R \
      -i result/kraken2/taxonomy_count.norm.spf \
      -t Genus \
      -n 25 \
      -o result/kraken2/heatmap_Genus

    # 绘制属水平Top30箱线图
    Rscript db/script/metaphlan_boxplot.R \
      -i result/kraken2/taxonomy_count.norm.spf \
      -t Genus \
      -n 30 \
      -o result/kraken2/boxplot_Genus
    # 绘制门水平Top10箱线图
    Rscript db/script/metaphlan_boxplot.R \
      -i result/kraken2/taxonomy_count.norm.spf \
      -t Phylum \
      -n 10 -w 6 -e 4 \
      -o result/kraken2/boxplot_Phylum
  


# 附录

# 命令行正常对照 
    # Rscript /c/amplicon/34Evolution/table2itol/table2itol.R -h
    
    dos2unix db/script/metaphlan_hclust_heatmap.R    # cat -A db/script/metaphlan_hclust_heatmap.R | head
    