
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/19.Venn_UTMOST_FUSION_MAGMA_FOCUS")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#加载包
if (!requireNamespace("VennDiagram", quietly = TRUE))install.packages("VennDiagram")
library(VennDiagram) 
library(Fast2TWAS)

#输入参数
trait="Prostate cancer"
P_method = "fdr" #可选择矫正方式
# 选取GTEx V8中的Whole_Blood组织
tissue=FUSION_tissue()[c(44),"Tissue"]
#维恩图长宽
width=7
height=7

# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

#读入数据
FUSION_file=paste0("../12.sig_FUSION/",trait,".",tissue,".fusion_twas_symble.",P_method,".sig.txt")
UTMOST_fle=paste0("../09.sig_cross_UTMOST/",trait,"_GTEXv8_49_cross_symble.",P_method,".sig.txt")
FOCUS_fle=paste0("../18.FOCUS/",trait,"_",tissue,"_focus.sig.txt")
MAGMA_fle=paste0("../15.plot_MAGMA/",trait,"_magma_gene_based_plot_symble.",P_method,".sig.txt")

data_fusion=data.table::fread(FUSION_file) 
data_UTMOST=data.table::fread(UTMOST_fle) 
data_FOCUS=data.table::fread(FOCUS_fle) 
data_MAGMA=data.table::fread(MAGMA_fle) 


#创建list
ven_list=list(FUSION=unique(data_fusion$ID),
              UTMOST=unique(data_UTMOST$symble),
             FOCUS=unique(data_FOCUS$mol_name),
              MAGMA=unique(data_MAGMA$SYMBOL) )

#四元#
#e8736b/#3a7c38/#459093/#e8c2aa
venn.plot=venn.diagram(ven_list,
                       filename=NULL,
                       main="Venn_UTMOST_FUSION_MAGMA_FOCUS",
                       main.cex = 2,
                       fill=c('#e8736b',
                              '#3a7c38',
                              '#e8c2aa',
                              '#459093' ),
                       cat.cex=1.2)

pdf(file=paste0(trait,"_Venn_UTMOST_FUSION_MAGMA_FOCUS.pdf"),width=height,height=height)
grid.draw(venn.plot)
dev.off()


#指定GTExv8基因转换参考文件
GTExv8_ref="../03.reference/GTEX_V8_gene_info.txt"
#读入参考文件
GTEX_V8_gene_id=data.table::fread(GTExv8_ref)

#交集基因
intersect_Genes=Reduce(intersect,ven_list) 

#合并
sig_gene_GTExv8_ref=data.frame(merge(data.frame(Gene=intersect_Genes),GTEX_V8_gene_id,by.x = "Gene", by.y = "symble",all.x = T)) 

#保存
data.table::fwrite(sig_gene_GTExv8_ref,
                   file =paste0(trait,"_Venn_UTMOST_FUSION_MAGMA_FOCUS.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

