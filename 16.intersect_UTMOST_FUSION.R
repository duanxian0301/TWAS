

#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/16.intersect_UTMOST_FUSION")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()

#加载包
library(Fast2TWAS)

#输入参数
trait="Prostate cancer"
# 选取GTEx V8中的Whole_Blood组织
tissue=FUSION_tissue()[c(44),"Tissue"] 
P_method = "fdr" #可选择矫正方式
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
P_fdr=0.05

# 求FUSION和UTMOST的交集基因 ---------------------------------------------------
FUSION_file=paste0("../12.sig_FUSION/",trait,".",tissue,".fusion_twas_symble.",P_method,".txt")
UTMOST_fle=paste0("../09.sig_cross_UTMOST/",trait,"_GTEXv8_49_cross_symble.",P_method,".txt")

#读入数据并筛选
data_fusion=data.table::fread(FUSION_file) %>% subset(.,FDR < P_fdr)
data_UTMOST=data.table::fread(UTMOST_fle) %>% subset(.,FDR < P_fdr)

#求交集
intersect_UTMOST_fusion_gene=intersect(data_fusion$ID, data_UTMOST$symble)


#保存交集基因
data.table::fwrite(data.frame(ID=intersect_UTMOST_fusion_gene),
                   file =paste0(trait,"_intersect_UTMOST_fusion_gene_symble.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

#筛选data_fusion中的intersect_UTMOST_fusion_gene基因
sig_FUSION=data_fusion[data_fusion$ID %in% intersect_UTMOST_fusion_gene,]


#保存FUSION结果中的交集基因
data.table::fwrite(sig_FUSION,
                   file =paste0(trait,".",tissue,".fusion_twas_symble_intersect_UTMOST_fusion_",P_method,".sig.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)
