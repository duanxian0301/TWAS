
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/17.Joint_conditional")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#加载包
library(Fast2TWAS)

#输入参数
trait="Prostate cancer"
# 选取GTEx V8中的Whole_Blood组织
tissue=FUSION_tissue()[c(44),"Tissue"] 
P_method = "fdr" #可选择矫正方式
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")


# Conditional and joint analysis -----------------------------------------------
data=TWAS_fusion_conditional_plot_multiprocess(sigtwasfile = paste0("../16.intersect_UTMOST_FUSION/",trait,".",tissue,".fusion_twas_symble_intersect_UTMOST_fusion_",P_method,".sig.txt"),
                                               sumstatsfile = paste0("../10.FUSION/",trait,".sumstats.gz"),
                                               LDREF = "../03.reference/FUSION/LDREF/1000G.EUR",
                                               outname = trait,
                                               nthread = 4)

#保存
data.table::fwrite(data,
                   file =paste0(trait,".",tissue,".fusion_twas_symble_intersect_UTMOST_fusion_",P_method,".sig.Jointgenes.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

# ------------------------------------------------------------------------------
file_list=list.files("./result/Joint_condition_plot/",
                     pattern = ".joint_included.dat", 
                     full.names = T)

# 使用lapply()读取文件
data_list <- lapply(file_list, data.table::fread)

# 如果需要合并所有数据，你可以使用do.call()和rbind()（假设所有文件具有相同的列名）
combined_data <- do.call(rbind, data_list)


#保存
data.table::fwrite(combined_data,
                   file =paste0(trait,".",tissue,".fusion_twas_symble_intersect_UTMOST_fusion_",P_method,".sig.Joint.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

