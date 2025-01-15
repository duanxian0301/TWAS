
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/09.sig_cross_UTMOST")

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

P_method = "fdr" #可选择矫正方式
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

#读入合并之后的cross_tissue_TWAS结果
data_cross=data.table::fread(paste0("../08.UTMOST_ensTOsymble/result/",trait,"_GTEXv8_49_cross_symble.txt"))

#FDR矫正
data_cross$FDR=p.adjust(p = data_cross$p_value, method = P_method)
data_cross=data_cross[order(data_cross$FDR),]

#保存
data.table::fwrite(data_cross,
                   file = paste0(trait,"_GTEXv8_49_cross_symble.",P_method,".txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

#筛选fdr<0.05
data_cross_sig=data_cross[data_cross$FDR<0.05,]
#保存
data.table::fwrite(data_cross_sig,
                   file = paste0(trait,"_GTEXv8_49_cross_symble.",P_method,".sig.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

