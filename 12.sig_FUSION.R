
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/12.sig_FUSION")

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
# 选取GTEx V8中的Whole_Blood组织,可调整数值，改换其他组织，但是要另行下载数据
tissue=FUSION_tissue()[c(44),"Tissue"] 
P_method = "fdr" #可选择矫正方式

# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

#读入FUSION结果
data_single=data.table::fread(paste0("../11.FUSION_ensTOsymble/",trait,".",tissue,".fusion_twas_symble.txt"))

#FDR矫正------------------------------------------------------------------------
data_single$FDR=p.adjust(p = data_single$TWAS.P, method = P_method)
data_single=data_single[order(data_single$FDR),]
#保存
data.table::fwrite(data_single,
                   file =paste0(trait,".",tissue,".fusion_twas_symble.",P_method,".txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

#筛选fdr<0.05
data_single_sig=subset(data_single,FDR<0.05)
#保存
data.table::fwrite(data_single_sig,
                   file =paste0(trait,".",tissue,".fusion_twas_symble.",P_method,".sig.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)


