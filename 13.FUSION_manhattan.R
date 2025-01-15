
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/13.FUSION_manhattan")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#加载包
if (!requireNamespace("qqman", quietly = TRUE))install.packages("qqman")
library(Fast2TWAS)
library(qqman)

#输入参数
trait="Prostate cancer"
#曼哈顿图长宽
width=14
height=6
# 选取GTEx V8中的Whole_Blood组织
tissue=FUSION_tissue()[c(44),"Tissue"] 

# 绘制FUSION结果的曼哈顿图------------------------------------------------------
#读入FUSION结果
data=data.table::fread(paste0("../11.FUSION_ensTOsymble/",trait,".",tissue,".fusion_twas_symble.txt"))
sig_p=0.05/nrow(data)

#设置颜色
colorset <- c('#FF0000', '#FFF999', '#2E8a57', '#7FFDDA', '#6595ED', '#000025', '#FF99FF')
#绘图
pdf(file = paste0(trait,".",tissue,".fusion_twas_symble_Manhattan.pdf"), width = width, height = height)
manhattan(
  data,
  chr = "CHR",
  bp = "P0",
  p = "TWAS.P",
  # col = chr_col,
  col=colorset,
  snp = "ID",
  main = "",
  ylab = '-log10(P-value)',
  # ylim = c(0, max(-log10(na.omit(data$TWAS.P)))+1),
  ylim = c(0, 10*ceiling(max(-log10(na.omit(data$TWAS.P)))/10)), 
  genomewideline = -log10(sig_p),                #【红线显示的是bonferroni矫正阈值=0.05/nrow(data)】
  suggestiveline = -log10(0.05),
  annotatePval = sig_p,
  # annotateTop = FALSE,
  logp = TRUE
)
dev.off()

