
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/15.plot_MAGMA")

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
library(ggplot2)
library(qqman)

#输入参数
trait="Prostate cancer"
P_method = "fdr" #可选择矫正方式

# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

#设置MAGMA Tissue Expression Analysis阈值
#曼哈顿图长宽
width=14
height=6
#条形图长宽
width_bar=8
height_bar=6

# gene_ensembl TO gene_symbol---------------------------------------------------
inputfile=paste0("../14.MAGMA/result/",trait,"_magma.genes.out.txt")
#读入数据
data=data.table::fread(inputfile)
data=subset(data,CHR!="X")
#读入ENSGv110.coding.genes.txt文件
reference=data.table::fread("../03.reference/magma/NCBI37/ENSGv110.coding.genes.txt")

#合并数据
mergedata=merge(data,reference,by.x = "GENE", by.y = "V1")
mergedata=mergedata[,c(1:9,14)]
mergedata=mergedata[order(mergedata$P),]
colnames(mergedata)[10]="SYMBOL"
mergedata=mergedata[order(mergedata$P),]

#保存数据
data.table::fwrite(mergedata,
                   file=paste0(trait,"_magma_gene_based_plot_symble.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

#FDR矫正------------------------------------------------------------------------
mergedata$FDR=p.adjust(p = mergedata$P, method = P_method)

#保存
data.table::fwrite(mergedata,
                   file =paste0(trait,"_magma_gene_based_plot_symble.",P_method,".txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

#筛选fdr<0.05
mergedata_sig=subset(mergedata,FDR<0.05)
#保存
data.table::fwrite(mergedata_sig,
                   file =paste0(trait,"_magma_gene_based_plot_symble.",P_method,".sig.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

#绘制MAGMA的gene_based曼哈顿图--------------------------------------------------
#设置颜色
colorset <- c('#FF0000', '#FFD700', '#2E8B57', '#7FFFAA', '#6495ED', '#0000FF', '#FF00FF')
#绘图
mergedata$CHR=as.numeric(mergedata$CHR)
pdf(file = paste0(trait,"_magma_gene_based_plot_symble.pdf"),width = width, height = height)
manhattan(
  mergedata,
  chr = "CHR",
  bp = "START",
  p = "P",
  # col = chr_col,
  col=colorset,
  snp = "SYMBOL",
  main = "",
  ylab = '-log10(P-value)',
  ylim = c(0, 10*ceiling(max(-log10(na.omit(mergedata$P)))/10)),
  genomewideline = -log10(0.05/nrow(mergedata)),         #【红线显示的是bonferroni矫正阈值=0.05/nrow(data)】
  suggestiveline = -log10(0.05),
  annotatePval = 0.05/nrow(mergedata),
  # annotateTop = FALSE,
  logp = TRUE
)
dev.off()


#绘制MAGMA的gene_set图 ---------------------------------------------------------
data2=data.table::fread(paste0("../14.MAGMA/result/",trait,"_magma.gsa.out.txt"))
#FDR矫正
data2$FDR=p.adjust(p = data2$P, method = P_method)
# 计算Negative Log10 FDR
data2$logFDR <- -log10(data2$FDR)
# 按照FDR从小到大对数据框进行排序
data2 <- data2[order(data2$FDR),]

# num <- nrow(subset(data2,FDR<0.05))
# if (num > 50) {
#   data_num <- data2[c(1:50),]
# }else{
#   data_num=data_num[order(data2$FDR),]%>% subset(.,FDR<0.05)
# }

#修改为只画前50条通路
data_num<-data2[c(1:50),]

threshold=0.05
# 使用ggplot2创建柱状图
p <- ggplot(data_num, aes(x = reorder(VARIABLE, -logFDR), y = logFDR, fill = ifelse(logFDR > -log10(threshold), "red", "blue"))) +
  geom_col() +
  scale_fill_identity(guide = FALSE) +
  labs(x = "", y = "-log10(FDR)",
       title = "                                                MAGMA Gene Pathway Analysis",
       subtitle = paste("FDR < ", 0.05)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "black")

#保存pdf
ggsave(paste0(trait,"_magma_gene_pathway.pdf"), plot = p,width = width_bar, height = height_bar)

#保存
data.table::fwrite(data2,
                   file =paste0(trait,"_magma_gene_pathway.",P_method,".txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

data2_sig=subset(data2,FDR<0.05)
data.table::fwrite(data2_sig,
                   file =paste0(trait,"_magma_gene_pathway.",P_method,".sig.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)


#绘制MAGMA组织特异性富集图------------------------------------------------------
data1=data.table::fread(paste0("../14.MAGMA/result/",trait,"_tissue_specific_magma_result.gsa.out.txt"))
#FDR矫正
data1$FDR=p.adjust(p = data1$P, method = P_method)
# 计算Negative Log10 FDR
data1$logFDR <- -log10(data1$FDR)
# 按照FDR从小到大对数据框进行排序
data1 <- data1[order(data1$FDR),]


threshold=0.05
# 使用ggplot2创建柱状图
p1 <- ggplot(data1, aes(x = reorder(FULL_NAME, -logFDR), y = logFDR, fill = ifelse(logFDR > -log10(threshold), "red", "blue"))) +
  geom_col() +
  scale_fill_identity(guide = FALSE) +
  labs(x = "", y = "-log10(FDR)",
       title = "                                                MAGMA Tissue Expression Analysis",
       subtitle = paste("FDR < ", 0.05)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "black")

#保存pdf
ggsave(paste0(trait,"_magma_GTEx54_tissue_specific.pdf"), plot = p1,width = width_bar, height = height_bar)

#保存
data.table::fwrite(data1,
                   file =paste0(trait,"_magma_GTEx54_tissue_specific.",P_method,".txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

data1_sig=subset(data1,FDR<0.05)
#保存
data.table::fwrite(data1_sig,
                   file =paste0(trait,"_magma_GTEx54_tissue_specific.",P_method,".sig.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)


