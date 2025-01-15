
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/24.plotMR")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#加载包
if (!requireNamespace("forestploter", quietly = TRUE))install.packages("forestploter")
library(forestploter)
library(ggplot2)
library(Fast2TWAS)

#参数设置
tissue=FUSION_tissue()[c(44),"Tissue"]%>% gsub("GTExv8.ALL.","",.) # 选取GTEx V8中的Whole_Blood组织
trait="Prostate cancer"
#图片的长宽
width=11
height=9

#读入基因
intersect_gene=data.table::fread(paste0("../19.Venn_UTMOST_FUSION_MAGMA_FOCUS/",trait,"_Venn_UTMOST_FUSION_MAGMA_FOCUS.txt"))

#构建数据文件-------------------------------------------------------------------
path="../23.mendelian_randomization/"

#循环
all_data=data.frame()
for (i in 1:length(intersect_gene$gene_id)) {
  gene=intersect_gene$Gene[i]
  message("当前第",i,"个基因: ",gene)
  print(intersect_gene[i,])
  #判断是否存在.Mendelian_IVs.txt文件
  iv_file=paste0(path,gene,"_",trait,".Mendelian_IVs.txt")
  if (file.exists(iv_file)) {
    data1=data.table::fread(paste0(path,gene,"_",trait,".Mendelian_IVs.txt"))
  }else{
    message(gene," 与 ",paste0("../06.pre_data/",trait,".txt"),"进行孟德尔随机化不成功，直接跳过 ",gene," 基因")
    next}
  
  # 判断工具变量的SNPs数量
  if (nrow(data1)==1) {
    data2=data.table::fread(paste0(path,gene,"_",trait,"_final_MRresult.txt"))
    mr_data=data.frame(
      Gene=gene,
      method="Wald ratio",
      Instruments=data1$SNP,
      b=data2$b,
      se=data2$se,
      pval=data2$pval,
      b_lo_ci=data2$lo_ci,
      b_up_ci=data2$up_ci,
      or=data2$or,
      or_lci95=data2$or_lci95,
      or_uci95=data2$or_uci95
    )
    all_data=rbind(mr_data,all_data)
  }else if(nrow(data1)>1){
    data2=data.table::fread(paste0(path,gene,"_",trait,"_final_MRresult.txt"))
    mr_data=data.frame(
      Gene=gene,
      method="Inverse variance weighted",
      Instruments=data2[data2$Method=="Inverse variance weighted",]$SNP,
      # Instruments=paste0(unlist(strsplit(data2[data2$Method=="Inverse variance weighted",]$SNP, "; ")),collapse = "\n"),
      b=data2[data2$Method=="Inverse variance weighted",]$b,
      se=data2[data2$Method=="Inverse variance weighted",]$se,
      pval=data2[data2$Method=="Inverse variance weighted",]$pval,
      b_lo_ci=data2[data2$Method=="Inverse variance weighted",]$lo_ci,
      b_up_ci=data2[data2$Method=="Inverse variance weighted",]$up_ci,
      or=data2[data2$Method=="Inverse variance weighted",]$or,
      or_lci95=data2[data2$Method=="Inverse variance weighted",]$or_lci95,
      or_uci95=data2[data2$Method=="Inverse variance weighted",]$or_uci95
    )
    all_data=rbind(mr_data,all_data)
  }else{
    stop("运行错误，请检查第23步孟德尔随机化")
  }
}  

#保存数据
#保存
data.table::fwrite(all_data,
                   file =paste0(trait,"_",tissue,"_UTMOST_FUSION_MAGMA_FOCUS_MR_result.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)



# 绘制森林图 -------------------------------------------------------------------
for (i in 1:nrow(all_data)) {
  all_data$Instruments[i]=paste0(unlist(strsplit(all_data$Instruments[i], "; ")),collapse = "\n")
}

data <- all_data[,c(1:6)]
#计算置信区间
data$low=data$b - data$se * 1.96
data$hi=data$b + data$se * 1.96
#添加绘制森领图空白行
data$` ` <- paste(rep(" ", 12), collapse = " ")
data$`OR (95% CI)` <- ifelse(is.na(data$se), "",
                             sprintf("%.2f (%.2f to %.2f)",
                                     data$b, data$low, data$hi))

data$pval<-ifelse(data$pval<0.05,"< 0.05",sprintf("%.3f",data$pval))

#设置主题
tm <- forest_theme(base_size = 15, # 基础大小
                   
                   # 可信区间点的形状，线型、颜色、宽度
                   ci_pch = 16,
                   ci_col = "#4575b4", # #762a83
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # 可信区间两端加短竖线
                   
                   # 参考线宽度、形状、颜色
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   
                   # 汇总菱形的填充色和边框色
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   
                   # 脚注大小、字体、颜色
                   footnote_cex = 0.6,
                   footnote_fontface = "italic",
                   footnote_col = "blue")

tm <- forest_theme(core = list(bg_params=list(fill = c("white"))),#森林图背景填充
                   base_size = 15,#主题字符大小
                   base_family = "Arial",#字体类型
                   summary_col = "#377EB8",#合并置信区间（菱形）的颜色
                   refline_col = "black",#无效线颜色
                   refline_lwd =2,#无效线粗细
                   ci_col = "#F781BF",#置信区间颜色
                   arrow_fill ="#F781BF",#箭头填充颜色
                   arrow_lwd =2,#箭头粗细
                   ci_lwd =2,#置信区间粗细
                   ci_Theight =unit(0.2,'cm'),#置信区间两侧竖线高度
                   arrow_label_just = "end",#箭头文本与箭头对齐方式
                   arrow_type = "closed")#箭头类型

p <- forest(data[,c(1:3,9:10,6)],
            est = data$b ,
            lower = data$low,
            upper = data$hi,
            sizes = data$se,
            ci_column = 4,
            ref_line = 0.05,
            xlim = c(min(data$low), max(data$hi)),
            theme = tm)
# Print plot
plot(p)

p1 <- forest(data = data[,c(1:3,9:10,6)],
             lower = data$low,
             upper = data$hi,
             est = data$b,
             ci_column = 4,
             
             #is_summary = c(rep(FALSE,nrow(tabletext)-1), TRUE), # 最后一列是汇总行
             ref_line = 0, # 把竖线放到1的位置
             xlim = c(-0.4,0.5), # x轴范围
             ticks_at = c(-0.35,0,0.2,0.4), # x轴刻度显示
             theme = tm
)

print()

#保存pdf
ggsave('MR', plot = p1,width = width, height = height)

