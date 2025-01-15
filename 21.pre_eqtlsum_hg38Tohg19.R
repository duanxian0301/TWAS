#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/21.pre_eqtlsum_hg38Tohg19")

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
tissue=FUSION_tissue()[c(44),"Tissue"]%>% gsub("GTExv8.ALL.","",.)


#读入基因
sig_gene_GTExv8_ref=data.table::fread(paste0("../19.Venn_UTMOST_FUSION_MAGMA_FOCUS/",trait,"_Venn_UTMOST_FUSION_MAGMA_FOCUS.txt"))


# 保存位置
path <- "../20.pre_eqtlsum_from_GTEx/pre_coloc/" 
path_hg19 <- "./pre_coloc_hg19/"  
if (!dir.exists(path_hg19)) {
  dir.create(path_hg19)
}



#gwas_liftover函数运行时需要下载hg38ToHg19.over.chain.gz文件，注意网速（关闭魔法/科学上网），若运行失败可多试几次或增大timeout参数
options(timeout = 1200)
for (i in 1:length(sig_gene_GTExv8_ref$gene_id)) {
  gene=sig_gene_GTExv8_ref$gene_id[i]
  message("当前第",i,"个基因: ",gene)
  print(sig_gene_GTExv8_ref[i,])
  gwas_liftover(
    gwasfile = paste0(path,gene,"_eqtlsum_from_GTEx_",tissue,"_hg38.txt"),
    build_from = "hg38",
    build_to = "hg19",
    savefile =  paste0(path_hg19,sig_gene_GTExv8_ref$Gene[i],"_",gene,"_eqtlsum_from_GTEx_",tissue,"_hg19.txt")
  )
}


