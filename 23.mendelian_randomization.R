

#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/23.mendelian_randomization")

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

#参数设置
tissue=FUSION_tissue()[c(44),"Tissue"]%>% gsub("GTExv8.ALL.","",.) # 选取GTEx V8中
trait="Prostate cancer"
clump_kb=10000
r2=0.001

#读入基因
intersect_gene=data.table::fread(paste0("../19.Venn_UTMOST_FUSION_MAGMA_FOCUS/",trait,"_Venn_UTMOST_FUSION_MAGMA_FOCUS.txt"))

#循环进行MR分析
for (i in 1:length(intersect_gene$gene_id)) {
  gene=intersect_gene$Gene[i]
  message("当前第",i,"个基因: ",gene)
  print(intersect_gene[i,])
  genechr=paste0(intersect_gene$Gene[i],"_",intersect_gene$gene_id[i])
  mr_result=mendelian_randomization(exposure_dat=paste0("../21.pre_eqtlsum_hg38Tohg19/pre_coloc_hg19/",genechr,"_eqtlsum_from_GTEx_",tissue,"_hg19.txt"),
                                    outcome_dat = paste0("../06.pre_data/",trait,".txt"),
                                    pval = 5e-8,
                                    exposure_trait = gene,
                                    outcome_trait = trait,
                                    clump_local = T,
                                    clump_kb = clump_kb,
                                    r2 = r2,
                                    bfile = "../03.reference/1000Gv3EUR/EUR",
                                    method_list =c("mr_ivw","mr_wald_ratio"),
                                    F_method = 2,
                                    Ffilter = 10,
                                    proxies = TRUE,
                                    action = 2,
                                    mr_report = F,
                                    plot=F)
  
  #判断是否存在_final_MRresult.txt文件
  iv_file=paste0(gene,"_",trait,"_final_MRresult.txt")
  if (file.exists(iv_file)) {
  }else{
    warning(gene," 与 ",paste0("../06.pre_data/",trait,".txt"),"进行孟德尔随机化不成功，直接跳过 ",gene," 基因")
    next}

}

