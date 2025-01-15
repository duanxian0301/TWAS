

#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/14.MAGMA")

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


# 调用MAGMA软件实现gene-based & gene-set-based关联分析 -------------------------
################################################################################
####慎重！慎重！慎重！magma_gene_based()，######################################
####这一步消耗较多时间（根据电脑性能的不同，时长可达到数个小时）################
################################################################################
magma_gene_based(bfile="../03.reference/magma/g1000_eur/g1000_eur",
                 gwasfile = paste0("../06.pre_data/",trait,".txt"),
                 gene_loc = "../03.reference/magma/NCBI37/ENSGv110.coding.genes.txt",
                 set_annot = "../03.reference/magma/MSigDB_20231Hs_MAGMA.txt",
                 SNP_P_col = c(3,9),
                 outname = paste0(trait,"_magma"),
                 flag_1_options = NULL,
                 flag_2_options = NULL,
                 flag_3_options = NULL,
                 N_col = "N") 

#调用MAGMA软件进行Gene property analysis for tissue specificity分析-------------
magma_gene_property(raw_file=paste0("./result/",paste0(trait,"_magma"),".genes.raw"),
                    gene_covar = "../03.reference/magma/gtex_v8_ts_avg_log2TPM.txt",
                    flag_options = NULL,
                    outname = paste0(trait,"_tissue_specific_magma_result"))


