

#设置当前工作路径
setwd("D:/5.TWAS_UTMOST/07.UTMOST")

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


#UTMOST软件执行TWAS分析，获取GTExV8中单个组织中基因与疾病之间的TWAS分析结果-----
################################################################################
###########################该步骤很耗时间请耐心等待#############################
################################################################################
UTMOST_single=TWAS_UTMOST_single_tissue(gwasfile='../06.pre_data/Prostate cancer.txt',
                                        db_path_dir= "../03.reference/UTMOST/database_normalized_pruned/",
                                        covariance_path_dir   = "../03.reference/UTMOST/covariance_GTEx8_normalized_pruned/",
                                        outname = paste0(trait,"_GTExv8"),
                                        tissue="GTExv8",
                                        SNP = "SNP",
                                        effect_allele = "effect_allele",
                                        other_allele = "other_allele",
                                        BETA = "beta",
                                        pval = "pval",
                                        flag_options=NULL,
                                        force = F,
                                        nthread = 8)


#UTMOST软件执行TWAS分析,计算GTExV8中跨组织联合测试关联结果----------------------
################################################################################
###########windows系统中该程序可能多次运行都不成功，请耐心多试几次##############
################################################################################
UTMOST_cross=TWAS_UTMOST_cross_tissue(db_path_dir =paste0(dirname(getwd()),"/03.reference/UTMOST/database_normalized_pruned/"),
                                      covariance_joint_path_dir  = paste0(dirname(getwd()),"/03.reference/UTMOST/covariance_joint_GTEx8_normalized_pruned"),
                                      gene_info_file = "./result/gene_info_GTExv8.txt",
                                      flag_options=NULL,
                                      outname  = paste0(trait,"_GTEXv8_49_cross"),
                                      force = F)


