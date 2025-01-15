
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/10.FUSION")

#加载包
library(Fast2TWAS)

#输入参数
trait="Prostate cancer"
# 选取GTEx V8中的Whole_Blood组织,可调整数值，改换其他组织，但是要另行下载数据
tissue=FUSION_tissue()[c(44),"Tissue"] 

#处理GWAS数据生成LDSC格式的数据用于FUSION软件的输入文件-------------------------
ldsc=format_LDSC_data(file = paste0("../06.pre_data/",trait,".txt"),
                      SNP="SNP",
                      BETA="beta",
                      SE="se",
                      effect_allele="effect_allele",
                      other_allele = "other_allele",
                      pval = "pval",
                      samplesize = "N",
                      savefile=trait,
                      sep = "\t") 


#FUSION软件执行TWAS分析，获取基因与疾病之间的TWAS分析结果，可以多线程进行-------
#生成sh脚本调用FUSION中的R脚本进行TWAS/PWAS分析(window，linux，MAC均可运行)
allresult=TWAS_fusion_multiprocess(ref_weights_dir = "../03.reference/FUSION/",
                                   LDREF = "../03.reference/FUSION/LDREF/1000G.EUR",
                                   sumstatsfile = paste0(trait,".sumstats.gz"),
                                   outname = trait,
                                   tissue = tissue,
                                   startChr = 1,
                                   endChr = 22,
                                   nthread = 4)


