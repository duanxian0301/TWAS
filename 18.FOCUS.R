

#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/18.FOCUS")

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
pip_th=0.8  #设置pip阈值
# 选取GTEx V8中的Whole_Blood组织
tissue=FUSION_tissue()[c(44),"Tissue"]

focusdata=TWAS_focus_format_clean(gwasfile=paste0("../06.pre_data/",trait,".txt"),
                                  flag_options="--chunksize 500000",
                                  outname = trait,
                                  force = F) 

#注意！在windows系统中TWAS_focus_fine_mapping_multiprocess函数中的参数路径中不能出现冒号(:)，
#即是说【不能写全部路径】，不然会出错！
#locations 指定BED文件路径，可指定'37:EUR', '37:AFR', '37:EAS', '38:EUR', '38:AFR', '38:EAS'，但这样部分染色体会报错
################################################################################
####这一步消耗较多时间（根据电脑性能的不同，时长可达到数个小时）################
################################################################################
focusall=TWAS_focus_fine_mapping_multiprocess(gwasfile=paste0("./result/",trait,".focus.sumstats.gz"),
                                              eQTL_db_weight = paste0("../03.reference/FOCUS/focus_GTEx_V8_from_fusion_",tissue,".db"),
                                              LDREF = "../03.reference/FOCUS/1000G_EUR_Phase3_plink/1000G.EUR.QC",
                                              outname = paste0(trait,"_",tissue,"_focus"),
                                              startChr = 1,
                                              endChr = 22,
                                              plot = TRUE,
                                              flag_options=NULL,
                                              locations  = "../03.reference/FOCUS/grch37.eur.loci - 副本.bed",
                                              pval_threshold = 5e-8,
                                              nthread = 8,
                                              force = F)


#筛选位于90% 置信区间且pip>0.8的基因
focusall_sig=focusall[focusall$in_cred_set_pop1==1 & focusall$pips_pop1 > pip_th & focusall$ens_gene_id!="NULL.MODEL",]


#保存
data.table::fwrite(focusall_sig,
                   file =paste0(trait,"_",tissue,"_focus.sig.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

