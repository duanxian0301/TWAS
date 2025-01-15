#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/08.UTMOST_ensTOsymble")

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

#搜寻UTMOST中的跨组织结果
cross_tissue_result=list.files("../07.UTMOST/result/",
                     pattern = paste0(trait,"_GTEXv8_49_cross"), 
                     full.names = T)

#创建result文件夹
if (!file.exists("result")) {
  dir.create("result")
}


#指定GTExv8基因转换参考文件
GTExv8_ref="../03.reference/GTEX_V8_gene_info.txt"
#读入参考文件
GTEX_V8_gene_id=data.table::fread(GTExv8_ref)

#合并all_single_result结果------------------------------------------------------
data_single=data.table::fread(paste0("../07.UTMOST/result/",trait,"_GTExv8.UTMOST_all.txt"))
data_single$gene=substring(data_single$gene, 1, 15)

#合并
data_single_merge=data.frame(merge(data_single,GTEX_V8_gene_id,by.x = "gene", by.y = "gene_id",all.x = T)) 

#保存
data.table::fwrite(data_single_merge,
                   file =paste0("./result/",trait,"_GTExv8.UTMOST_all_symble.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

#分组织保存
TISSUE_GTEx=FUSION_tissue()[c(6:54),"Tissue"] %>% gsub("GTExv8.ALL.","",.)

for (tissue in TISSUE_GTEx) {
  #筛选
  data_single_tissue=data_single_merge[data_single_merge$tissue==tissue,]
  
  #保存
  data.table::fwrite(data_single_tissue,
                     file =paste0("./result/",trait,"_GTExv8_",tissue,"_symble.txt"),
                     sep = "\t",
                     na = "NA",
                     quote = FALSE)
}


#合并cross_tissue_result结果----------------------------------------------------
data_cross=data.table::fread(cross_tissue_result)
data_cross$gene=substring(data_cross$gene, 1, 15)

#合并
data_cross_merge=merge(data_cross,GTEX_V8_gene_id,by.x = "gene", by.y = "gene_id",all.x = T) 

#保存
data.table::fwrite(data_cross_merge,
                   file =paste0("./result/",trait,"_GTEXv8_49_cross_symble.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)


