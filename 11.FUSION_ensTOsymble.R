
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/11.FUSION_ensTOsymble/")

#加载包
library(Fast2TWAS)

#输入参数
trait="Prostate cancer"
# 选取GTEx V8中的Whole_Blood组织,可调整数值，改换其他组织，但是要另行下载数据
tissue=FUSION_tissue()[c(44),"Tissue"] 


#读入数据
inputfile=paste0("../10.FUSION/result/Gene_disease_association/",trait,".",tissue,".fusion_twas.txt")
data=data.table::fread(inputfile)

#指定GTExv8基因转换参考文件
GTExv8_ref="../03.reference/GTEX_V8_gene_info.txt"
#读入参考文件
GTEX_V8_gene_id=data.table::fread(GTExv8_ref)


#合并all_single_result结果------------------------------------------------------
data$ID=substring(data$ID, 1, 15)

#合并
data_single_merge=data.frame(merge(data,GTEX_V8_gene_id,by.x = "ID", by.y = "gene_id",all.x = T)) 
data_single_merge$ensembl=data_single_merge$ID
data_single_merge$ID=data_single_merge$symble
# 去除指定列
data_single_merge <- subset(data_single_merge, select = -c(symble,tissue))

data.table::fwrite(data_single_merge,
                   file =paste0(trait,".",tissue,".fusion_twas_symble.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)

