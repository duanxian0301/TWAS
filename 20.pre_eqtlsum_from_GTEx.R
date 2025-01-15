
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/20.pre_eqtlsum_from_GTEx")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#加载包
if (!requireNamespace("memuse", quietly = TRUE))install.packages("memuse")
library(Fast2TWAS)

#输入参数
trait="Prostate cancer"
# 选取GTEx V8中的Whole_Blood组织
tissue=FUSION_tissue()[c(44),"Tissue"]%>% gsub("GTExv8.ALL.","",.)


#读入基因
sig_gene_GTExv8_ref=data.table::fread(paste0("../19.Venn_UTMOST_FUSION_MAGMA_FOCUS/",trait,"_Venn_UTMOST_FUSION_MAGMA_FOCUS.txt"))

# 保存位置
path <- "./pre_coloc/"  
if (!dir.exists(path)) {
  dir.create(path)
}


# R语言调用shell脚本从GTExV8中的全血中获取指定基因共定位的eQTL summary data数据-
for (i in 1:length(sig_gene_GTExv8_ref$gene_id)) {
  gene=sig_gene_GTExv8_ref$gene_id[i]
  message("当前第",i,"个基因: ",gene)
  print(sig_gene_GTExv8_ref[i,])
  sh_file="../02.software/get_eqlsum_form_GTExV8_local.sh"
  tissue_file=paste0("../03.reference/GTExv8/",tissue,".tsv.gz")
  cmd=paste("sh",sh_file,gene,tissue_file,path)
  status=system(cmd)
}



# R语言从GTExV8中的全血中获取指定基因共定位的eQTL summary data数据--------------
################################################################################
################GTExV8本地全血eQTL文件太大，需要电脑有较大内存##################
################################################################################
library(memuse)
Freeram=as.numeric(gsub(" GiB","",as.character(Sys.meminfo()[["freeram"]])))

if(Freeram>30){
#读入数据
GTEx_data=data.table::fread(paste0("../03.reference/GTExv8/",tissue,".tsv.gz"))

for (i in 1:length(sig_gene_GTExv8_ref$gene_id)) {
  gene=sig_gene_GTExv8_ref$gene_id[i]
  message("当前第",i,"个基因: ",gene)
  print(sig_gene_GTExv8_ref[i,])
  GTEx_data_gene_id=GTEx_data[GTEx_data$gene_id==gene,]
  GTEx_data_gene_id=subset(GTEx_data,gene_id==gene)
  
  #添加样本量列
  GTEx_data_gene_id$N=GTEx_data_gene_id$an/2
  
  #转换为Fast2TWAS::show_demo_gwas()格式
  pre_data=GTEx_data_gene_id[,c("chromosome","position","rsid","alt","ref","maf","beta","se","pvalue","N")]
  colnames(pre_data)=colnames(show_demo_gwas())
  
  #保存
  data.table::fwrite(data.frame(pre_data),
                     file =paste0(gene,"_eqtlsum_from_GTEx_Whole_Blood_hg38.txt"),
                     sep = "\t",
                     na = "NA",
                     quote = FALSE)
  }  
}
