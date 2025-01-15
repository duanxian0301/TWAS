
#设置当前工作路径
setwd("D:/1.全基因转录组/5.TWAS_UTMOST/22.Colocalization")

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
tissue=FUSION_tissue()[c(44),"Tissue"]%>% gsub("GTExv8.ALL.","",.) # 选取GTEx V8中的Whole_Blood组织
P_method = "fdr" #可选择矫正方式
trait="Prostate cancer"
wind_size=500000
type1="quant"         #eqtl的性状类型，则type1="quant"
type2="cc"            #GWAS的性状类型：cc/quant；若是GWAS性状是case-control，则type2="cc"；若是GWAS性状是连续性变量，则type2="quant"
case_pro1=NULL        #type1参数为"quant",则case_pro1可以输入NULL
case_pro2=122188 /(726828)   #type2参数为"cc",则必须添加case_pro2，输入GWAS2数据集中病例与样本量的比值
width = 7
height = 7
#133,384 women affected by breast cancer and 113,789 control 
#读入基因
intersect_gene=data.table::fread(paste0("../19.Venn_UTMOST_FUSION_MAGMA_FOCUS/",trait,"_Venn_UTMOST_FUSION_MAGMA_FOCUS.txt"))


# eqtl数据位置
path_hg19 <- "../21.pre_eqtlsum_hg38Tohg19/pre_coloc_hg19/"
if (!dir.exists("./result/")) {
  dir.create("./result/")
}

#读入数据
FUSION_file=paste0("../12.sig_FUSION/",trait,".GTExv8.ALL.",tissue,".fusion_twas_symble.",P_method,".sig.txt")
data_fusion=data.table::fread(FUSION_file)

#循环计算共定位结果并绘制图
for (i in 1:length(intersect_gene$gene_id)) {
  gene=intersect_gene$gene_id[i]
  message("当前第",i,"个基因: ",gene)
  print(intersect_gene[i,])
  coloc_local_gwas2gwas(gwasfile=c(paste0(path_hg19,intersect_gene$Gene[i],"_",gene,"_eqtlsum_from_GTEx_",tissue,"_hg19.txt"),
                                   paste0("../06.pre_data/",trait,".txt")),
                        wind_size = wind_size,
                        gene = gene,
                        SNP=unique(data_fusion[data_fusion$ensembl==gene,]$BEST.GWAS.ID),
                        bfile = "../03.reference/1000Gv3EUR/EUR",
                        plink_bin = friendlypostGWASsoft::plink_binary(),
                        outfile = paste0("./result/",intersect_gene$Gene[i],"_",wind_size/1000,"kb.local"),
                        width = width,
                        height = height,
                        build = 37,
                        EUR = "EUR",
                        plot_locuscompare = TRUE,
                        plot_locuscompare_title1="eQTL",
                        plot_locuscompare_title2="GWAS",
                        plot_gassocplot2 = TRUE,
                        type1 = type1,
                        type2 = type2,
                        trait_name1 = intersect_gene$Gene[i],
                        trait_name2 = trait,
                        case_pro1 = case_pro1,
                        case_pro2 = case_pro2)
  
}

#合并共定位结果-----------------------------------------------------------------
all_coloc=data.frame()
for (i in 1:length(intersect_gene$gene_id)) {
  gene=intersect_gene$gene_id[i]
  message("当前第",i,"个基因: ",gene)
  print(intersect_gene[i,])
  
  #入读文件
  data=data.table::fread(paste0("./result/",intersect_gene$Gene[i],"_",wind_size/1000,"kb.local.coloc_summary.txt"))
  data$Gene=intersect_gene$Gene[i]
  all_coloc=rbind(all_coloc,data)
  
}

# 合并
data_coloc_merge=merge(all_coloc,intersect_gene,by = "Gene")
data_coloc_merge=data_coloc_merge[order(data_coloc_merge$PP.H4.abf,decreasing = T),]

#保存
data.table::fwrite(data_coloc_merge,
                   file =paste0(trait,".",tissue,"_Venn_UTMOST_FUSION_MAGMA_FOCUS_coloc.txt"),
                   sep = "\t",
                   na = "NA",
                   quote = FALSE)



