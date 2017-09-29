create_log_tcga_files=function() {
  log2_it=function(data){
    gene=data[,1]
    dat=log2(data[,-1]+ 1)
    dat= cbind(gene,dat)
    return(dat)
  }
  
  coding = read.delim("protein-coding_gene.txt")
  all.files.short.path=list.files("//10.47.223.100/data2/users/nbalanis/Input_Files_and_Standards/TCGA_TOIL_RECOMPUTE/TOIL_TCGA_norm",pattern="_norm.txt",full.names=F)
  all.names=as.character(sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][3]))
  for(name in all.names){
    data=read.delim(paste0("//10.47.223.100/data2/users/nbalanis/Input_Files_and_Standards/TCGA_TOIL_RECOMPUTE/TOIL_TCGA_norm/TOIL_TCGA_",name,"_norm.txt"), header =T, stringsAsFactors = F,check.names=F)
    cancer_samples=which(as.numeric(sapply(colnames(data),function(x) strsplit(x,"\\.")[[1]][4])) <= 9) 
    data=cbind(data[,1],data[,cancer_samples]) 
    colnames(data)[1]="gene"  
    data=data[data$gene %in% coding$symbol,]
    data=data[order(data$gene),]
    data=log2_it(data)
    write.table(data,paste0(name,"_rsem_genes_upper_norm_counts_coding_log2.txt"),sep="\t",quote=F,row.names = F)
  }
}