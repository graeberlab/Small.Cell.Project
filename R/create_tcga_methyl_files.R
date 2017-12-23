#' Create methylation files for TCGA cancers using defined sets of sites
#'
#' Takes illumina methylation files and creates new individualfiles using features that exist in the mysites files 
#'
#' @param file input methylation file , file name is of form "ACC_etc_etc...."
#' @param subset Whether to restrict to sites in mysites file,default T
#' @param mysites vector of site names
#' @param out.string If subset=T Should be the string explaining the sites used 
#' @param normals If T only returns normals, default is false so returns only cancer samples 
#'
#' @export
#'
create_tcga_methyl_files=function(file,subset=TRUE,mysites=NULL,out.string="blah",normals=F,output_folder="./"){
  methyl.pre=fread(file,sep="\t",header =T, stringsAsFactors = F,check.names=F,nrows=0)
  methyl.colnums=c(1,seq(2,ncol(methyl.pre),4))
  methyl=fread(file,sep="\t",header =T, stringsAsFactors = F,check.names=F,select=methyl.colnums)
  colnames(methyl)[-1]=substr(colnames(methyl)[-1],1,15)
  colnames(methyl)[-1]=make.names(colnames(methyl)[-1])
  methyl=methyl[-1,]
  my_samps.normals=c(1,which(as.numeric(sapply(colnames(methyl),function(x) strsplit(x,"\\.")[[1]][4])) > 9))
  my_samps.cancers=c(1,which(as.numeric(sapply(colnames(methyl),function(x) strsplit(x,"\\.")[[1]][4])) <= 9))
  methyl=as.data.frame(methyl)
  if(length(my_samps.normals)==0){
    methyl.cancers=methyl[ ,my_samps.cancers,drop=F]
    colnames(methyl.cancers)[1]="Site"
    methyl.cancers=na.omit(methyl.cancers)
    
    if(subset==T){
      methyl.cancers=methyl.cancers[methyl.cancers[,1] %in% mysites,]
    } 
    nam=strsplit(strsplit(file,"/")[[1]][length(strsplit(file,"/")[[1]])],"\\.")[[1]][1]
    write.table(methyl.cancers,paste0(output_folder,nam,"_cancers_methyl_","450K_",out.string,".txt"),sep="\t",quote=F,row.names=F)
    
  } else{
    methyl.normals=methyl[ ,my_samps.normals,drop=F]
    methyl.cancers=methyl[ ,my_samps.cancers,drop=F]
    colnames(methyl.normals)[1]="Site"
    colnames(methyl.cancers)[1]="Site"
    methyl.normals=na.omit(methyl.normals)
    methyl.cancers=na.omit(methyl.cancers)
    if(subset==T){
      methyl.cancers=methyl.cancers[methyl.cancers[,1] %in% mysites,]
      methyl.normals=methyl.normals[methyl.normals[,1] %in% mysites,]
    }
    
    write.table(methyl.normals,paste0(output_folder,nam,"_normals_methyl_","450K_",out.string,".txt"),sep="\t",quote=F,row.names=F)
    write.table(methyl.cancers,paste0(output_folder,nam,"_cancers_methyl_","450K_",out.string,".txt"),sep="\t",quote=F,row.names=F)
  }
}



