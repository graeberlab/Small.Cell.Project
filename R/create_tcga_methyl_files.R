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
create_tcga_methyl_files=function(file,subset=TRUE,mysites=NULL,out.string="blah",normals=F){
  methyl.pre=fread(file,sep="\t",header =T, stringsAsFactors = F,check.names=F,nrows=0)
  colnames(methyl)[-1]=substr(colnames(methyl)[-1],1,15)
  colnames(methyl)[-1]=make.names(colnames(methyl)[-1])
  methyl.colnums=c(1,seq(2,ncol(methyl.pre),4))
  methyl=fread(file,sep="\t",header =T, stringsAsFactors = F,check.names=F,select=methyl.colnums)
  methyl=methyl[-1,]
  if(normals==T) {
    my_samps=which(as.numeric(sapply(colnames(methyl.pre),function(x) strsplit(x,"\\.")[[1]][4])) > 9)
  } else {
  my_samps=c(1,which(as.numeric(sapply(colnames(methyl),function(x) strsplit(x,"\\.")[[1]][4])) <= 9))
  }
  methyl=as.data.frame(methyl)
  methyl=methyl[ ,my_samps]
  if(subset==T){
    methyl=methyl[methyl[,1] %in% mysites,]
  } 
  colnames(methyl)[1]="Site"
  methyl=na.omit(methyl)
  nam=strsplit(strsplit(file,"/")[[1]][length(strsplit(file,"/")[[1]])],"\\.")[[1]][1]
  write.table(methyl,paste0(nam,"_methyl_","450K_",out.string,".txt"),sep="\t",quote=F,row.names=F)
  
}


