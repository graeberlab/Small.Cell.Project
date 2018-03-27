#' PCA from file
#' 
#' Reads file with samples in columns and variables in rows, and does PCA. Writes to file scores, loadings, eigenvalues.
#' 
#' @param file Filepath/filename of data matrix
#' @param center default=T
#' @param scale default=F
#' @param rank. how many pcs to calculate, default is equal to # of samples, useful when you only really need a few pcs
#' 
#' @importFrom stats prcomp screeplot
#' @importFrom utils read.delim read.table write.table
#' 
#' @export
#' 


PCA_from_file=function(file,center=TRUE,scale=FALSE,rank.=NULL,rotate=F) {
  data=read.delim(file, header = T, stringsAsFactors = F)
  
  #remove rows that are all 0
  data= data[rowSums((data[,-1]==0))<ncol(data[-1]),]
  
  #t.data = t(data) #if genenames in rownames
  t.data=t(data[,-1]) ##subtract off the gene name
  pca<-prcomp(t.data,scale=scale,center=center,rank.=rank.);
  pca_scores=pca$x
  pca_scores=cbind("Score"=rownames(pca_scores),pca_scores)
  pca_loadings=pca$rotation
  #pca_loadings=cbind("Loading"=rownames(pca_loadings),pca_loadings)#if genenames in rownames
  pca_loadings=cbind("Loading"=data[,1],pca_loadings)#if genenames not in rownames
  pca_evalues=pca$sdev
  if(rotate==T){
   pca_scores[,"PC1"]=-1*as.numeric(pca_scores[,"PC1"])
   pca_loadings[,"PC1"]=-1*as.numeric(pca_loadings[,"PC1"])
    
  }
  
  #save data
  name=sub(".txt","",file)
  savename=paste(name,"_prcomp_scores.txt",sep='');
  write.table(pca_scores,savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_prcomp_loadings.txt",sep='');
  write.table(pca_loadings,savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_prcomp_sdev.txt",sep='');
  write.table(pca_evalues,savename,sep='\t',row.names=FALSE,quote=FALSE);
  print(summary(pca))
  screeplot(pca)
  
}