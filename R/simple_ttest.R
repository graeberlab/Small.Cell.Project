#' does survival analysis
#' 
#' Reads file with features(PC's) in columns and samples in rows. And an annotation file.You pick a max # of features. Using caret
#' runs a prediction iteratively from 2 to max.features. Picks best predictor with least featurs. Returns samples correctly predicted
#' leaves out those samples incorrectly predicted.
#' 
#' @param dataset  data file, first column is features , columns samples, rows genes
#' @param fact  first column sample names, second column type, must have no header
#' @param file.prefix name of output file will be postended with ".txt"
#' @importFrom genefilter rowttests
#'
#' @export

simple_ttest<-function(dataset,fact,file.prefix){
  dat <- read.delim(dataset,stringsAsFactors=F,header=T,row.names=1)
  factr<-read.delim(fact,header=F)
  common=intersect(factr[,1],colnames(dat))
  factr=factr[factr[,1] %in% common,]
  dat=dat[,colnames(dat) %in% common]
  types=factr[match(colnames(dat),factr[,1]),2]
  colnames(dat)<-types
  out<- rowttests(as.matrix(dat),fac=as.factor(types))
  out<- na.omit(out)
  out$logpval<- -log10(out$p.value) *sign(out$statistic)
  out<-out[order(out$logpval,decreasing=T),]
  
  out<-out[,-c(1,2,3),drop=F]
  write.table(out,paste0(file.prefix,".txt"),col.names = F,quote=F,sep="\t")
}

