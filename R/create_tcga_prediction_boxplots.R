#' Partial Least Squares and predict second dataset
#'
#' Builds PLS model from training dataset and predicts second dataset
#' Writes out predicted.scores of the second dataset
#' Plots original PLS, projected samples only, and projected samples ontop of original PLS
#'
#' @param file file for X matirx
#' @param sample.names Vector of sample names in X matrix
#' @param sample.type vector of sample groups
#' @param y.response numeric vector of response values in same order as samples
#' @param comp number of components to compute
#' @param scale default=T
#' @param labels label the plot, default = T
#' @param comp.x,comp.y comps to display
#' @param title title of the plot
#' 
#' @param file2 file for test data matrix
#' @param sample.names2 Vector of sample names in 2nd dataset, if needed
#' @param sample.type2 Vector of sample types in 2nd dataset, if needed
#' @param train_string string of training data to insert in file name of predicted scores
#' @param test_string string of data being tested to insert in file name of predicted scores
#' @param saveplot whether to save the plot, default is F
#' @param savetype the type of plot to save,options are ".pdf" or ".png"
#' @param w is width of plot to be saved
#' @param h is height of plot to be saved
#' @param legendname is the legend name
#' @param plot_both if true plots both training and test set in color
#' @param colpalette allows you to put in a color palette of form c("#F87660", "#39B600",....etc) to manually assign colors
#' @param shape.palette allows you to put in a shape palette of form c(1, 3,....etc) to manually assign shapes
#' @param varimax If T performs Varimax rotation, 
#' @param varimax.comp # of varimax components, kind of hacky, keep this # the same as # of comps. Will fix later.
#' @param output_folder the folder to output to, default is ./ i.e. current folder
#' @param TCGA predicted files are from TCGA, barcodes separated by periods, so remove normal samples, default is FALSE
#'
#' @importFrom mixOmics pls
#' @export
#'

tcga_boxplots=function(dir,pattern,output_string,output_folder="./"){
  col_vector=c("darkseagreen4","indianred1","violetred","tomato2","darkcyan","springgreen3","steelblue1","skyblue","brown","royalblue3","rosybrown3",
               "red3","purple2","darkblue","limegreen","deepskyblue","darkviolet","purple4","burlywood1","magenta1","maroon3","indianred3",
               "cornflowerblue","darkorange1","darkseagreen2","dodgerblue1","firebrick2","cadetblue3","purple2","blue1","chocolate",
               "darkslategray1","aquamarine","palevioletred","chartreuse2","forestgreen","orange3","darkolivegreen2","mediumpurple4","gold")
  all.files.short.path=list.files(dir,pattern=pattern,full.names=F)
  all.files.long.path=list.files(dir,pattern=pattern,full.names=T)
  all.names=as.character(sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1]))
  list_of_dataframes<-lapply(all.files.long.path, read.table,header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
  names(list_of_dataframes)=all.names
  for(i in 1:length(all.files.short.path)) {
    list_of_dataframes[[i]]$type=all.names[i]
  }
  list_of_dataframes_scaled=lapply(list_of_dataframes,function (x) {
    cols=ncol(x)
    blah=cbind (x[,1],scale(x[,2:(cols-1)]),x[,cols,drop=F])
    colnames(blah)[[1]]="Sample"
    blah
  })
  output<-data.frame(rbindlist(list_of_dataframes))
  output.scaled<-data.frame(rbindlist(list_of_dataframes_scaled))
  
  
  output.scaled.by.top3.avg=output.scaled %>% dplyr::group_by(type) %>% dplyr::arrange(type,dplyr::desc(comp.1)) %>% dplyr::slice(1:3) %>% dplyr::summarise(avg=mean(comp.1)) %>% dplyr::arrange(dplyr::desc(avg))
  scaled.by.top.3= unique(output.scaled.by.top3.avg$type)
  scaled.bymean=levels(with(output.scaled, reorder(type, -comp.1, mean)))
  scaled.bymax=rev(levels(with(output.scaled, reorder(type, comp.1, max))))
  
  unscaled.bymean=levels(with(output, reorder(type, -comp.1, mean)))
  normalize01 <- function(x){
    return((x-min(x)) / (max(x)-min(x)))
  }
  ## TCGA data unscaled 
  set.seed(79)
  unscaled.pic<-ggplot(data=output,aes(x = factor(type,levels=unscaled.bymean), y = comp.1))+geom_boxplot(outlier.color=NA)+
    geom_jitter(width=0.1,aes(color=factor(type)))+theme(axis.text=element_text(size=10,face="bold",angle=90),axis.title=element_text(size=8,face="bold"),legend.position="none")+
    scale_color_manual(values =  sample(col_vector,size=40))
  ggsave(paste0(output_folder,output_string,"_TCGA_unscaled_by_mean.png"), 
         dpi = 300, plot = unscaled.pic, width = 8, height = 5)
  # TCGA scaled by top 3
  set.seed(79)
  scaled.pic<-ggplot(data=output.scaled,aes(x = factor(type,levels=scaled.by.top.3), y = comp.1))+geom_boxplot(outlier.color=NA)+  
    geom_jitter(width=0.1,aes(color=factor(type)))+theme(axis.text=element_text(size=10,face="bold",angle=90),axis.title=element_text(size=8,face="bold"),legend.position="none")+
    scale_color_manual(values =  sample(col_vector,size=40))
  ggsave(paste0(output_folder,output_string,"_TCGA_scaled_by_avg_top3.png"), 
         dpi = 300, plot = scaled.pic, width = 8, height = 5)
}