#' does survival analysis
#' 
#' Reads file with features(PC's) in columns and samples in rows. And an annotation file.You pick a max # of features. Using caret
#' runs a prediction iteratively from 2 to max.features. Picks best predictor with least featurs. Returns samples correctly predicted
#' leaves out those samples incorrectly predicted.
#' @param input_folder location of survival files
#' @param output_folder the output folder

#' @importFrom utils read.delim read.table write.table
#'
#' @export


pval_graph<-function(input_folder,output_folder){
  all.files.long.path=list.files(pattern="valorate",full.names=T,path = input_folder)
  all.files.short.path=list.files(pattern="valorate",full.names=F,path = input_folder)
  all.names=sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1])
  for(cancer in unique(all.names)) {
    cancer.logical=str_detect(all.files.long.path,cancer)
    my.files.long=all.files.long.path[cancer.logical]
    my.files.short=all.files.short.path[cancer.logical]
    #temp.vector=vector(length = length(my.files))
    list_of_dataframes<-lapply(my.files.long, read.table,header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
    pval_vector=unlist(lapply(list_of_dataframes,function(x) -1*log(x[,"pval"])*(sign(x[,"Z"])) ))
    threshold_vector=as.numeric(sapply(my.files.short,function(x) strsplit(x,"_")[[1]][4]))
    dat=data.frame(x=threshold_vector,y=pval_vector)
    g=ggplot(dat,aes(x=x,y=y))+geom_line(size=1.2)+geom_point()+
      scale_x_continuous(limits = c(-4,4),breaks=seq(-4,4,by=0.2))+geom_hline(yintercept=2.99,color="red",size=1.5,linetype="dotted")+
      geom_hline(yintercept=-2.99,color="red",size=1.5,linetype="dotted")+geom_text(aes(x=1.5,y=3.3),label="Worse Survival")+
      geom_text(aes(x=1.5,y=-3.25),label="Better Survival")+ theme(axis.text.x= element_text(size=7))
    ggsave(filename = paste0(output_folder,cancer,"_threshold_pval_graph.png"),plot=g,dpi=300,w=8,h=6)
  }
  
  all.files.long.path=list.files(pattern="pancan.summary.threshold_",full.names=T,path = input_folder)
  all.files.short.path=list.files(pattern="pancan.summary.threshold_",full.names=F,path = input_folder)
  all.names=sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1])
  cancer="pancan.summary.threshold"
  for(cancer in unique(all.names)) {
    cancer.logical=str_detect(all.files.long.path,cancer)
    my.files.long=all.files.long.path[cancer.logical]
    my.files.short=all.files.short.path[cancer.logical]
    #temp.vector=vector(length = length(my.files))
    list_of_dataframes<-lapply(my.files.long, read.table,header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F,skip=6)
    pval=as.numeric(strsplit(list_of_dataframes[[1]][1,]," ")[[1]][12])
    pval_vector=unlist(lapply(list_of_dataframes,function(x)  {
      y=strsplit(x[1,]," ")[[1]]
      y=as.numeric(y[y != ""])
      y= -1*sign(y[3]) *log(y[7])
      y
      #-1*sign(as.numeric(strsplit(x[1,]," ")[[1]][4]) )
      # *log(as.numeric(strsplit(x[1,]," ")[[1]][12]))
    }
    ))
    threshold_vector=as.numeric(sapply(my.files.short,function(x) strsplit(x,"_")[[1]][2]))
    dat=data.frame(x=threshold_vector,y=pval_vector)
    g=ggplot(dat,aes(x=x,y=y))+geom_line()+geom_point()+scale_x_continuous(limits = c(-4,4),breaks=seq(-4,4,by=0.2))+geom_hline(yintercept=2.99,color="red",size=1.5,linetype="dotted")+
      geom_hline(yintercept=-2.99,color="red",size=1.5,linetype="dotted")+geom_text(aes(x=1.5,y=3.3),label="Worse Survival")+
      geom_text(aes(x=1.5,y=-3.25),label="Better Survival")+ theme(axis.text.x= element_text(size=7))
    ggsave(filename = paste0(output_folder,cancer,".pval_graph.png"),plot=g,dpi=300,w=8,h=6)
  }
  
  
  all.files.long.path=list.files(pattern="continuous.txt",full.names=T,path = input_folder)
  all.files.short.path=list.files(pattern="continuous.txt",full.names=F,path = input_folder)
  all.names=sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1])
  list_of_dataframes<-lapply(all.files.long.path, read.table,header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F,skip=6)
  
  pval_vector=unlist(lapply(list_of_dataframes,function(x)  {
    y=strsplit(x[1,]," ")[[1]]
    y=as.numeric(y[y != ""])
    y= -1*sign(y[3]) *log10(y[7])
    y
  }))
  names(pval_vector)=all.names
  stuff=as.data.frame(pval_vector)
  stuff=stuff %>% tibble::rownames_to_column()
  colnames(stuff)[1]="sample"
  stuff=stuff[order(stuff[,2],decreasing=T),]
  stuff=rbind(stuff[!stuff$sample %in% "pancan",] , stuff[stuff$sample %in% "pancan",] )
  my.levels=stuff$sample
  stuff1= stuff[stuff$sample == "pancan",]
  stuff2= stuff[stuff$sample != "pancan",]
  
  
  col_vector=c("brown","darkorange1","maroon3","purple4","purple2","#B79F00","deepskyblue","skyblue","springgreen3","magenta1",
               "#9590FF","#556B2F","indianred3","cornflowerblue","limegreen","darkblue","burlywood1","rosybrown3","red3","darkseagreen2","darkviolet","dodgerblue1")
  
  
  
  h=ggplot(data=stuff2, aes(x=factor(sample,levels = my.levels), y=pval_vector, fill=factor(sample,levels = my.levels)))+
    geom_bar(color="black",stat="identity")+ theme_bw()+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),axis.text.x= element_text(size=18),axis.text.y= element_text(size=20),axis.title.y=element_text(size=20),axis.title.x=element_text(size=16),legend.position="none")+
    geom_hline(yintercept=1.3,color="red",size=1.5,linetype="dotted")+ylab("P-value")+
    geom_hline(yintercept=-1.3,color="red",size=1.5,linetype="dotted")+
    guides(fill=guide_legend(title="Cancer Type"))+ggpubr::rotate_x_text(55)+
    xlab("Cancer")+ylim(c(-5,12))+scale_fill_manual(values =  col_vector)
  
  
  

  i=ggplot(data=stuff1, aes(x=factor(sample,levels = my.levels), y=pval_vector, fill=factor(sample,levels = my.levels)))+
    geom_bar(color="black",stat="identity")+ theme_bw()+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),axis.text.x= element_text(size=18),axis.text.y= element_text(size=20),axis.title.y=element_text(size=16),axis.title.x=element_text(size=16),legend.position="none")+
    geom_hline(yintercept=1.3,color="red",size=1.5,linetype="dotted")+ylab("")+
    geom_hline(yintercept=-1.3,color="red",size=1.5,linetype="dotted")+
    guides(fill=guide_legend(title="Cancer Type"))+ggpubr::rotate_x_text(55)+
    xlab("")+ylim(c(-5,12))+scale_color_manual(values =  col_vector[length(col_vector)])
  
  j= ggarrange(plots = list(h,i),widths = c(14/16,2/16),heights=c(8,8),ncol=2,align="h")
  
  ggsave(filename = paste0(output_folder,"pancan.continuous.barchart.pval_graph.png"),plot=j,dpi=300,w=9,h=5)
  
}