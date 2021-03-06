#' create tcga boxplots
#'
#' Creates boxplots from prediction scores from pls or plsda
#' Makes some boxplots
#'
#' @param pattern pattern to grep for in files that have predictions
#' @param tcga_location location of predicted files of tcga
#' @param crpc_file location of predicted file of crpc 
#' @param nepc_file location of predicted file of nepc
#' @param sclc_file location of predicted file of sclc
#' @param output_folder the folder to output to, default is ./ i.e. current folder
#' 
#' @export
#'

do_full_boxplot_RNA=function(file_pattern,tcga_location,crpc_file,nepc_file,sclc_file,output_folder){
  library(data.table)
  all.files.short.path=list.files(tcga_location,pattern=file_pattern,full.names=F)
  all.files.long.path=list.files(tcga_location,pattern=file_pattern,full.names=T)
  all.names=as.character(sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1]))
  list_of_dataframes<-lapply(all.files.long.path, read.table,header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
  names(list_of_dataframes)=all.names
  for(i in 1:length(all.files.short.path)) {
    list_of_dataframes[[i]]$type=all.names[i]
  }
  
  crpc=read.delim(crpc_file)
  crpc.mean=mean(crpc$comp.1)
  crpc.sd=sd(crpc$comp.1)
  crpc.scaled=crpc
  crpc.scaled$comp.1= (crpc.scaled$comp.1-crpc.mean)/crpc.sd
  crpc.scaled$type="CRPC"
  nepc=read.delim(nepc_file)
  colnames(nepc)[1]="Sample"
  nepc$type="NEPC"
  nepc.scaled=nepc
  nepc.scaled$comp.1= (nepc$comp.1-crpc.mean)/crpc.sd
  output.scaled.beltran=rbind(nepc.scaled,crpc.scaled)
  output.scaled.beltran=mutate(output.scaled.beltran, size_n = ifelse(comp.1 >= 3, 5, 3))
  
  lung=rbind(list_of_dataframes[["LUAD"]],list_of_dataframes[["LUSC"]])
  lung$type="NSCLC"
  lung.mean=mean(lung$comp.1)
  lung.sd=sd(lung$comp.1)
  lung.scaled=lung
  lung.scaled$comp.1 = scale(lung.scaled$comp.1)
  george.15=read.delim(sclc_file)
  colnames(george.15)[1]="Sample"
  george.15$type="SCLC"
  george.15.scaled=george.15
  george.15.scaled$comp.1= (george.15$comp.1-lung.mean)/lung.sd
  output.scaled.george=rbind(george.15.scaled,lung.scaled)
  output.scaled.george=mutate(output.scaled.george, size_n = ifelse(comp.1 >= 3, 5, 3))
  
  list_of_dataframes_scaled=lapply(list_of_dataframes,function (x) {
    blah=cbind (x[,1],scale(x[,2:4]),x[,5,drop=F])
    colnames(blah)[[1]]="Sample"
    blah
  })
  output.tcga<-rbindlist(list_of_dataframes)
  output.scaled.tcga<-rbindlist(list_of_dataframes_scaled)
  output.tcga=output.tcga[,c(1,2,5)]
  output.scaled.tcga=output.scaled.tcga[,c(1,2,5)]
  output.scaled.tcga=mutate(output.scaled.tcga, size_n = ifelse(comp.1 >= 3, 5, 3))
  by.top.3=output.scaled.tcga %>% dplyr::group_by(type) %>% dplyr::arrange(type,dplyr::desc(comp.1)) %>% dplyr::slice(1:3) %>% dplyr::summarise(avg=mean(comp.1)) %>% dplyr::arrange(dplyr::desc(avg))
  by.top.3= unique(by.top.3$type)
  bymean=levels(with(output.scaled.tcga, reorder(type, -comp.1, mean)))
  bymax=rev(levels(with(output.scaled.tcga, reorder(type, comp.1, max))))
  col_vector=c("darkseagreen4","indianred1","violetred","tomato2","darkcyan","springgreen3","steelblue1","skyblue","brown","royalblue3","rosybrown3",
               "red3","purple2","darkblue","limegreen","deepskyblue","darkviolet","purple4","burlywood1","magenta1","maroon3","indianred3",
               "cornflowerblue","darkorange1","darkseagreen2","dodgerblue1","firebrick2","cadetblue3","purple2","blue1","chocolate",
               "darkslategray1","aquamarine","palevioletred","chartreuse2","forestgreen","orange3","darkolivegreen2","mediumpurple4","gold")
  set.seed(79)
  library(egg)
  plot1=ggplot(data=output.scaled.beltran,aes(x = factor(type,levels=c("NEPC","CRPC")), y = comp.1))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(shape=19,width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  c(col_vector[1:2],"black"))+ labs(x="",y="Neuroendocrine Score")+scale_x_discrete(labels=c("NEPC","CRPC-Adeno"))+coord_cartesian(ylim = c(-5, 10)) 
  plot2=ggplot(data=output.scaled.george,aes(x = factor(type,levels=c("SCLC","NSCLC")), y = comp.1))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+ 
    geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[3:4])+ labs(x="",y=element_blank()) +coord_cartesian(ylim = c(-5, 10))+
    scale_x_discrete(labels=c("SCLC","NSCLC"))
  plot3=ggplot(data=output.scaled.tcga,aes(x = factor(type,levels=bymax), y = comp.1))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+
    theme_bw()+theme(axis.text=element_text(size=28,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[-c(1:4)])+ labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) 
  
  png(paste0(output_folder,"TCGA_zscore_with_color_bymax.png"),width=2000,height=900)
  ggarrange(plot1,plot2, plot3, widths = c(1/12,1/12, 10/12))
  dev.off()
}
