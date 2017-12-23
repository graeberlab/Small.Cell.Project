#' create tcga boxplots
#'
#' Creates boxplots from prediction scores from pls or plsda
#' Makes some boxplots
#'
#' @param pattern pattern to grep for in files that have predictions
#' @param tcga_location location of predicted files of tcga
#' @param iorio.adeno_file location of predicted file of iorio.adeno 
#' @param iorio.sclc_file location of predicted file of iorio.sclc
#' @param sclc_file location of predicted file of sclc
#' @param output_folder the folder to output to, default is ./ i.e. current folder
#' 
#' @export
#'

do_full_boxplot_meth=function(file_pattern,tcga_location,iorio.adeno_file,iorio.sclc_file,sclc_file,output_folder){
  library(data.table)
  all.files.short.path=list.files(tcga_location,pattern=file_pattern,full.names=F)
  all.files.long.path=list.files(tcga_location,pattern=file_pattern,full.names=T)
  all.names=as.character(sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1]))
  list_of_dataframes<-lapply(all.files.long.path, read.table,header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
  names(list_of_dataframes)=all.names
  for(i in 1:length(all.files.short.path)) {
    list_of_dataframes[[i]]$type=all.names[i]
  }
  
  iorio.adeno=read.delim(iorio.adeno_file)
  iorio.adeno.mean=mean(iorio.adeno$comp.1)
  iorio.adeno.sd=sd(iorio.adeno$comp.1)
  iorio.adeno.scaled=iorio.adeno
  iorio.adeno.scaled$comp.1= (iorio.adeno.scaled$comp.1-iorio.adeno.mean)/iorio.adeno.sd
  iorio.adeno.scaled$type="lung_adeno_line"
  iorio.sclc=read.delim(iorio.sclc_file)
  colnames(iorio.sclc)[1]="Sample"
  iorio.sclc$type="lung_small_line"
  iorio.sclc.scaled=iorio.sclc
  iorio.sclc.scaled$comp.1= (iorio.sclc$comp.1-iorio.adeno.mean)/iorio.adeno.sd
  output.scaled.iorio=rbind(iorio.sclc.scaled,iorio.adeno.scaled)
  output.scaled.iorio=mutate(output.scaled.iorio, size_n = ifelse(comp.1 >= 3, 5, 3))
  
  luad=rbind(list_of_dataframes[["LUAD"]])
  luad$type="LUAD"
  luad.mean=mean(luad$comp.1)
  luad.sd=sd(luad$comp.1)
  luad.scaled=luad
  luad.scaled$comp.1 = scale(luad.scaled$comp.1)
  sclc.tumor=read.delim(sclc_file)
  colnames(sclc.tumor)[1]="Sample"
  sclc.tumor$type="SCLC"
  sclc.tumor.scaled=sclc.tumor
  sclc.tumor.scaled$comp.1= (sclc.tumor$comp.1-luad.mean)/luad.sd
  output.scaled.sclc.tumor=rbind(sclc.tumor.scaled,luad.scaled)
  output.scaled.sclc.tumor=mutate(output.scaled.sclc.tumor, size_n = ifelse(comp.1 >= 3, 5, 3))
  
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
  plot1=ggplot(data=output.scaled.iorio,aes(x = factor(type,levels=c("lung_small_line","lung_adeno_line")), y = comp.1))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(shape=19,width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  c(col_vector[1:2],"black"))+ labs(x="",y="Neuroendocrine Score")+scale_x_discrete(labels=c("lung_small_line","lung_adeno_line"))+coord_cartesian(ylim = c(-5, 10)) 
  plot2=ggplot(data=output.scaled.sclc.tumor,aes(x = factor(type,levels=c("SCLC","NSCLC")), y = comp.1))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+ 
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
