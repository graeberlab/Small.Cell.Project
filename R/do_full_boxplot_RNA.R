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

do_full_boxplot_RNA=function(file_pattern,tcga_location,crpc_file,nepc_file,sclc_file,output_folder,comp.x="comp.1",components=3,threshold=3){
  library(data.table)
  
  epithelial=c('ACC','BLCA','BRCA','CESC','CHOL','COADREAD','ESCA','HNSC','KICH','KIRC','KIRP','LIHC','LUAD','LUSC','OV','PAAD','PRAD',
               'STAD','THCA','UCEC','UCS')

  
  all.files.short.path=list.files(tcga_location,pattern=file_pattern,full.names=F)
  all.files.long.path=list.files(tcga_location,pattern=file_pattern,full.names=T)
  all.names=as.character(sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1]))
  list_of_dataframes<-lapply(all.files.long.path, read.table,header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
  names(list_of_dataframes)=all.names
  for(i in 1:length(all.files.short.path)) {
    list_of_dataframes[[i]]$type=all.names[i]
  }
  
  crpc=read.delim(crpc_file)
  crpc$type="CRPC"
  crpc.mean=mean(crpc[,comp.x])
  crpc.sd=sd(crpc[,comp.x])
  crpc.scaled=crpc
  crpc.scaled[,comp.x]= (crpc.scaled[,comp.x]-crpc.mean)/crpc.sd
  crpc.scaled$type="CRPC"
 
  nepc=read.delim(nepc_file)
  colnames(nepc)[1]="Sample"
  nepc$type="NEPC"
  nepc.scaled=nepc
  nepc.scaled[,comp.x]= (nepc[,comp.x]-crpc.mean)/crpc.sd
  output.scaled.beltran=rbind(nepc.scaled,crpc.scaled)
  output.beltran=rbind(nepc,crpc)
  output.scaled.beltran$size_n=ifelse(output.scaled.beltran[,comp.x] >=threshold,5,3 )
 # output.beltran$size_n=ifelse(output.beltran[,comp.x] >=threshold,5,3 )
  
  lung=rbind(list_of_dataframes[["LUAD"]],list_of_dataframes[["LUSC"]])
  lung$type="NSCLC"
  lung.mean=mean(lung[,comp.x])
  lung.sd=sd(lung[,comp.x])
  lung.scaled=lung
  lung.scaled[,comp.x] = scale(lung.scaled[,comp.x])
  george.15=read.delim(sclc_file)
  colnames(george.15)[1]="Sample"
  george.15$type="SCLC"
  george.15.scaled=george.15
  george.15.scaled[,comp.x]= (george.15[,comp.x]-lung.mean)/lung.sd
  output.scaled.george=rbind(george.15.scaled,lung.scaled)
  output.george=rbind(george.15,lung)
  output.scaled.george$size_n=ifelse(output.scaled.george[,comp.x] >=threshold,5,3 )
 # output.george$size_n=ifelse(output.george[,comp.x] >=threshold,5,3 )
  
  list_of_dataframes_scaled=lapply(list_of_dataframes,function (x) {
    blah=cbind (x[,1],scale(x[,2:(2+components-1)]),x[,(2+components),drop=F])
    colnames(blah)[[1]]="Sample"
    blah
  })
  output.tcga<-as.data.frame(rbindlist(list_of_dataframes))
  output.scaled.tcga<-as.data.frame(rbindlist(list_of_dataframes_scaled))
  
  output.tcga.length=ncol(output.tcga)
  output.tcga=output.tcga[,c(1,2,ncol(output.tcga))]
  output.scaled.tcga=output.scaled.tcga[,c(1,2,ncol(output.scaled.tcga))]
  output.scaled.tcga$size_n=ifelse(output.scaled.tcga[,comp.x] >=threshold,5,3 )
  output.tcga$size_n=ifelse(output.tcga[,comp.x] >=threshold,5,3 )
  sorting <- paste0('desc(', comp.x, ')') 
  
  by.top.3=output.scaled.tcga %>% dplyr::group_by(type) %>% dplyr::arrange_("type",sorting) %>% dplyr::slice(1:3) %>% dplyr::summarise_at(2,mean) %>% dplyr::arrange_at(2,desc)
  by.top.3= unique(by.top.3$type)
  
  by.top.3.nonscaled=output.tcga %>% dplyr::group_by(type) %>% dplyr::arrange_("type",sorting) %>% dplyr::slice(1:3) %>% dplyr::summarise_at(2,mean) %>% dplyr::arrange_at(2,desc)
  by.top.3.nonscaled= unique(by.top.3.nonscaled$type)
  
  #bymean=levels(with(output.scaled.tcga, reorder(type, -comp.x, mean)))
  bymax=unique(output.scaled.tcga[order(output.scaled.tcga[,comp.x],decreasing = T),]$type)
  bymax.notscaled=unique(output.tcga[order(output.tcga[,comp.x],decreasing = T),]$type)
  col_vector=c("darkseagreen4","indianred1","violetred","tomato2","darkcyan","springgreen3","steelblue1","skyblue","brown","royalblue3","rosybrown3",
               "red3","purple2","darkblue","limegreen","deepskyblue","darkviolet","purple4","burlywood1","magenta1","maroon3","indianred3",
               "cornflowerblue","darkorange1","darkseagreen2","dodgerblue1","firebrick2","cadetblue3","purple2","blue1","chocolate",
               "darkslategray1","aquamarine","palevioletred","chartreuse2","forestgreen","orange3","darkolivegreen2","mediumpurple4","gold","indianred2","springgreen1","darkorange3","coral4")
  set.seed(79)
  library(egg)
  
  ##### Z-scored all cancers
  plot1=ggplot(data=output.scaled.beltran,aes(x = factor(type,levels=c("NEPC","CRPC")), y = output.scaled.beltran[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(shape=19,width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  c(col_vector[1:2],"black"))+ labs(x="",y="Neuroendocrine Score")+scale_x_discrete(labels=c("NEPC","CRPC-Adeno"))+coord_cartesian(ylim = c(-5, 10)) 
  plot2=ggplot(data=output.scaled.george,aes(x = factor(type,levels=c("SCLC","NSCLC")), y = output.scaled.george[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+ 
    geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[3:4])+ labs(x="",y=element_blank()) +coord_cartesian(ylim = c(-5, 10))+
    scale_x_discrete(labels=c("SCLC","NSCLC"))
  plot3=ggplot(data=output.scaled.tcga,aes(x = factor(type,levels=bymax), y = output.scaled.tcga[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+
    theme_bw()+theme(axis.text=element_text(size=28,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[-c(1:4)])+ labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) 
  plot4=ggplot(data=output.scaled.tcga,aes(x = factor(type,levels=by.top.3), y = output.scaled.tcga[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+
    theme_bw()+theme(axis.text=element_text(size=28,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[-c(1:4)])+ labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) 
  
  png(paste0(output_folder,"TCGA_zscore_with_color_bymax.png"),width=2000,height=900)
  ggarrange(plot1,plot2, plot3, widths = c(1/12,1/12, 10/12))
  dev.off()
  png(paste0(output_folder,"TCGA_zscore_with_color_bytop3.png"),width=2000,height=900)
  ggarrange(plot1,plot2, plot4, widths = c(1/12,1/12, 10/12))
  dev.off()
  
 
  
  ##### Z-scored epithelial cancers
  output.scaled.tcga.epithelial=output.scaled.tcga[output.scaled.tcga$type %in% epithelial,]
  bymax.scaled.tcga.epithelial=unique(output.scaled.tcga.epithelial[order(output.scaled.tcga.epithelial[,comp.x],decreasing = T),]$type)
  by.top.3.scale.tcga.epithelial=output.scaled.tcga.epithelial %>% dplyr::group_by(type) %>% dplyr::arrange_("type",sorting) %>% dplyr::slice(1:3) %>% dplyr::summarise_at(2,mean) %>% dplyr::arrange_at(2,desc)
  by.top.3.scale.tcga.epithelial= unique(by.top.3.scale.tcga.epithelial$type)
  
  
  plot1=ggplot(data=output.scaled.beltran,aes(x = factor(type,levels=c("NEPC","CRPC")), y = output.scaled.beltran[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(shape=19,width=0.1,aes(color=factor(type)),size=size_n)+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  c(col_vector[1:2],"black"))+ labs(x="",y="Neuroendocrine Score")+scale_x_discrete(labels=c("NEPC","CRPC-Adeno"))+coord_cartesian(ylim = c(-5, 10)) 
  plot2=ggplot(data=output.scaled.george,aes(x = factor(type,levels=c("SCLC","NSCLC")), y = output.scaled.george[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+ 
    geom_jitter(width=0.1,aes(color=factor(type)),size=size_n)+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[3:4])+ labs(x="",y=element_blank()) +coord_cartesian(ylim = c(-5, 10))+
    scale_x_discrete(labels=c("SCLC","NSCLC"))
  plot3=ggplot(data=output.scaled.tcga.epithelial,aes(x = factor(type,levels=bymax.scaled.tcga.epithelial), y = output.scaled.tcga.epithelial[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type)),size=size_n)+
    theme_bw()+theme(axis.text=element_text(size=28,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[-c(1:4)])+ labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) 
  plot4=ggplot(data=output.scaled.tcga.epithelial,aes(x = factor(type,levels=by.top.3.scale.tcga.epithelial), y = output.scaled.tcga.epithelial[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type)),size=size_n)+
    theme_bw()+theme(axis.text=element_text(size=28,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[-c(1:4)])+ labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) 
  
  png(paste0(output_folder,"TCGA_epithelial_zscore_with_color_bymax.png"),width=2000,height=900)
  ggarrange(plot1,plot2, plot3, widths = c(1/12,1/12, 10/12))
  dev.off()
  png(paste0(output_folder,"TCGA_epithelial_zscore_with_color_bytop3.png"),width=2000,height=900)
  ggarrange(plot1,plot2, plot4, widths = c(1/12,1/12, 10/12))
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  ## NOT ZSCORED ALL
  by.avg.noz.all=output.tcga %>% dplyr::group_by(type) %>% dplyr::arrange_("type",sorting) %>% dplyr::summarise_at(2,mean) %>% dplyr::arrange_at(2,desc)
  by.avg.noz.all= unique(by.avg.noz.all$type)
  
  
   plot1=ggplot(data=output.beltran,aes(x = factor(type,levels=c("NEPC","CRPC")), y = output.beltran[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(shape=19,width=0.1,aes(color=factor(type)))+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  c(col_vector[1:2],"black"))+ labs(x="",y="Neuroendocrine Score")+scale_x_discrete(labels=c("NEPC","CRPC-Adeno"))
  plot2=ggplot(data=output.george,aes(x = factor(type,levels=c("SCLC","NSCLC")), y = output.george[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+ 
    geom_jitter(width=0.1,aes(color=factor(type)))+theme_bw()+
    theme(axis.text=element_text(size=30,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[3:4])+ labs(x="",y=element_blank()) +
    scale_x_discrete(labels=c("SCLC","NSCLC"))
  plot3=ggplot(data=output.tcga,aes(x = factor(type,levels=bymax.notscaled), y = output.tcga[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type)))+
    theme_bw()+theme(axis.text=element_text(size=28,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[-c(1:4)])+ labs(x="Cancer",y=element_blank()) 
  plot4=ggplot(data=output.tcga,aes(x = factor(type,levels=by.avg.noz.all), y = output.tcga[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type)))+
    theme_bw()+theme(axis.text=element_text(size=28,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[-c(1:4)])+ labs(x="Cancer",y=element_blank()) 
  

  png(paste0(output_folder,"TCGA_no_zscore_with_color_bymax.png"),width=2000,height=900)
  ggarrange(plot1,plot2, plot3, widths = c(1/12,1/12, 10/12))
  dev.off()
  png(paste0(output_folder,"TCGA_no_zscore_with_color_byavg.png"),width=2000,height=900)
  ggarrange(plot1,plot2, plot4, widths = c(1/12,1/12, 10/12))
  dev.off()
}
