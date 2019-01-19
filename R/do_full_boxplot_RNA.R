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
  lung$type="LUAD"
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
  col_vector1 = c("#F8766D","#00BFC4","black")  #NEPC CRPC
  col_vector2 = c("#B79F00","#F564E3","black")
  col_vector3=c("#556B2F","springgreen3","steelblue1","skyblue","brown","#9590FF","rosybrown3",
                "red3","purple2","darkblue","limegreen","deepskyblue","darkviolet","purple4","burlywood1","magenta1","maroon3","indianred3",
                "cornflowerblue","darkorange1","darkseagreen2","dodgerblue1","firebrick2","cadetblue3","purple2","blue1","chocolate",
                "darkslategray1","aquamarine","palevioletred","chartreuse2","forestgreen","orange3","darkolivegreen2","mediumpurple4","gold","indianred2","springgreen1","darkorange3","coral4")
  
  
  #colpalette.final <- c("PRAD.norm"="#619CFF", "LUAD.norm"="#00BA38", "NEPC"="#00BFC4", "SCLC"="#F564E3", 
  #                "CRPC"="#F8766D", "LUAD"="#B79F00", "PRAD"="#9590FF", "BLCA.norm"="#DC143C", "BLCA"="#556B2F", "BLCA.SCN"="#A52A2A")
  col_vector4=c("LUAD"="#B79F00", "BLCA"="#556B2F","BRCA" ="springgreen3","PAAD"="steelblue1","THCA"="skyblue","KIRC"="brown",
                "PRAD"="#9590FF","STAD"="rosybrown3","HNSC"="red3","LIHC"="purple2","CESC"="darkblue",
                "COADREAD"="limegreen","UCEC"="deepskyblue","LUSC"="darkviolet","KIRP"="purple4",
                "OV"="burlywood1","ESCA"="magenta1","KICH"="maroon3","UCS"="indianred3",
                "CHOL"="cornflowerblue","ACC"="darkorange1","NBL.TARGET"="darkseagreen2","WT.TARGET"="dodgerblue1",
                "ALL.TARGET"="firebrick2","LGG"="cadetblue3","PCPG"="purple2","TGCT"="blue1","GBM"="chocolate",
                "LAML"="darkslategray1","DLBC"="aquamarine","THYM"="palevioletred","SKCM"="chartreuse2",
                "UVM"="forestgreen","SARC"="orange3","MESO"="darkolivegreen2","AML.TARGET"="mediumpurple4",
                "AMLIF.TARGET"="gold","indianred2","ACC"="springgreen1","darkorange3","coral4")
  
  
  
  set.seed(79)
  library(egg)
  
  ##### Z-scored all cancers
  plot1=ggplot(data=output.scaled.beltran,aes(x = factor(type,levels=c("NEPC","CRPC")), y = output.scaled.beltran[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(shape=19,width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=27,face="bold",angle=90,hjust=1),axis.title=element_text(size=30,face="bold"),
          axis.title.y = element_text(size = 37,face="bold",hjust=0.5,vjust=1),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  c(col_vector1))+ labs(x="",y="Small Cell Neuroendocrine       \n(SCN) Score")+
    scale_x_discrete(labels=c("NEPC","CRPC-Adeno"))+coord_cartesian(ylim = c(-5, 10))+ggpubr::rotate_x_text(65) 
  plot2=ggplot(data=output.scaled.george,aes(x = factor(type,levels=c("SCLC","LUAD")), y = output.scaled.george[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+ 
    geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=27,face="bold",angle=90,hjust=1),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector2)+ labs(x="",y=element_blank()) +coord_cartesian(ylim = c(-5, 10))+
    scale_x_discrete(labels=c("SCLC","LUAD"))+ggpubr::rotate_x_text(65)
  
  output.scaled.tcga3= output.scaled.tcga
  output.scaled.tcga3$type=factor(x=output.scaled.tcga3$type,levels=bymax)
  output.scaled.tcga3=output.scaled.tcga3[order(output.scaled.tcga3$type),]
  plot3=ggplot(data=output.scaled.tcga,aes(x = factor(type,levels=bymax), y = output.scaled.tcga[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+
    theme_bw()+theme(axis.text=element_text(size=27,face="bold",angle=90,hjust=1),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector3)+ 
    labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) +ggpubr::rotate_x_text(65)
  
  output.scaled.tcga4= output.scaled.tcga
  output.scaled.tcga4$type=factor(x=output.scaled.tcga4$type,levels=by.top.3)
  output.scaled.tcga4=output.scaled.tcga4[order(output.scaled.tcga4$type),]
  plot4=ggplot(data=output.scaled.tcga4,aes(x = factor(type,levels=by.top.3), y = output.scaled.tcga4[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+
    theme_bw()+theme(axis.text=element_text(size=28,face="bold",angle=90,hjust=1),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector3)+ labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) +ggpubr::rotate_x_text(65)
  
  png(paste0(output_folder,"TCGA_zscore_with_color_bymax.png"),width=1800,height=900)
  ggarrange(plot1,plot2, plot3, widths = c(1/12,1/12, 10/12))
  dev.off()
  png(paste0(output_folder,"TCGA_zscore_with_color_bytop3.png"),width=1600,height=600)
  ggarrange(plot1,plot2, plot4, widths = c(1/12,1/12, 10/12))
  dev.off()
  
  
  
  ##### Z-scored epithelial cancers
  output.scaled.tcga.epithelial=output.scaled.tcga[output.scaled.tcga$type %in% epithelial,]
  bymax.scaled.tcga.epithelial=unique(output.scaled.tcga.epithelial[order(output.scaled.tcga.epithelial[,comp.x],decreasing = T),]$type)
  by.top.3.scale.tcga.epithelial=output.scaled.tcga.epithelial %>% dplyr::group_by(type) %>% dplyr::arrange_("type",sorting) %>% dplyr::slice(1:3) %>% dplyr::summarise_at(2,mean) %>% dplyr::arrange_at(2,desc)
  by.top.3.scale.tcga.epithelial= unique(by.top.3.scale.tcga.epithelial$type)
  
  
  plot1=ggplot(data=output.scaled.beltran,aes(x = factor(type,levels=c("NEPC","CRPC")), y = output.scaled.beltran[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(shape=19,width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=36,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  c(col_vector1,"black"))+ labs(x="",y="Neuroendocrine Score")+
    scale_x_discrete(labels=c("NEPC","CRPC-Adeno"))+coord_cartesian(ylim = c(-5, 10))+
    ggpubr::rotate_x_text(65) 
  
  
  plot2=ggplot(data=output.scaled.george,aes(x = factor(type,levels=c("SCLC","LUAD")), y = output.scaled.george[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+ 
    geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+theme_bw()+
    theme(axis.text=element_text(size=36,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector2)+ labs(x="",y=element_blank()) +coord_cartesian(ylim = c(-5, 10))+
    scale_x_discrete(labels=c("SCLC","LUAD"))+
    ggpubr::rotate_x_text(65)
  
  plot3=ggplot(data=output.scaled.tcga.epithelial,aes(x = factor(type,levels=bymax.scaled.tcga.epithelial), y = output.scaled.tcga.epithelial[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+
    theme_bw()+theme(axis.text=element_text(size=36,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector[-c(1:4)])+ 
    labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) +
    ggpubr::rotate_x_text(65)
  
  plot4=ggplot(data=output.scaled.tcga.epithelial,aes(x = factor(type,levels=by.top.3.scale.tcga.epithelial), y = output.scaled.tcga.epithelial[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(width=0.1,aes(color=factor(type),size=size_n))+
    theme_bw()+theme(axis.text=element_text(size=36,face="bold",angle=90),axis.title=element_text(size=30,face="bold"),legend.position="none",,panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector3)+ labs(x="Cancer",y=element_blank()) +coord_cartesian(ylim = c(-5, 10)) +
    ggpubr::rotate_x_text(65)
  
  png(paste0(output_folder,"TCGA_epithelial_zscore_with_color_bymax.png"),width=2000,height=900)
  ggarrange(plot1,plot2, plot3, widths = c(1/12,1/12, 10/12),ncol=3)
  dev.off()
  png(paste0(output_folder,"TCGA_epithelial_zscore_with_color_bytop3.png"),width=2000,height=800)
  ggarrange(plot1,plot2, plot4, widths = c(1/12,1/12, 10/12),ncol=3)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  ## NOT ZSCORED ALL
  by.avg.noz.all=output.tcga %>% dplyr::group_by(type) %>% dplyr::arrange_("type",sorting) %>% dplyr::summarise_at(2,mean) %>% dplyr::arrange_at(2,desc)
  by.avg.noz.all= unique(by.avg.noz.all$type)
  
  
  plot1=ggplot(data=output.beltran,aes(x = factor(type,levels=c("NEPC","CRPC")), y = output.beltran[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  
    geom_jitter(shape=19,width=0.1,aes(color=factor(type)))+theme_bw()+
    theme(axis.text.y=element_text(size=30,face="bold",angle=90,hjust=1),
          axis.text.x=element_text(size=30,face="bold",angle=90,hjust=0.5),
          axis.title.y=element_text(size=48,face="bold",vjust=12,margin = margin(l=40),hjust=0.5),
          legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector1)+ labs(x="",y="Small Cell Neuroendocrine\n(SCN)Score")+
    scale_x_discrete(labels=c("NEPC","CRPC-Adeno"))
  
  plot2=ggplot(data=output.george,aes(x = factor(type,levels=c("SCLC","LUAD")), y = output.george[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+ 
    geom_jitter(width=0.1,aes(color=factor(type)))+theme_bw()+
    theme(axis.text.y=element_text(size=30,face="bold",angle=90,hjust=1),
          axis.text.x=element_text(size=30,face="bold",angle=90,hjust=0.5),
          axis.title=element_text(size=30,face="bold"),legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector2)+ labs(x="",y=element_blank()) +
    scale_x_discrete(labels=c("SCLC","LUAD"))
  
  plot3=ggplot(data=output.tcga,aes(x = factor(type,levels=bymax.notscaled), y = output.tcga[,comp.x]))+geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+  geom_jitter(width=0.1,aes(color=factor(type)))+
    theme_bw()+theme(axis.text=element_text(size=30,face="bold",angle=90,hjust=1),
                     axis.text.x=element_text(size=30,face="bold",angle=90,hjust=0.5),
                     axis.title=element_text(size=34,face="bold"),
                     legend.position="none",panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector3)+ labs(x="Cancer",y=element_blank()) 
 
  
  output.tcga4= output.tcga
  output.tcga4$type=factor(x=output.tcga4$type,levels=by.avg.noz.all)
  output.tcga4=output.tcga4[order(output.tcga4$type),]
  
  
   plot4=ggplot(data=output.tcga4,aes(x = factor(type,levels=by.avg.noz.all), y = output.tcga4[,comp.x]))+
     geom_boxplot(lwd=1.25,outlier.color=NA,aes(color=factor(type)))+
     geom_jitter(width=0.1,aes(color=factor(type)))+theme_bw()+
     theme(axis.text.y=element_text(size=32,face="bold",angle=90,hjust=1),
           axis.text.x=element_text(size=32,face="bold",angle=90,vjust=0.5,hjust=1),
           axis.title=element_text(size=34,face="bold"),legend.position="none",
           panel.border = element_rect(colour = "black", fill=NA, size=2))+
    scale_color_manual(values =  col_vector4)+ labs(x="Cancer",y=element_blank()) 
  
  png(paste0(output_folder,"TCGA_no_zscore_with_color_bymax.png"),width=10.75,height=3.5,units = "in",res = 600)
  ggarrange(plots=list(plot1,plot2, plot3), widths = c(1/12,1/12, 10/12),ncol=3,align="h")
  dev.off()
  png(paste0(output_folder,"TCGA_no_zscore_with_color_byavg.png"),width=1950,height=900)
  ggarrange(plots=list(plot1,plot2, plot4), widths = c(1/12,1/12, 10/12),ncol=3,align="h")
  dev.off()
  
  
  
}
