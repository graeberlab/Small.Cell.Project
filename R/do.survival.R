#' does survival analysis
#' 
#' Reads file with features(PC's) in columns and samples in rows. And an annotation file.You pick a max # of features. Using caret
#' runs a prediction iteratively from 2 to max.features. Picks best predictor with least featurs. Returns samples correctly predicted
#' leaves out those samples incorrectly predicted.
#' 
#' @param threshold threshold to use to pick samples
#' @param pattern pattern of prediction files
#' @param path location of prediction files
#' @param output_folder the output folder
#' @param component component to use
#' @param pred_file using scores or a file with predictions
#' @importFrom utils read.delim read.table write.table
#'
#' @export

do.survival=function(threshold=0,pattern,path,output_folder,component= 'comp.1',pred_file=F){
  library(valorate)
  all.files.long.path=list.files(pattern=pattern,full.names=T,path = path)
  all.files.short.path=list.files(pattern=pattern,full.names=F,path = path)
  all.names=sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1])
  list_of_dataframes<-lapply(all.files.long.path, read.table,header=T,row.names=NULL,sep = '\t',check.names=F,stringsAsFactors = F)
  names(list_of_dataframes)=all.names
  for(i in 1:length(all.files.short.path)) {
    list_of_dataframes[[i]]$type=all.names[i]
  }
  if(pred_file==F) {
  list_of_dataframes_scaled=lapply(list_of_dataframes,function (x) {
    end=ncol(x)-1
    blah=cbind (x[,1],scale(x[,2:end]),x[,ncol(x),drop=F])
    colnames(blah)[[1]]="Sample"
    blah
  })
  
  output<-as.data.frame(data.table:::rbindlist(list_of_dataframes_scaled))
  output.predictions=output
  output.predictions$prediction =ifelse(output.predictions[,component] >=threshold, 1,0)
  write.table(output.predictions,paste0(output_folder,"Predictions_threshold_",threshold,"_sd.txt"),quote=F,sep="\t",row.names=F)
  output.scaled.tcga=output[,c(1,2,ncol(output))]
  } else{
    output<-as.data.frame(data.table:::rbindlist(list_of_dataframes))
    output.scaled.tcga=output
    
  }
 
  
  if(pred_file==F){
  surv.anno=read.delim("//10.47.223.100/data2/users/nbalanis/SmallCell/Annotation/pancancer.with.typecensoring.table.epithelial.txt")
  epithelial=unique(surv.anno$type)
  output.scaled.tcga.epithelial= output.scaled.tcga[output.scaled.tcga$type %in% epithelial,]
  
  output.surv=output.scaled.tcga.epithelial %>% filter(!output.scaled.tcga.epithelial$type %in% c("PAAD","LIHC"))
  #output.surv=output.scaled.tcga.epithelial
  output.surv$prediction=0
  threshold=threshold
  output.surv$prediction =ifelse(output.surv[,component] >=threshold, 1,0)
  output.surv.continuous=output.surv[,c(1,3,2)]
  output.surv.discrete=output.surv[,c(1,3,4)]
  sum(output.surv.discrete$prediction)
  divide_by_avg<- function(x){ x/mean(x)}
  #surv.anno=read.delim("//10.47.223.100/data2/users/nbalanis/SmallCell/Annotation/pancancer.allcancers.censoring.table.txt")
  surv.anno$patientID=make.names(surv.anno$patientID)
  NE.surv.frame.discrete=inner_join(output.surv.discrete,surv.anno,by=c("Sample"="patientID"))
  NE.surv.frame.discrete$type=NE.surv.frame.discrete$type.x
  NE.surv.frame.discrete.OS=NE.surv.frame.discrete %>% dplyr::select(Sample,survival_time=OS,censoring_status=OS_IND,type,prediction) %>% na.omit(.)
  NE.surv.frame.discrete.OS.normalized=NE.surv.frame.discrete.OS %>% group_by(type) %>% mutate(survival_time=divide_by_avg(survival_time))
  NE.surv.frame.discrete.OS.normalized=data.frame(NE.surv.frame.discrete.OS.normalized)
  NE.surv.frame.discrete.RFS=NE.surv.frame.discrete %>% dplyr::select(Sample,survival_time=RFS,censoring_status=RFS_IND,type,prediction) %>% na.omit(.)
  
  NE.surv.frame.continuous=inner_join(output.surv.continuous,surv.anno,by=c("Sample"="patientID"))
  NE.surv.frame.continuous$type=NE.surv.frame.continuous$type.x
  NE.surv.frame.continuous.OS=NE.surv.frame.continuous %>% dplyr::select(Sample,survival_time=OS,censoring_status=OS_IND,type,component) %>% na.omit(.)
  NE.surv.frame.continuous.OS.normalized=NE.surv.frame.continuous.OS %>% group_by(type) %>% mutate(survival_time=divide_by_avg(survival_time))
  NE.surv.frame.continuous.OS.normalized=data.frame(NE.surv.frame.continuous.OS.normalized)
  NE.surv.frame.continuous.RFS=NE.surv.frame.continuous %>% dplyr::select(Sample,survival_time=RFS,censoring_status=RFS_IND,type,component) %>% na.omit(.)
  } else {
    surv.anno=read.delim("//10.47.223.100/data2/users/nbalanis/SmallCell/Annotation/pancancer.with.typecensoring.table.epithelial.txt")
    epithelial=unique(surv.anno$type)
    output.scaled.tcga.epithelial= output.scaled.tcga[output.scaled.tcga$type %in% epithelial,]
    output.surv=output.scaled.tcga.epithelial %>% filter(!output.scaled.tcga.epithelial$type %in% c("PAAD","LIHC")) 
    #output.surv=output.scaled.tcga.epithelial
    output.surv.discrete=output.surv
    sum(output.surv.discrete$prediction)
    divide_by_avg<- function(x){ x/mean(x)}
    #surv.anno=read.delim("//10.47.223.100/data2/users/nbalanis/SmallCell/Annotation/pancancer.allcancers.censoring.table.txt")
    surv.anno$patientID=make.names(surv.anno$patientID)
    NE.surv.frame.discrete=inner_join(output.surv.discrete,surv.anno,by=c("sample"="patientID"))
    NE.surv.frame.discrete$type=NE.surv.frame.discrete$type.x
    NE.surv.frame.discrete.OS=NE.surv.frame.discrete %>% dplyr::select(sample,survival_time=OS,censoring_status=OS_IND,type,prediction) %>% na.omit(.)
    NE.surv.frame.discrete.OS.normalized=NE.surv.frame.discrete.OS %>% group_by(type) %>% mutate(survival_time=divide_by_avg(survival_time))
    NE.surv.frame.discrete.OS.normalized=data.frame(NE.surv.frame.discrete.OS.normalized)
    NE.surv.frame.discrete.RFS=NE.surv.frame.discrete %>% dplyr::select(sample,survival_time=RFS,censoring_status=RFS_IND,type,prediction) %>% na.omit(.)
    
    
    
    
  }
  
  if(pred_file==F){
  stats=NE.surv.frame.discrete.OS %>%count(type,prediction) %>% as.data.frame(.)
  write.table(stats,paste0(output_folder,"stats_",threshold,"_","sd.txt"),quote=F,sep="\t",row.names=F)
  } else{
    stats=NE.surv.frame.discrete.OS %>%count(type,prediction) %>% as.data.frame(.)
    write.table(stats,paste0(output_folder,"stats_prediction.txt"),quote=F,sep="\t",row.names=F)
    
  }
  remove_genes_that_same2<-function(exalt_frame){
    exalt_clinical<-exalt_frame[,1:3]
    exalt_data<-exalt_frame[,4:ncol(exalt_frame),drop=F]
    remove_col <- !sapply( exalt_data, function(x) length(unique(x[!is.na(x)]))<=1 )
    exalt_out<-cbind(exalt_clinical,exalt_data[,remove_col])
    return(exalt_out)
  }
  #debug(remove_genes_that_same2)
  
  
  #valorate.p.value.sampling(blah,7,attributes(val.surv)[[1]]["LR"])
  #valorate.plot.kaplan(blah,exalt_mut_OS$C4BPA)
  
  
  graph_logrank_OS_cox<-function(data,prefx=prefix,pval_cutoff=1,pct_alt=0.005,output_folder,start.column,output.string){
    require(survival)
    # input<-as.formula(paste0("Surv(survival_time,censoring_status) ~ ",gene_name))
    # gene.survival <-survdiff(input,data=data)
    #data= data[!is.na(data$stage),]
    #data=remove_genes_that_same2(data)
    
    for(gene in colnames(data)[start.column:ncol(data)]) {
      gene_name=gene
      if(length(unique(data[,gene_name])) <=1)
      {
        return(NULL)
      }
      gene.survival <-survdiff(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      fit <-           survfit(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      # val.obj=new.valorate(time = data$survival_time,status = data$censoring_status,tails=2,sampling.size = 10000,min.sampling.size = 1000)
      #  val.surv=valorate.survdiff(val.obj,data[,gene_name])
      
      cfit= tryCatch({cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] + data$type )},
                     warning = function(w) {
                       print (w); 
                       print (gene_name)
                       cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] + data$type )
                       return(cfit)
                     },
                     error = function(e) {
                       print(e);
                       print(gene_name);
                     }
      ) #+ data$type + data$n_stage +data$m_stage +data$t_stage)  #)
      #cpval = coef(summary(cfit))['data[, paste0(sig, ".HiLo")]','Pr(>|z|)']
      #fit<-survfit(input,data=data)
      surv_gene_name<-gene_name
      capture.output(summary(cfit),file=paste0(output_folder,"pancan.summary.threshold","_",threshold,"_sd.txt"))
      p.val <-  1 - pchisq(gene.survival$chisq, length(gene.survival$n) - 1)
      cp.val =coef(summary(cfit))[1,5]
      num_not_altered=as.numeric(gene.survival$n[1])
      num_altered=as.numeric(gene.survival$n[2])
      # #pval less than 0.1,percent altered
      if(cp.val <= pval_cutoff & ((num_altered/(num_altered+num_not_altered))>=pct_alt)){
        png(paste0(output_folder,prefx,"_",surv_gene_name,"_",threshold,"_",output.string,"_cox_twoway_logrank_survival.png"))
        plot(fit, lty = c(1, 2), xlab = "Days", ylab = "Overall Survival", lwd = 2, main = paste0("Overall Survival in ",prefx))
        legend("topright",c(paste0("Non-Small-Like"," (",num_not_altered,")"),paste0("Small-Like"," (",num_altered,")") ),title=toupper(surv_gene_name))
        text(x= 0.8*max(fit$time), y=.8*(max(fit$surv)),paste0("p.val = ", cp.val))
        dev.off()
        
      }
    }
  }
  
  
  graph_logrank_OS_cox_continuous<-function(data,prefx=prefix,pval_cutoff=1,pct_alt=0.005,output_folder,start.column,output.string){
    require(survival)
    # input<-as.formula(paste0("Surv(survival_time,censoring_status) ~ ",gene_name))
    # gene.survival <-survdiff(input,data=data)
    #data= data[!is.na(data$stage),]
    #data=remove_genes_that_same2(data)
    for(gene in colnames(data)[start.column:ncol(data)]) {
      gene_name=gene
      # gene.survival <-survdiff(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      # fit <-           survfit(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      # # val.obj=new.valorate(time = data$survival_time,status = data$censoring_status,tails=2,sampling.size = 10000,min.sampling.size = 1000)
      #  val.surv=valorate.survdiff(val.obj,data[,gene_name])
      
      cfit= tryCatch({cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] + data$type )},
                     warning = function(w) {
                       print (w); 
                       print (gene_name)
                       cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] + data$type )
                       return(cfit)
                     },
                     error = function(e) {
                       print(e);
                       print(gene_name);
                     }
      ) #+ data$type + data$n_stage +data$m_stage +data$t_stage)  #)
      #cpval = coef(summary(cfit))['data[, paste0(sig, ".HiLo")]','Pr(>|z|)']
      #fit<-survfit(input,data=data)
      #surv_gene_name<-gene_name
      capture.output(summary(cfit),file=paste0(output_folder,"pancan_coxph_summary_continuous.txt"))
      # p.val <-  1 - pchisq(gene.survival$chisq, length(gene.survival$n) - 1)
      #cp.val =coef(summary(cfit))[1,5]
      #num_not_altered=as.numeric(gene.survival$n[1])
      # num_altered=as.numeric(gene.survival$n[2])
      # #pval less than 0.1,percent altered
      # if(cp.val <= pval_cutoff & ((num_altered/(num_altered+num_not_altered))>=pct_alt)){
      #   png(paste0(output_folder,prefx,"_",surv_gene_name,"_",threshold,"_",output.string,"_cox_twoway_logrank_survival.png"))
      #   plot(fit, lty = c(1, 2), xlab = "Days", ylab = "Overall Survival", lwd = 2, main = paste0("Overall Survival in ",prefx))
      #   legend("topright",c(paste0("Non-Small-Like"," (",num_not_altered,")"),paste0("Small-Like"," (",num_altered,")") ),title=toupper(surv_gene_name))
      #   text(x= 0.8*max(fit$time), y=.8*(max(fit$surv)),paste0("p.val = ", cp.val))
      #   dev.off()
      #   
      # }
    }
  }
  
  
  graph_logrank_RFS_cox<-function(data,prefx=prefix,pval_cutoff=1,pct_alt=0.005,type,output_folder){
    require(survival)
    # input<-as.formula(paste0("Surv(survival_time,censoring_status) ~ ",gene_name))
    # gene.survival <-survdiff(input,data=data)
    #data= data[!is.na(data$stage),]
    #data=remove_genes_that_same2(data)
    
    for(gene in colnames(data)[5:ncol(data)]) {
      gene_name=gene
      if(length(unique(data[,gene_name])) <=1)
      {
        return(NULL)
      }
      gene.survival <-survdiff(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      fit <-           survfit(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      
      cfit= tryCatch({cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] +data$type )},
                     warning = function(w) {
                       print (w); 
                       print (gene_name)
                       cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name]  + data$type)
                       return(cfit)
                     },
                     error = function(e) {
                       print(e);
                       print(gene_name);
                     }
      ) #+ data$type + data$n_stage +data$m_stage +data$t_stage)  #)
      #cpval = coef(summary(cfit))['data[, paste0(sig, ".HiLo")]','Pr(>|z|)']
      #fit<-survfit(input,data=data)
      surv_gene_name<-gene_name
      p.val <-  1 - pchisq(gene.survival$chisq, length(gene.survival$n) - 1)
      cp.val =coef(summary(cfit))[1,5]
      num_not_altered=as.numeric(gene.survival$n[1])
      num_altered=as.numeric(gene.survival$n[2])
      # #pval less than 0.1,percent altered
      if(cp.val <= pval_cutoff & ((num_altered/(num_altered+num_not_altered))>=pct_alt)){
        png(paste0(output_folder,prefx,"_",surv_gene_name,"_",type,"__cox_twoway_logrank_survival.png"))
        plot(fit, lty = c(1, 2), xlab = "Days", ylab = "Recurrence Free Survival", lwd = 2, main = paste0("Recurrence Free Survival in ",prefx))
        legend("topright",c(paste0("Non-Small-Like"," (",num_not_altered,")"),paste0("Small-Like"," (",num_altered,")") ),title=toupper(surv_gene_name))
        text(x= 0.8*max(fit$time), y=.8*(max(fit$surv)),paste0("p.val = ", cp.val))
        dev.off()
        
      }
    }
  }
  
  graph_logrank_OS_cox_indiv<-function(data,prefx=prefix,pval_cutoff=1,pct_alt=0.005,output_folder,output.string=""){
    require(survival)
    # input<-as.formula(paste0("Surv(survival_time,censoring_status) ~ ",gene_name))
    # gene.survival <-survdiff(input,data=data)
    #data= data[!is.na(data$stage),]
    #data=remove_genes_that_same2(data)
    if(length(unique(data[,5])) <=1) {
      return("Only one group")
    }
    for(gene in colnames(data)[5:ncol(data)]) {
      gene_name=gene
      gene.survival <-survdiff(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      fit <-           survfit(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      val.obj=new.valorate(time = data$survival_time,status = data$censoring_status,tails=2,sampling.size = 10000,min.sampling.size = 1000)
      val.surv=valorate.survdiff(val.obj,data[,gene_name])
      val.frame.length=length(attributes(val.surv))
      if(val.frame.length>1){
        val.frame=cbind( data.frame(pval=val.surv[1]),as.data.frame(t(attributes(val.surv)[[1]] )),as.data.frame(t( attributes(val.surv)[[2]])))
      } else {
        val.frame=cbind( data.frame(pval=val.surv[1]),as.data.frame(t(attributes(val.surv)[[1]] )))
      }
      write.table(val.frame,paste0(output_folder,prefx,"_valorate_threshold_",threshold,"_sd.txt"),quote=F,sep="\t",row.names=F)
      
      cfit= tryCatch({cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name])},
                     warning = function(w) {
                       print (w); 
                       print (gene_name)
                       cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
                       return(cfit)
                     },
                     error = function(e) {
                       print(e);
                       print(gene_name);
                     }
      ) #+ data$type + data$n_stage +data$m_stage +data$t_stage)  #)
      #cpval = coef(summary(cfit))['data[, paste0(sig, ".HiLo")]','Pr(>|z|)']
      #fit<-survfit(input,data=data)
      surv_gene_name<-gene_name
      capture.output(summary(cfit),file=paste0(output_folder,prefx,"_coxph_threshold_",threshold,"_sd.txt"))
      p.val <-  1 - pchisq(gene.survival$chisq, length(gene.survival$n) - 1)
      cp.val =coef(summary(cfit))[1,5]
      num_not_altered=as.numeric(gene.survival$n[1])
      num_altered=as.numeric(gene.survival$n[2])
      # #pval less than 0.1,percent altered
      
      png(paste0(output_folder,prefx,"_",surv_gene_name,"_",threshold,"_",output.string,"_cox_twoway_logrank_survival.png"))
      plot(fit, lty = c(1, 2), xlab = "Days", ylab = "Overall Survival", lwd = 2, main = paste0("Overall Survival in ",prefx))
      legend("topright",c(paste0("Not Mutated"," (",num_not_altered,")"),paste0("Mutated"," (",num_altered,")") ),title=toupper(surv_gene_name))
      text(x= 0.8*max(fit$time), y=.8*(max(fit$surv)),paste0("p.val = ", cp.val))
      dev.off()
      return(val.frame$pval)
    }
  }
  graph_logrank_OS_cox_indiv_continuous<-function(data,prefx=prefix,pval_cutoff=1,pct_alt=0.005,output_folder,output.string=""){
    require(survival)
    # input<-as.formula(paste0("Surv(survival_time,censoring_status) ~ ",gene_name))
    # gene.survival <-survdiff(input,data=data)
    #data= data[!is.na(data$stage),]
    #data=remove_genes_that_same2(data)
    if(length(unique(data[,5])) <=1) {
      return("Only one group")
    }
    for(gene in colnames(data)[5:ncol(data)]) {
      gene_name=gene
      #gene.survival <-survdiff(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      #fit <-           survfit(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
      # val.obj=new.valorate(time = data$survival_time,status = data$censoring_status,tails=2,sampling.size = 10000,min.sampling.size = 1000)
      #val.surv=valorate.survdiff(val.obj,data[,gene_name])
      #val.frame.length=length(attributes(val.surv))
      #if(val.frame.length>1){
      #  val.frame=cbind( data.frame(pval=val.surv[1]),as.data.frame(t(attributes(val.surv)[[1]] )),as.data.frame(t( attributes(val.surv)[[2]])))
      # } else {
      #   val.frame=cbind( data.frame(pval=val.surv[1]),as.data.frame(t(attributes(val.surv)[[1]] )))
      # }
      #  write.table(val.frame,paste0(output_folder,prefx,"_valorate_threshold_",threshold,"_sd.txt"),quote=F,sep="\t",row.names=F)
      
      cfit= tryCatch({cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name])},
                     warning = function(w) {
                       print (w); 
                       print (gene_name)
                       cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
                       return(cfit)
                     },
                     error = function(e) {
                       print(e);
                       print(gene_name);
                     }
      ) #+ data$type + data$n_stage +data$m_stage +data$t_stage)  #)
      #cpval = coef(summary(cfit))['data[, paste0(sig, ".HiLo")]','Pr(>|z|)']
      #fit<-survfit(input,data=data)
      # surv_gene_name<-gene_name
      capture.output(summary(cfit),file=paste0(output_folder,prefx,"_coxph_summary_continuous.txt"))
      #p.val <-  1 - pchisq(gene.survival$chisq, length(gene.survival$n) - 1)
      #cp.val =coef(summary(cfit))[1,5]
      #num_not_altered=as.numeric(gene.survival$n[1])
      #num_altered=as.numeric(gene.survival$n[2])
      # #pval less than 0.1,percent altered
      # 
      # png(paste0(output_folder,prefx,"_",surv_gene_name,"_",threshold,"_",output.string,"_cox_twoway_logrank_survival.png"))
      # plot(fit, lty = c(1, 2), xlab = "Days", ylab = "Overall Survival", lwd = 2, main = paste0("Overall Survival in ",prefx))
      # legend("topright",c(paste0("Not Mutated"," (",num_not_altered,")"),paste0("Mutated"," (",num_altered,")") ),title=toupper(surv_gene_name))
      # text(x= 0.8*max(fit$time), y=.8*(max(fit$surv)),paste0("p.val = ", cp.val))
      # dev.off()
      
      
    }
  }
  
  graph_logrank_OS_cox_pval<-function(data,prefx=prefix,pval_cutoff=1,pct_alt=0.005,type,perms=1000){
    require(survival)
    resample_vect<-rep(NA,perms)
    for(i in 1:perms){
      # input<-as.formula(paste0("Surv(survival_time,censoring_status) ~ ",gene_name))
      # gene.survival <-survdiff(input,data=data)
      #data= data[!is.na(data$stage),]
      #data=remove_genes_that_same2(data)
      data[,5]= sample(data[,5])
      for(gene in colnames(data)[5:ncol(data)]) {
        gene_name=gene
        gene.survival <-survdiff(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
        fit <-           survfit(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] )
        
        cfit= tryCatch({cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] +data$type)},
                       warning = function(w) {
                         print (w); 
                         print (gene_name)
                         cfit = coxph(Surv(data$survival_time,data$censoring_status) ~ data[,gene_name] + data$type )
                         return(cfit)
                       },
                       error = function(e) {
                         print(e);
                         print(gene_name);
                       }
        ) #+ data$type + data$n_stage +data$m_stage +data$t_stage)  #)
        #cpval = coef(summary(cfit))['data[, paste0(sig, ".HiLo")]','Pr(>|z|)']
        #fit<-survfit(input,data=data)
        surv_gene_name<-gene_name
        p.val <-  1 - pchisq(gene.survival$chisq, length(gene.survival$n) - 1)
        cp.val =coef(summary(cfit))[1,5]
        resample_vect[i]=cp.val
        
        
      }
    }
    return(resample_vect)
  }
  
  
  
  
  
  
  #undebug(graph_logrank_OS_cox)
  if(pred_file==F){
  graph_logrank_OS_cox(prefx="pancan",start.column=5,data = NE.surv.frame.discrete.OS,pct_alt=0,pval_cutoff=1,output_folder=output_folder,output.string="epithelial_discrete")
  graph_logrank_OS_cox_continuous(prefx="pancan",start.column=5,data = NE.surv.frame.continuous.OS,pct_alt=0,pval_cutoff=1,output_folder=output_folder,output.string="epithelial_continuous")
  } else{
    graph_logrank_OS_cox(prefx="pancan",start.column=5,data = NE.surv.frame.discrete.OS,pct_alt=0,pval_cutoff=1,output_folder=output_folder,output.string="epithelial_discrete")
    
  }
  
  #debug(graph_logrank_OS_cox_indiv)
  for(cancer in unique(NE.surv.frame.discrete.OS$type)){
    NE.surv.frame.OS.discrete.temp= NE.surv.frame.discrete.OS %>% filter(type==cancer)
    graph_logrank_OS_cox_indiv(prefx=cancer,data = NE.surv.frame.OS.discrete.temp,pct_alt=0,pval_cutoff=1,output_folder=output_folder,output.string="epithelial_discrete")
    
  }
  
  if(pred_file==F){
  for(cancer in unique(NE.surv.frame.continuous.OS$type)){
    NE.surv.frame.OS.continuous.temp= NE.surv.frame.continuous.OS %>% filter(type==cancer)
    graph_logrank_OS_cox_indiv_continuous(prefx=cancer,data = NE.surv.frame.OS.continuous.temp,pct_alt=0,pval_cutoff=1,output_folder=output_folder,output.string="epithelial_continous")
  }
  }
  
  
  
  
}