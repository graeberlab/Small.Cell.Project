#' Partial Least Squares and predict second dataset
#'
#' Builds PLS model from training dataset and Writes out loadings.
#' Predicts classifications of test dataset and returns dataframe of predictions
#'
#' @param file file for X matirx
#' @param sample.names Vector of sample names in X matrix
#' @param response.values Vector of response values in same order matching sample.names
#' @param comps number of components to compute
#' @param scale default=T
#' @param ind.names Labels the samples, default =F
#' @param file2 file for test data matrix
#' @param output_folder the output_folder to write files to
#' @param sample.names2 Vector of sample names in 2nd dataset, if needed
#' @param response.values Vector of response values in same order matching sample.names2, if available
#' @param train_string string of training data to insert in file name of predicted scores
#' @param test_string string of data being tested to insert in file name of predicted scores
#' @param TCGA default is false, true if test samples are TCGA, removes normal samples
#' @param plot_both if true plots both training and test set in color
#' @param colpalette allows you to put in a color palette of form c("#F87660", "#39B600",....etc) to manually assign colors
#' @param shape.palette allows you to put in a shape palette of form c(1, 3,....etc) to manually assign shapes
#' @importFrom mixOmics plsda plotIndiv
#' @importFrom data.table fread
#' 
#' @export
#'


PLSDA_from_file_and_predict_second_dataset = function(file, file2, sample.names, sample.type, y.response,comps = 3, scale = F, ind.names = F,
                                                      output_folder="./",train_string="",test_string="",comp.x = "comp.1", comp.y = "comp.2",TCGA=F,plot_both = T, 
                                                      colpalette = NULL, shape.palette = NULL,labels = F,legendname = "default",ellipses=F,saveplot=T,savetype=".png",w = 8, 
                                                      h = 6,do.legend=T,title="PLSDA"){
  require(mixOmics)
  data = fread(file, sep = "\t", header = T, stringsAsFactors = FALSE, 
                    quote = "",na.strings="NA")
  data=as.data.frame(data)
  data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), 
              ]
  data2 = fread(file2, sep = "\t", header = T, stringsAsFactors = FALSE, 
                     quote = "",na.strings="NA")
  print(paste0("Done reading ",file2))
  data2=as.data.frame(data2)
  data2=na.omit(data2)
  data = data[!duplicated(data[, 1]), ]
  data2 = data2[!duplicated(data2[, 1]), ]
  if (TCGA == T) {
    temp_name = colnames(data2)[1]
    cancer_samples = which(as.numeric(sapply(colnames(data2)[-1], 
                                             function(x) strsplit(x, "\\.")[[1]][4])) <= 9)
    data2 = cbind(data2[, 1], data2[, -1][, cancer_samples])
    colnames(data2)[1] = temp_name
  }
  common.genes = intersect_all(data[, 1], data2[, 1])
  data = data[data[, 1] %in% common.genes, ]
  data2 = data2[data2[, 1] %in% common.genes, ]
  data = data[order(data[, 1]), ]
  data2 = data2[order(data2[, 1]), ]
  
  rownames(data) = make.names(data[, 1], unique = TRUE)
  t.data = data.frame(t(data[, -1]))
  y.response = (data.frame(y.response)[match(rownames(t.data), 
                                             as.character(sample.names)), ])
  y.response = as.factor(y.response)
  pls.fit = mixOmics::plsda(X = t.data, Y = y.response, scale = scale, ncomp = comps)
  plotIndiv(pls.fit, legend = T,ind.names = ind.names)
  x.variates = data.frame(pls.fit$variates$X) %>% tibble::rownames_to_column('sample')
  x.loadings = data.frame(pls.fit$loadings$X)  %>% tibble::rownames_to_column('loadings')
  
  
  write.table(x.loadings,paste0(output_folder,train_string,  "_PLSDA_Xloadings.txt"), sep = "\t", row.names = F, 
              quote = F)
  write.table(x.variates,paste0(output_folder,train_string,  "_PLSDA_XScores.txt"), sep = "\t", row.names = F,
              quote = F)
  
  x.variates$type = sample.type[match(x.variates[,1],sample.names)]
  pc.pred = ggplot(data = x.variates, aes_string(x = comp.x, 
                                                 y = comp.y)) + geom_point(size = I(2), aes(color = factor(type))) + 
    theme(legend.position = "right", plot.title = element_text(size = 30), 
          legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
          axis.title = element_text(size = 30), legend.background = element_rect(), 
          axis.text.x = element_text(margin = margin(b = -2)), 
          axis.text.y = element_text(margin = margin(l = -14))) + 
    labs(title = title) + theme_bw() + if (labels == TRUE) {
      geom_text(data = x.variates, mapping = aes(label = (rownames(x.variates))), 
                check_overlap = TRUE, size = 2.5)
    }
  pc.pred
  rownames(data2) = make.names(data2[, 1], unique = TRUE)
  t.data2 = data.frame(t(data2[,-1])) 
  
  test.predict <- predict(pls.fit, t.data2, method = "max.dist")
  test.predict.scores=as.data.frame(test.predict$variates) %>% tibble::rownames_to_column()
  colnames(test.predict.scores)[1]="sample"
  colnames(test.predict.scores)[2:ncol(test.predict.scores)]=  paste0(rep( "comp.",(ncol(test.predict.scores)-1)),1:(ncol(test.predict.scores)-1))
  write.table(test.predict.scores,paste0(output_folder,test_string,"_projected_onto_",train_string,"_",comps,"_comps_PLSDA_Xscores.txt"),col.names=T,quote=F,sep="\t",row.names=F)
  test.predict.scores$type = sample.type[match(test.predict.scores[,1],sample.names)]
  
  
  
  prediction <- as.data.frame(test.predict$class$max.dist[, comps])
  prediction <- prediction %>% tibble::rownames_to_column()
  colnames(prediction)= c("sample","prediction")
  prediction$type = sample.type[match(prediction$sample,sample.names)]
  write.table(prediction,paste0(output_folder,test_string,"_projected_onto_",train_string,"_",comps,"_comps_PLSDA_prediction.txt"),col.names=T,quote=F,sep="\t",row.names=F)
  print(table(prediction$prediction, prediction$type))
  
    
    pc.pred <- ggplot(test.predict.scores, aes_string(x = comp.x, y = comp.y)) + 
      geom_point(size = I(2), aes(color = factor(type))) + 
      theme(legend.position = "right", plot.title = element_text(size = 30), 
            legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
            axis.title = element_text(size = 30), legend.background = element_rect(), 
            axis.text.x = element_text(margin = margin(b = -2)), 
            axis.text.y = element_text(margin = margin(l = -14))) + 
      guides(color = guide_legend(title = "Type")) + labs(title = title) + 
      theme_bw() + if (labels == TRUE) {
        geom_text(data = test.predict.scores, mapping = aes(label = (rownames(test.predict.scores))), 
                  check_overlap = TRUE, size = 2.3)
      }
    pc.pred
    if (plot_both == T) {
      comb = rbind(test.predict.scores, x.variates)
      comb$type=as.factor(as.character(comb$type))
      
      if(do.legend==F) {
        pc.pred3 = ggplot(data = comb, aes_string(x = comp.x,  y = comp.y) ) + geom_point(size = I(3), aes(color = factor(type),shape = factor(type))) +
          theme(legend.position = "none",  plot.title = element_text(size = 30), legend.text = element_text(size = 22),legend.title = element_text(size = 20), axis.title = element_text(size = 30), 
                legend.background = element_rect(), axis.text.x = element_text(margin = margin(b = -2)),axis.text.y = element_text(margin = margin(l = -14))) + 
          labs(title = title) 
        
      }else{
        pc.pred3 = ggplot(data = comb, aes_string(x = comp.x,  y = comp.y) ) + geom_point(size = I(3), aes(color = factor(type), shape = factor(type))) + 
          theme(legend.position = "right", plot.title = element_text(size = 30), legend.text = element_text(size = 2), legend.title = element_text(size = 20), axis.title = element_text(size = 30), 
                legend.background = element_rect(), axis.text.x = element_text(margin = margin(b = -2)), axis.text.y = element_text(margin = margin(l = -14))) + 
          labs(title = title) +theme_bw() +theme(legend.text = element_text(size = 8))+guides(shape=F)# +theme(legend.position="none")
      }
      
      if (labels ==  TRUE) {
        pc.pred3<-pc.pred3 +geom_text(data = comb, mapping = aes(label = (rownames(comb))), 
                                      check_overlap = TRUE, size = 2.5)
      }
    
      if (!is.null(shape.palette)) {
        pc.pred3 <- pc.pred3 + scale_shape_manual(legendname, values = shape.palette)
      }
      if (!is.null(colpalette)) {
        pc.pred3 <- pc.pred3 + scale_color_manual(legendname, values = colpalette)
      }
      if (ellipses == T) {
        pc.pred3 <- pc.pred3 + stat_ellipse(aes(color = factor(type)), 
                                            level = conf)
      }
      
      if (saveplot == T) {
        ggsave(paste0(output_folder, test_string, "_projected_onto_", 
                      train_string, "_", comp.x, "_vs_", comp.y,"_PLSDA", savetype), 
               dpi = 300, plot = pc.pred3, width = w, height = h)
      }
      pc.pred3
    } else {
      pc.pred = pc.pred + geom_point(data = x.variates, aes_string(x = comp.x,
                                                                   y = comp.y)) + geom_point(size = I(1.3), aes(color = factor(type))) +
        theme(legend.position = "right", plot.title = element_text(size = 30),
              legend.text = element_text(size = 22), legend.title = element_text(size = 20),
              axis.title = element_text(size = 30), legend.background = element_rect(),
              axis.text.x = element_text(margin = margin(b = -2)),
              axis.text.y = element_text(margin = margin(l = -14))) +
        labs(title = title) + theme_bw()
      if (saveplot == T) {
        ggsave(paste0(output_folder, test_string, "_projected_onto_",
                      train_string, "_", comp.x, "_vs_", comp.y,"_PLSDA", savetype),
               dpi = 300, plot = pc.pred, width = w, height = h)
      }
      pc.pred
    }
  
}


