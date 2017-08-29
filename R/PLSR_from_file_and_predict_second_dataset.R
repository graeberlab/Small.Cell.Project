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
#' @param train_string string to insert in file name of predicted scores
#' @param saveplot whether to save the plot, default is F
#' @param savename the name to save plot under
#' @param savetype the type of plot to save,options are ".pdf" or ".png"
#' @param w is width of plot to be saved
#' @param h is height of plot to be saved
#' @param legendname is the legend name
#' @param plot_both if true plots both training and test set in color
#' @param colpalette allows you to put in a color palette of form c("#F87660", "#39B600",....etc) to manually assign colors
#' @param shape.palette allows you to put in a shape palette of form c(1, 3,....etc) to manually assign shapes
#' @param varimax If T performs Varimax rotation, 
#' @param varimax.comp # of varimax components, kind of hacky, keep this # the same as # of comps. Will fix later.
#' @export
#'

#sample input
# file = "filename.txt"
# sample.names = info$sample
# sample.type = info$type
# y.response = ifelse(info$type=="TYPE", 1, 0)
# scale = F
# labels = T
# comps = 3
# comp.x = "comp.1"
# comp.y = "comp.2"
# title = "Title"
# #file2 = "filename2.txt"
# sample.names2 = info$sample
# sample.type2 = info$type
# train_string = "tofile1"
PLSR_from_file_and_predict_second_dataset<-function (file, file2, sample.names, sample.type, y.response, 
                                                     sample.names2 = NULL, sample.type2 = NULL, train_string, 
                                                     title = "PLSR", comp.x = "comp.1", comp.y = "comp.2", comps = 2, 
                                                     labels = F,saveplot=F,savename="default",
                                                     savetype = ".pdf", w = 8, h = 6, legendname = "default",scale=F,
                                                     plot_both=F,colpalette=NULL,shape.palette=NULL,ellipses=T,conf=0.9,
                                                     varimax=F,varimax.comp=2) {
  require(mixOmics)
  data = read.table(file, sep = "\t", header = T, stringsAsFactors = FALSE, 
                    quote = "")
  data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), 
              ]
  data2 = read.table(file2, sep = "\t", header = T, stringsAsFactors = FALSE, 
                     quote = "")
  data = data[!duplicated(data[, 1]), ]
  data2 = data2[!duplicated(data2[, 1]), ]
  common.genes = intersect_all(data[, 1], data2[, 1])
  data = data[data[, 1] %in% common.genes, ]
  data2 = data2[data2[, 1] %in% common.genes, ]
  data = data[order(data[, 1]), ]
  data2 = data2[order(data2[, 1]), ]
  rownames(data) = make.names(data[, 1], unique = TRUE)
  t.data = data.frame(t(data[, -1]))
  y.response = (data.frame(y.response)[match(rownames(t.data), 
                                             as.character(sample.names)), ])
  y.response = as.matrix(y.response)
  pls.fit = pls(X = t.data, Y = y.response, scale = scale, 
                ncomp = comps)
  x.variates = data.frame(pls.fit$variates$X)
  x.loadings = data.frame(pls.fit$loadings$X)
  
  if(varimax==T){
    
    rotation = varimax(as.matrix(x.loadings[, c(1:(varimax.comp))]), normalize = F) 
    #rot.load = as.data.frame(matrix(unlist(rotation$loadings), nrow=nrow(x.loadings), byrow=T)) 
    #rot.load = cbind(Loading = x.loadings[,1], rot.load)
    
    
    scores <- as.matrix(x.variates[,1:(varimax.comp)]) %*% rotation$rotmat
    scores = as.data.frame(scores)
    colnames(scores)=colnames(x.variates)
    x.variates=scores
    
    #scores = cbind(Score = rownames(x.variates), scores)
    #write.table(scores, paste0(gsub(".txt", "", file.scores), "_VARIMAX.txt"), row.names = F, sep = "\t", quote = F)
    
    
  }
  
  
  
  
  
  x.variates$type = sample.type[match(rownames(x.variates), 
                                      sample.names)]
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
  t.data2 = data.frame(t(data2[, -1]))
  test.predict <- predict(pls.fit, t.data2)
  prediction <- as.data.frame(test.predict$variates)
  
  colnames(prediction) <- colnames(x.variates)[-ncol(x.variates)]
  
  
  if(varimax==T){
    predit <- as.matrix(prediction[,1:(varimax.comp)]) %*% rotation$rotmat
    colnames(predit)=colnames(prediction)
    prediction=as.data.frame(predit)
    
  }
  write.table(cbind(Sample = rownames(prediction), (prediction)), 
              paste0(gsub(".txt", "", file2), "_", train_string, "_PLSR_predicted.scores.txt"), 
              sep = "\t", row.names = F, quote = F)
  prediction$type = sample.type2[match(rownames(prediction), 
                                       sample.names2)]
  pc.pred <- ggplot(prediction, aes_string(x = comp.x, y = comp.y)) + 
    geom_point(size = I(2), aes(color = factor(type))) + 
    theme(legend.position = "right", plot.title = element_text(size = 30), 
          legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
          axis.title = element_text(size = 30), legend.background = element_rect(), 
          axis.text.x = element_text(margin = margin(b = -2)), 
          axis.text.y = element_text(margin = margin(l = -14))) + 
    guides(color = guide_legend(title = "Type")) + labs(title = title) + 
    theme_bw() + if (labels == TRUE) {
      geom_text(data = prediction, mapping = aes(label = (rownames(prediction))), 
                check_overlap = TRUE, size = 2.3)
    }
  pc.pred
  if(plot_both==T ){
    comb=rbind(prediction,x.variates)
    pc.pred3 = ggplot(data = comb, aes_string(x = comp.x, 
                                              y = comp.y),) + geom_point(size = I(3), aes(color = factor(type),shape=factor(type))) + 
      theme(legend.position = "right", plot.title = element_text(size = 30), 
            legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
            axis.title = element_text(size = 30), legend.background = element_rect(), 
            axis.text.x = element_text(margin = margin(b = -2)), 
            axis.text.y = element_text(margin = margin(l = -14))) + 
      labs(title = title) + theme_bw() +
      if (labels == TRUE) {
        geom_text(data = comb, mapping = aes(label = (rownames(comb))), 
                  check_overlap = TRUE, size = 2.5)
      }
    if(!is.null(shape.palette)){
      pc.pred3<- pc.pred3 + scale_shape_manual(legendname,values=shape.palette)
    }
    if(!is.null(colpalette)){
      pc.pred3<- pc.pred3 + scale_color_manual(legendname,values=colpalette)
    }
    if(ellipses==T){
      pc.pred3<- pc.pred3 + stat_ellipse(aes(color=factor(type)),level=conf)
    }
    
    if(saveplot==T){
      ggsave(paste0(savename, "_", comp.x, "_vs_", comp.y, savetype), 
             dpi = 300, plot = pc.pred3, width = w, height = h)
    }
    pc.pred3
  } else{
    
    pc.pred = pc.pred + geom_point(data = x.variates, aes_string(x = comp.x, 
                                                                 y = comp.y)) + geom_point(size = I(1.3), aes(color = factor(type))) + 
      theme(legend.position = "right", plot.title = element_text(size = 30), 
            legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
            axis.title = element_text(size = 30), legend.background = element_rect(), 
            axis.text.x = element_text(margin = margin(b = -2)), 
            axis.text.y = element_text(margin = margin(l = -14))) + 
      labs(title = title) + theme_bw()
    if(saveplot==T){
      ggsave(paste0(savename, "_", comp.x, "_vs_", comp.y, savetype), 
             dpi = 300, plot = pc.pred, width = w, height = h)
    }
    pc.pred
  }
}
