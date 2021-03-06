#' MINT Partial Least Squares and predict second dataset
#'
#' Builds MINT PLS model from training dataset and predicts second dataset
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
#' @param study.train.names vector sample types  that are in group 1  in the training studies e.g. c("NEPC","CRPC") (max 2)
#' @param study.test.names a string of the sample type that is in the test group , e.g. "PAAD" 
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
#' @param threshold threshold for predictions
#'
#' @importFrom mixOmics pls
#' @export
#'
MINT.PLSR_from_file_and_predict_second_dataset<-function (..., file2 ,sample.names, sample.type, y.response, 
                                                         sample.names2 = NULL, sample.type2 = NULL, train_string, 
                                                         test_string, title = "PLSR", comp.x = "comp.1", comp.y = "comp.2", 
                                                         comps = 2, labels = F, saveplot = T, savetype = ".png", w = 8, 
                                                         h = 6, legendname = "default", scale = F, plot_both = T, 
                                                         colpalette = NULL, shape.palette = NULL, ellipses = T, conf = 0.9, 
                                                         varimax = F, varimax.comp = 2, output_folder = "./", TCGA = F,study.train.names,study.test.names,threshold=3) {
  require(mixOmics)
  datasets=list(...)
  list_of_dataframes<-lapply(datasets, read.table,sep = "\t",header=T,stringsAsFactors = FALSE,quote="")
  data=Reduce( function(x,y) inner_join(x,y,by="gene"), list_of_dataframes)
  
  data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), 
              ]
  data2 = read.table(file2, sep = "\t", header = T, stringsAsFactors = FALSE, 
                     quote = "")
  if (TCGA == T) {
    temp_name = colnames(data2)[1]
    cancer_samples = which(as.numeric(sapply(colnames(data2)[-1], 
                                             function(x) strsplit(x, "\\.")[[1]][4])) <= 9)
    data2 = cbind(data2[, 1], data2[, -1][, cancer_samples])
    colnames(data2)[1] = temp_name
  }
  data = data[!duplicated(data[, 1]), ]
  data2 = data2[!duplicated(data2[, 1]), ]
  common.genes = intersect_all(data[, 1], data2[, 1])
  data = data[data[, 1] %in% common.genes, ]
  data2 = data2[data2[, 1] %in% common.genes, ]
  data = data[order(data[, 1]), ]
  data2 = data2[order(data2[, 1]), ]
  combined.data=Reduce( function(x,y) inner_join(x,y,by="gene"), list(data,data2))
  
  rownames(data) = make.names(data[, 1], unique = TRUE)
  t.data = data.frame(t(data[, -1]))
  train.type.vector=sample.type[match(rownames(t.data),as.character(sample.names))]
  my.train.study=ifelse(train.type.vector %in% study.train.names,1,2)
  y.response = (data.frame(y.response)[match(rownames(t.data), 
                                             as.character(sample.names)), ])
  y.response = as.matrix(y.response)
  pls.fit = mint.pls(X = t.data, Y = y.response, scale = scale, 
                     ncomp = comps,study=my.train.study)
  
  x.variates = data.frame(pls.fit$variates$X)
  x.loadings = data.frame(pls.fit$loadings$X)
  
  if (varimax == T) {
    rotation = varimax(as.matrix(x.loadings[, c(1:(varimax.comp))]), 
                       normalize = F)
    scores <- as.matrix(x.variates[, 1:(varimax.comp)]) %*% 
      rotation$rotmat
    scores = as.data.frame(scores)
    colnames(scores) = colnames(x.variates)
    x.variates = scores
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
  test.type.vector=sample.type[match(rownames(t.data2),as.character(sample.names))]
  my.test.study=ifelse(test.type.vector %in% study.test.names,3,4)
  test.predict <- predict(pls.fit, t.data2,study.test=my.test.study,dist="max.dist")
  prediction <- as.data.frame(test.predict$variates)
  colnames(prediction) <- colnames(x.variates)[-ncol(x.variates)]
  if (varimax == T) {
    predit <- as.matrix(prediction[, 1:(varimax.comp)]) %*% 
      rotation$rotmat
    colnames(predit) = colnames(prediction)
    prediction = as.data.frame(predit)
  }
  write.table(cbind(Sample = rownames(prediction), (prediction)), 
              paste0(output_folder, test_string, "_projected_onto_", 
                     train_string, "_MINT.PLSR_predicted.scores.txt"), sep = "\t", 
              row.names = F, quote = F)
  prediction.to.write=scale(prediction)
  prediction.to.write=as.data.frame(prediction.to.write)
  prediction.to.write$type = sample.type[match(rownames(prediction.to.write),sample.names)]
  prediction.to.write$prediction =ifelse(prediction.to.write$comp.1 >=threshold, 1,0)
  write.table(prediction.to.write,paste0(output_folder,test_string,"_projected_onto_",train_string,"_",comps,"_comps_MINT.PLSR_prediction.txt"),col.names=NA,quote=F,sep="\t",row.names=T)
  
  
  
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
  if (plot_both == T) {
    comb = rbind(prediction, x.variates)
    pc.pred3 = ggplot(data = comb, aes_string(x = comp.x, 
                                              y = comp.y), ) + geom_point(size = I(3), aes(color = factor(type), 
                                                                                           shape = factor(type))) + theme(legend.position = "right", 
                                                                                                                          plot.title = element_text(size = 30), legend.text = element_text(size = 22), 
                                                                                                                          legend.title = element_text(size = 20), axis.title = element_text(size = 30), 
                                                                                                                          legend.background = element_rect(), axis.text.x = element_text(margin = margin(b = -2)), 
                                                                                                                          axis.text.y = element_text(margin = margin(l = -14))) + 
      labs(title = title) + theme_bw() + if (labels == 
                                             TRUE) {
        geom_text(data = comb, mapping = aes(label = (rownames(comb))), 
                  check_overlap = TRUE, size = 2.5)
      }
    if (!is.null(shape.palette)) {
      pc.pred3 <- pc.pred3 + scale_shape_manual(legendname, 
                                                values = shape.palette)
    }
    if (!is.null(colpalette)) {
      pc.pred3 <- pc.pred3 + scale_color_manual(legendname, 
                                                values = colpalette)
    }
    if (ellipses == T) {
      pc.pred3 <- pc.pred3 + stat_ellipse(aes(color = factor(type)), 
                                          level = conf)
    }
    if (saveplot == T) {
      ggsave(paste0(output_folder, test_string, "_projected_onto_", 
                    train_string, "_", comp.x, "_vs_", comp.y, savetype), 
             dpi = 300, plot = pc.pred3, width = w, height = h)
    }
    pc.pred3
  }
  else {
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
                    train_string, "_", comp.x, "_vs_", comp.y,"_MINT_PLSR", savetype), 
             dpi = 300, plot = pc.pred, width = w, height = h)
    }
    pc.pred
  }
}