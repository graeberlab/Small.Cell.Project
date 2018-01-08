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
#' 
#' @importFrom mixOmics plsda plotIndiv
#' 
#' @export
#'


PLSDA_from_file_and_predict_second_dataset = function(file, file2, sample.names, sample.type, y.response, 
                                                      sample.names2 = NULL, sample.type2 = NULL, comps = 3, scale = F, ind.names = F,output_folder="./",train_string="",test_string=""){
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
  y.response = as.factor(y.response)
  pls.fit = mixOmics::plsda(X = t.data, Y = y.response, scale = scale, ncomp = comps)
  plotIndiv(pls.fit, legend = T,ind.names = ind.names)
  write.table(as.data.frame(pls.fit$loadings$X),paste0(output_folder,train_string,  "_PLSDA_Xloadings.txt"), sep = "\t", row.names = T, 
              quote = F)
  write.table(as.data.frame(pls.fit$variates$X),paste0(output_folder,train_string,  "_PLSDA_XScores.txt"), sep = "\t", row.names = T, 
              quote = F)
  rownames(data2) = data2[,1]
  t.data2 = data.frame(t(data2[,-1])) 
  
  test.predict <- predict(pls.fit, t.data2, method = "max.dist")
  test.predict.scores=as.data.frame(test.predict$variates) %>% tibble::rownames_to_column()
  colnames(test.predict.scores)[1]="sample"
  colnames(test.predict.scores)[2:ncol(test.predict.scores)]=  paste0(rep( "comp.",(ncol(test.predict.scores)-1)),1:(ncol(test.predict.scores)-1))
  write.table(test.predict.scores,paste0(output_folder,test_string,"_projected_onto_",train_string,"_",comps,"_comps_PLSDA_Xscores.txt"),col.names=T,quote=F,sep="\t",row.names=F)
  
  prediction <- as.data.frame(test.predict$class$max.dist[, comps])
  prediction <- prediction %>% tibble::rownames_to_column()
  colnames(prediction)= c("sample","prediction")
  write.table(prediction,paste0(output_folder,test_string,"_projected_onto_",train_string,"_",comps,"_comps_PLSDA_prediction.txt"),col.names=T,quote=F,sep="\t",row.names=F)

  if (!is.null(sample.names2)) {
    prediction$actual = (response.values2[match(rownames(prediction), sample.names2)])
    print(table(prediction$prediction, droplevels(prediction$actual)))
  }
  return(prediction)
  
}


