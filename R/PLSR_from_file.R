#'Partial Least Squares Reg from file
#'
#' Writes out loadings
#'
#' @param file file for X matrix
#' @param sample.names vector of names of the samples
#' @param sample.type vector of sample groups
#' @param y.response numeric vector of response values in same order as samples
#' @param comp number of components to compute
#' @param scale default=F
#' @param labels label the plot, default = F
#' @param comp.x,comp.y comps to display
#' @param title title of the plot
#' @param varimax performs varimax 
#' 
#' @export
#'
PLSR_from_file<-function (file, sample.names, sample.type, y.response, title = "PLSR", 
          comps = 5, scale = F, comp.x = "comp.1", comp.y = "comp.2", 
          labels = F,varimax=F,varimax.comp=5) {
  require(mixOmics)
  require(ggplot2)
  data = read.table(file, sep = "\t", header = T, stringsAsFactors = FALSE, 
                    quote = "")
  data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), 
              ]
  rownames(data) = make.names(data[, 1], unique = TRUE)
  t.data = data.frame(t(data[, -1]))
  y.response = (data.frame(y.response)[match(rownames(t.data), 
                                             sample.names), ])
  y.response = as.matrix(y.response)
  pls.fit = pls(X = t.data, Y = y.response, scale = scale, 
                ncomp = comps)
  print(pls.fit$explained_variance$X)
  x.variates = data.frame(pls.fit$variates$X)
  x.loadings = data.frame(pls.fit$loadings$X)
  if(varimax==T){
  rotation = varimax(as.matrix(x.loadings[, c(1:(varimax.comp))]), 
                     normalize = F)
  scores <- as.matrix(x.variates[, 1:(varimax.comp)]) %*% 
    rotation$rotmat
  scores = as.data.frame(scores)
  colnames(scores) = colnames(x.variates)
  x.variates = scores
}
  
  x.exp_variance = data.frame(pls.fit$explained_variance$X)
  variates.X = cbind(Score = rownames(pls.fit$variates$X), 
                     x.variates)
  loadings.X = cbind(Loading = rownames(pls.fit$loadings$X), 
                     x.loadings)
  rownames(x.exp_variance) = paste0("comp.", seq(1, nrow(x.exp_variance)))
  write.table(as.data.frame(variates.X), paste0(gsub(".txt", 
                                                     "", file), "_PLSR_Xscores.txt"), sep = "\t", row.names = F, 
              quote = F)
  write.table(as.data.frame(loadings.X), paste0(gsub(".txt", 
                                                     "", file), "_PLSR_Xloadings.txt"), sep = "\t", row.names = F, 
              quote = F)
  write.table(as.data.frame(x.exp_variance), paste0(gsub(".txt", 
                                                         "", file), "_PLSR_Xpve.txt"), sep = "\t", row.names = T, 
              quote = F)
}
