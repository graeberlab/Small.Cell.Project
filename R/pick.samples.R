#' Pick samples
#' 
#' Reads file with features(PC's) in columns and samples in rows. And an annotation file.You pick a max # of features. Using caret
#' runs a prediction iteratively from 2 to max.features. Picks best predictor with least featurs. Returns samples correctly predicted
#' leaves out those samples incorrectly predicted.
#' @param data.file Filepath/filename of data matrix
#' @param max.features the maximum number of features to use
#' @param anno.file annotation file, has headers, first column is sample names, 2nd is type
#' @param train.method the method used for training, default is leave one out cross validation, check caret for options
#' @param technique the method used for model building, default is linear discriminant analysis, check caret for options (e.g. svmLinear,
#' @importFrom utils read.delim read.table write.table
#' 
#' @export
#' 


pick.samples<-function(data.file,max.features=5,anno.file,train.method="LOOCV",technique="lda"){
  require(caret)
  find.simplest.model.with.best.accuracy <- function(x) {
    maxX <- x[1] ## simplest case 2 features
    feat=2
    for (i in (seq_along(x[-1]) + 1)) { ## i take values 2, 3, ..., length(x)
      if(x[i] > maxX){ ## do comparison with current target
        maxX <- x[i]  ## if element i is smaller than target, update target
        feat=i+1
      }
    }
    feat
  }

  data=read.delim(data.file,row.names=1)[1:max.features]
  anno=read.delim(anno.file)
  accuracy_vect=resample_vect<-rep(NA,max.features-1)
  train_control<- trainControl(method=train.method)
  for(i in 2:max.features){
    data.temp=data[,1:i]
    data.temp$type=anno[match(rownames(data.temp),anno[,1]),2]
    
    model<- train(type~., data=data.temp, trControl=train_control, method=technique)
    accuracy=as.numeric(model$results[2])
    accuracy_vect[i-1]=accuracy
    print(paste0("Using ",i," features we have ",accuracy," Accuracy"))
  }
  features=find.simplest.model.with.best.accuracy(accuracy_vect)
  data=data[,1:features]
  data.temp$type=anno[match(rownames(data.temp),anno[,1]),2]
  model<- train(type~., data=data.temp, trControl=train_control, method="lda")
  accuracy=as.numeric(model$results[2])
  print(paste0("Best model with fewest featues has ",features," features and ",accuracy," Accuracy"))
  samples= rownames(model$pred[which(model$pred$pred==model$pred$obs),])
  samples
}

