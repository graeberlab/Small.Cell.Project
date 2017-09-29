\name{pick.samples}
\alias{pick.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pick samples 
}
\description{
Takes a data matrix of samples and features, usually PCA components are the features. Additionally takes a annotation file where samples are classified in two classes. Then using a model specified by train.method and technique. Tries to predict all individual samples class based on the data instead of the annotation file. In any case where the annotation and the data driven classification do not agree the sample is removed. Returns a list of the sample names where the annotation file and the data driven model agree.
Refines the data by perfor
}
\usage{
pick.samples(data.file, max.features = 5, anno.file, train.method = "LOOCV", technique = "lda")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.file}{
   the data file
}
  \item{max.features}{
the maximum number of features to use
}
  \item{anno.file}{
annotation file, has headers, first column is sample names, 2nd is type
}
  \item{train.method}{
 the method used for training, default is leave one out cross validation, check caret for options
}
  \item{technique}{
the method used for model building, default is linear discriminant analysis, check caret for options (e.g. svmLinear,
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data.file, max.features = 5, anno.file, train.method = "LOOCV", 
    technique = "lda") 
{
    require(caret)
    find.simplest.model.with.best.accuracy <- function(x) {
        maxX <- x[1]
        feat = 2
        for (i in (seq_along(x[-1]) + 1)) {
            if (x[i] > maxX) {
                maxX <- x[i]
                feat = i + 1
            }
        }
        feat
    }
    data = read.delim(data.file, row.names = 1)[1:max.features]
    anno = read.delim(anno.file)
    accuracy_vect = resample_vect <- rep(NA, max.features - 1)
    train_control <- trainControl(method = train.method)
    for (i in 2:max.features) {
        data.temp = data[, 1:i]
        data.temp$type = anno[match(rownames(data.temp), anno[, 
            1]), 2]
        model <- train(type ~ ., data = data.temp, trControl = train_control, 
            method = technique)
        accuracy = as.numeric(model$results[2])
        accuracy_vect[i - 1] = accuracy
        print(paste0("Using ", i, " features we have ", accuracy, 
            " Accuracy"))
    }
    features = find.simplest.model.with.best.accuracy(accuracy_vect)
    data.temp = data[, 1:features]
    data.temp$type = anno[match(rownames(data.temp), anno[, 1]), 
        2]
    model <- train(type ~ ., data = data.temp, trControl = train_control, 
        method = "lda")
    accuracy = as.numeric(model$results[2])
    print(paste0("Best model with fewest featues has ", features, 
        " features and ", accuracy, " Accuracy"))
    samples = rownames(model$pred[which(model$pred$pred == model$pred$obs), 
        ])
    samples
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
