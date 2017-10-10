\name{PLSR_from_file_and_predict_second_dataset_graph}
\alias{PLSR_from_file_and_predict_second_dataset_graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PLSR from file and then project then output a graph
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PLSR_from_file_and_predict_second_dataset_graph(file, file2, sample.names, sample.type, y.response, sample.names2 = NULL, sample.type2 = NULL, train_string, title = "PLSR", comp.x = "comp.1", comp.y = "comp.2", comps = 3, scale = F, labels = F, output_folder = "./", test_string)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{file2}{
%%     ~~Describe \code{file2} here~~
}
  \item{sample.names}{
%%     ~~Describe \code{sample.names} here~~
}
  \item{sample.type}{
%%     ~~Describe \code{sample.type} here~~
}
  \item{y.response}{
%%     ~~Describe \code{y.response} here~~
}
  \item{sample.names2}{
%%     ~~Describe \code{sample.names2} here~~
}
  \item{sample.type2}{
%%     ~~Describe \code{sample.type2} here~~
}
  \item{train_string}{
%%     ~~Describe \code{train_string} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{comp.x}{
%%     ~~Describe \code{comp.x} here~~
}
  \item{comp.y}{
%%     ~~Describe \code{comp.y} here~~
}
  \item{comps}{
%%     ~~Describe \code{comps} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{output_folder}{
%%     ~~Describe \code{output_folder} here~~
}
  \item{test_string}{
%%     ~~Describe \code{test_string} here~~
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
function (file, file2, sample.names, sample.type, y.response, 
    sample.names2 = NULL, sample.type2 = NULL, train_string, 
    title = "PLSR", comp.x = "comp.1", comp.y = "comp.2", comps = 3, 
    scale = F, labels = F, output_folder = "./", test_string) 
{
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
    data = data[data[, 1] \%in\% common.genes, ]
    data2 = data2[data2[, 1] \%in\% common.genes, ]
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
            check_overlap = TRUE, size = 3)
    }
    pc.pred
    rownames(data2) = make.names(data2[, 1], unique = TRUE)
    t.data2 = data.frame(t(data2[, -1]))
    test.predict <- predict(pls.fit, t.data2)
    prediction <- as.data.frame(test.predict$variates)
    colnames(prediction) <- colnames(x.variates)[-ncol(x.variates)]
    write.table(cbind(Sample = rownames(prediction), (prediction)), 
        paste0(output_folder, test_string, "_projected_onto_", 
            train_string, "_PLSR_predicted.scores.txt"), sep = "\t", 
        row.names = F, quote = F)
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
            check_overlap = TRUE, size = 3)
    }
    pc.pred
    pc.pred = pc.pred + geom_point(data = x.variates, aes_string(x = comp.x, 
        y = comp.y)) + geom_point(size = I(2), aes(color = factor(type))) + 
        theme(legend.position = "right", plot.title = element_text(size = 30), 
            legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
            axis.title = element_text(size = 30), legend.background = element_rect(), 
            axis.text.x = element_text(margin = margin(b = -2)), 
            axis.text.y = element_text(margin = margin(l = -14))) + 
        labs(title = title) + theme_bw()
    png(paste0(output_folder, test_string, "_projected_onto_", 
        train_string, ".png"))
    print(pc.pred)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
