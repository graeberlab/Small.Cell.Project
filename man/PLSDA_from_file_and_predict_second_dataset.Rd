% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLSDA_from_file_and_predict_second_dataset.R
\name{PLSDA_from_file_and_predict_second_dataset}
\alias{PLSDA_from_file_and_predict_second_dataset}
\title{Partial Least Squares and predict second dataset}
\usage{
PLSDA_from_file_and_predict_second_dataset(file, file2, sample.names,
  sample.type, y.response, comps = 3, scale = F, ind.names = F,
  output_folder = "./", train_string = "", test_string = "",
  comp.x = "comp.1", comp.y = "comp.2", TCGA = F, plot_both = T,
  colpalette = NULL, shape.palette = NULL, labels = F,
  legendname = "default", ellipses = F, saveplot = T,
  savetype = ".png", w = 8, h = 6, do.legend = T,
  title = "PLSDA")
}
\arguments{
\item{file}{file for X matirx}

\item{file2}{file for test data matrix}

\item{sample.names}{Vector of sample names in X matrix}

\item{comps}{number of components to compute}

\item{scale}{default=T}

\item{ind.names}{Labels the samples, default =F}

\item{output_folder}{the output_folder to write files to}

\item{train_string}{string of training data to insert in file name of predicted scores}

\item{test_string}{string of data being tested to insert in file name of predicted scores}

\item{TCGA}{default is false, true if test samples are TCGA, removes normal samples}

\item{plot_both}{if true plots both training and test set in color}

\item{colpalette}{allows you to put in a color palette of form c("#F87660", "#39B600",....etc) to manually assign colors}

\item{shape.palette}{allows you to put in a shape palette of form c(1, 3,....etc) to manually assign shapes}

\item{response.values}{Vector of response values in same order matching sample.names}

\item{sample.names2}{Vector of sample names in 2nd dataset, if needed}

\item{response.values}{Vector of response values in same order matching sample.names2, if available}
}
\description{
Builds PLS model from training dataset and Writes out loadings.
Predicts classifications of test dataset and returns dataframe of predictions
}
