#' Predict all files in a folder using PLSDA of a training dataset, file names have to be of form name_.......txt
#' Predicts classifications of test dataset and writes out predictions
#'
#' @param prediction_file file  with training data features rows , genes columns, features are in column 1
#' @param test_files_folder folder of test files
#' @param anno.file Annotation file, first column is sample names, 2nd is annotation
#' @param comps number of components to compute
#' @param output_folder output folder
#' @param train_string string of training data to insert in file name of predicted scores
#' @param test_string string of data being tested to insert in file name of predicted scores
#' @param train_pattern annotation type of main group in your comparison ,taken from 2nd column in anno.
#' 
#' @importFrom mixOmics plsda plotIndiv
#' 
#' @export
#'



predict_PLSDA_graph=function(prediction_file,test_files_folder,anno.file,output_folder="./",train_string="",test_pattern="",train_pattern="",comps=3){
  
  all.files.short.path=list.files(test_files_folder,pattern=test_pattern,full.names=F)
  all.files.long.path=list.files(test_files_folder,pattern=test_pattern,full.names=T)
  all.names=as.character(sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1]))
  human.info=read.delim(anno.file)
  for(i in 1:length(all.names)) {
    nam=all.names[i]
    PLSDA_from_file_and_predict_second_dataset(file = prediction_file,
                                              file2=all.files.long.path[i],
                                              sample.names= human.info$sample, 
                                              response.values =  ifelse(human.info$type==train_pattern,1,0),
                                              sample.names2 = human.info$sample, response.values2=factor(human.info$type),
                                              test_string = nam,
                                              "Projection", comps = comps, scale = F,output_folder=output_folder,train_string=train_string)
    print(paste0("Cancer ",nam," is done!"))
  }
}