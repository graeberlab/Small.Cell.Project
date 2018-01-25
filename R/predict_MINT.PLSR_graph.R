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



predict_MINT.PLSR_graph=function(...,test_files_folder,anno.file,output_folder="./",train_string="",test_pattern="",comps=3,study.train.names,y.response){
  
  all.files.short.path=list.files(test_files_folder,pattern=test_pattern,full.names=F)
  all.files.long.path=list.files(test_files_folder,pattern=test_pattern,full.names=T)
  all.names=as.character(sapply(all.files.short.path, function(x) strsplit(x,"_")[[1]][1]))
  human.info=read.delim(anno.file)
  for(i in 1:length(all.names)) {
    nam=all.names[i]
    MINT.PLSR_from_file_and_predict_second_dataset(...,
                                              file2=all.files.long.path[i],
                                              sample.names= human.info$sample,sample.type=human.info$type,
                                              sample.names2 = human.info$sample,sample.type2 =  factor(human.info$type),
                                              y.response =  y.response,
                                              test_string = nam,comps = comps, scale = F,output_folder=output_folder,train_string=train_string,
                                              TCGA=T,study.train.names=study.train.names,study.test.names=nam,saveplot=T,plot_both=T,colpalette = cbbPalette4,shape.palette=shape.palt)
    print(paste0("Cancer ",nam," is done!"))
  }
}