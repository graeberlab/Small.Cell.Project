#' Create methylation files for TCGA cancers using defined sets of sites
#'
#' Takes illumina methylation files and creates new individualfiles using features that exist in the mysites files 
#'
#' @param file input methylation file , file name is of form "ACC_etc_etc...."
#' @param subset Whether to restrict to sites in mysites file,default T
#' @param mysites vector of site names
#' @param out.string If subset=T Should be the string explaining the sites used 
#'
#' @export
#'
create_tcga_methyl_files=function (file, subset = TRUE, mysites = NULL, out.string = "blah", 
                                   output_folder = "./") 
{
  methyl.pre = fread(file, sep = "\t", header = T, stringsAsFactors = F, 
                     check.names = F, nrows = 2)
  colnames(methyl.pre)[-1] = substr(colnames(methyl.pre)[-1], 1, 15)
  colnames(methyl.pre)[-1] = make.names(colnames(methyl.pre)[-1])
  methyl.colnums = c(1, seq(2, ncol(methyl.pre), 4))
  my_samps.normals = c(1, which(as.numeric(sapply(colnames(methyl.pre), 
                                                  function(x) strsplit(x, "\\.")[[1]][4])) > 9))
  my_samps.cancers = c(1, which(as.numeric(sapply(colnames(methyl.pre), 
                                                  function(x) strsplit(x, "\\.")[[1]][4])) <= 9))
  my_samps.normals=intersect(my_samps.normals,methyl.colnums)
  my_samps.cancers=intersect(my_samps.cancers,methyl.colnums)
  nam = strsplit(strsplit(file, "/")[[1]][length(strsplit(file, 
                                                          "/")[[1]])], "\\.")[[1]][1]
  if (length(my_samps.normals) == 0) {
    methyl.cancers = fread(file, sep = "\t", header = T, stringsAsFactors = F, 
                           check.names = F, select = my_samps.cancers)
    methyl.cancers = methyl.cancers[-1, ]
    methyl.cancers = as.data.frame(methyl.cancers)
    colnames(methyl.cancers)[1] = "Site"
    methyl.cancers = na.omit(methyl.cancers)
    if (subset == T) {
      methyl.cancers = methyl.cancers[methyl.cancers[, 
                                                     1] %in% mysites, ]
    }
    
    write.table(methyl.cancers, paste0(output_folder, nam, 
                                       "_cancers_methyl_", "450K_", out.string, ".txt"), 
                sep = "\t", quote = F, row.names = F)
  }
  else {
    methyl.cancers = fread(file, sep = "\t", header = T, stringsAsFactors = F, 
                           check.names = F, select = my_samps.cancers)
    methyl.cancers = methyl.cancers[-1, ]
    methyl.cancers = as.data.frame(methyl.cancers)
    colnames(methyl.cancers)[1] = "Site"
    methyl.normals = fread(file, sep = "\t", header = T, stringsAsFactors = F, 
                           check.names = F, select = my_samps.normals)
    methyl.normals = methyl.normals[-1, ]
    methyl.normals = as.data.frame(methyl.normals)
    colnames(methyl.normals)[1] = "Site"
    methyl.normals = na.omit(methyl.normals)
    methyl.cancers = na.omit(methyl.cancers)
    if (subset == T) {
      methyl.cancers = methyl.cancers[methyl.cancers[, 
                                                     1] %in% mysites, ]
      methyl.normals = methyl.normals[methyl.normals[, 
                                                     1] %in% mysites, ]
    }
    write.table(methyl.normals, paste0(output_folder, nam, 
                                       "_normals_methyl_", "450K_", out.string, ".txt"), 
                sep = "\t", quote = F, row.names = F)
    write.table(methyl.cancers, paste0(output_folder, nam, 
                                       "_cancers_methyl_", "450K_", out.string, ".txt"), 
                sep = "\t", quote = F, row.names = F)
  }
}