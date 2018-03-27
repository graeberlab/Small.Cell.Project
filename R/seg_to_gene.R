#' does survival analysis
#' 
#' Reads file with features(PC's) in columns and samples in rows. And an annotation file.You pick a max # of features. Using caret
#' runs a prediction iteratively from 2 to max.features. Picks best predictor with least featurs. Returns samples correctly predicted
#' leaves out those samples incorrectly predicted.
#' 
#' @param seg_file  segment file
#' @param gene_file  gene to segment mapping
#' @param output_name name of output file will be postended with "_by_gene.txt"
#' @importFrom CNTools getRS CNSeg rs
#'
#' @export

#library(CNTools)


#data("sampleData")
#myGeneInfo <-read.table('mod_RefSeq_GrCh37.txt', header = TRUE)

seg_to_gene <- function(seg_file,gene_file,output_name){
options(scipen=999)
	seg <-read.table(seg_file, header = TRUE, sep='\t')         
	colnames(seg) <- colnames(sampleData)
	print(head(seg))
	cnseg <- CNSeg(seg)
	myGeneInfo <-read.table(gene_file, header = TRUE)
	rdByGene <- getRS(cnseg, by = "gene", imput = FALSE, XY = FALSE, geneMap = myGeneInfo, what = "mean")
	readByGene <- rs(rdByGene)
	return(readByGene)
	write.table(readByGene, file = paste0(output_name,"_by_gene.txt"), sep = "\t", row.names = FALSE, quote=FALSE)

}

# setwd("C:/Users/Nick/Desktop/ImmunoActiva")
# gene_reads <- map.to.gene("CCLE_copynumber_2013-12-03_lung.seg")
# rounded <- round(gene_reads[,-1:-5], 3)
# rounded_gene_reads <- cbind(gene_reads[,1:5], rounded)
# write.table(rounded_gene_reads, file = "CCLE_lung_CNA_gene.txt", sep = "\t", row.names = FALSE, quote=FALSE)

