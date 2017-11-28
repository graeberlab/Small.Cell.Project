#' Plot PCA
#' 
#' Plots PCA from scores file (output of PCA_from_file)
#' 
#' @param file File containing scores matrix
#' @param info.name Vector of sample names
#' @param info.type Vector of sample types in the same order
#' @param title Title of the plot
#' @param labels default=T, do you want labels
#' @param label_file Separate file that has annotations you want to use for labeling. Samples should be in a column called "Score". 
#' @param label_name Name of column you want to use to label points in plot 
#' @param PCx,PCy PCs to display
#' @param ellipses Construct confidence region based on groups in info.type, default = T
#' @param conf default = 0.95 
#' @param point_size size of points
#' @param saveplot do you want to save the plot
#' @param savetype default=".pdf", ".png" is also possible
#' @param w  width default=8
#' @param h  height default=6
#' @param colpalette  a vector of colors that changes the default color scheme, default is NULL, e.g. something like c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#' @param shape.palette  a vector of shapes that changes the default shapes, default is null,  e.g. something like c(16,16,15,17,17,15,17)

# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
# @importFrom ggrepel geom_text_repel
#' 
#' @export
#' 
plot_pca = function (file, info.name, info.type, title = "", labels = TRUE, 
          label_file = NULL, label_name, PCx = "PC1", PCy = "PC2", 
          ellipses = F, conf = 0.95, point_size = 3,savename="example",saveplot=F,savetype=".png",w=8,h=6,label_size=3,colpalette=NULL,shape.palette=NULL,legendname = "default") {
  require(ggplot2)
  require(vegan)
  table <- read.table(file, header = TRUE)
  table$type = info.type[match(table$Score, info.name)]
  if (!is.null(label_file)) {
    label_frame = read.delim(label_file)
    table = left_join(table, label_frame, by = "Score")
  }
  sdev = read.delim(paste0(gsub("scores.txt", "", file), "sdev.txt"))
  sdev$var = sdev^2
  sdev$pve = round(sdev$var/sum(sdev$var) * 100, digits = 2)
  rownames(sdev) = paste0("PC", seq(1, nrow(sdev)))
  pcx.y <- ggplot(table, aes_string(x = PCx, y = PCy)) + geom_point(size = I(point_size), 
                                                                    aes(color = factor(type))) + theme(legend.position = "right", 
                                                                                                       plot.title = element_text(size = 30), legend.text = element_text(size = 22), 
                                                                                                       legend.title = element_text(size = 20), axis.title = element_text(size = 30), 
                                                                                                       legend.background = element_rect(), axis.text.x = element_text(margin = margin(b = -2)), 
                                                                                                       axis.text.y = element_text(margin = margin(l = -14))) + 
    guides(color = guide_legend(title = "Type")) + labs(title = title, 
                                                        x = paste0(PCx, " (", sdev$pve[match(PCx, rownames(sdev))], 
                                                                   "%)"), y = paste0(PCy, " (", sdev$pve[match(PCy, 
                                                                                                               rownames(sdev))], "%)")) + theme_bw() + if (labels == 
                                                                                                                                                           TRUE && is.null(label_file)) {
                                                                                                                 geom_text_repel(data = table, mapping = aes(label = Score), 
                                                                                                                           check_overlap = TRUE, size = label_size)
                                                                                                               }
  else if (labels == TRUE && !is.null(label_file)) {
    geom_text(data = table, mapping = aes_string(label = label_name), 
              check_overlap = TRUE, size = label_size)
  }
  if(!is.null(shape.palette)){
    pcx.y<- pcx.y + scale_shape_manual(legendname,values=shape.palette)
  }
  if(!is.null(colpalette)){
    pcx.y<- pcx.y + scale_color_manual(legendname,values=colpalette)
  }
  if (ellipses == TRUE) {
    pcx.y<- pcx.y + stat_ellipse(aes(color=factor(type)),level=conf)
    if (saveplot == T) {
      ggsave(paste0(savename, "_", PCx, "_vs_", PCy, savetype), 
             dpi = 300, plot = pcx.y, width = w, height = h)
    pcx.y} else {pcx.y}
  }else {
    if (saveplot == T) {
      ggsave(paste0(savename, "_", PCx, "_vs_", PCy, savetype), 
             dpi = 300, plot = pcx.y, width = w, height = h)
    pcx.y } else { pcx.y }
  }
}
