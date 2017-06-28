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
#' @param ellipse Construct confidence region based on groups in info.type, default = T
#' @param conf default = 0.95 
#' @param point_size size of points
#' @param saveplot do you want to save the plot
#' @param savetype default=".pdf", ".png" is also possible
#' @param w  width default=8
#' @param h  height default=6
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' 
#' @export
#' 
plot_pca = function (file, info.name, info.type, title = "", labels = TRUE, 
          label_file = NULL, label_name, PCx = "PC1", PCy = "PC2", 
          ellipse = F, conf = 0.95, point_size = 3,savename="example",saveplot=F,savetype=".pdf",w=8,h=6) {
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
                                                                                                                 geom_text(data = table, mapping = aes(label = Score), 
                                                                                                                           check_overlap = TRUE, size = 3)
                                                                                                               }
  else if (labels == TRUE && !is.null(label_file)) {
    geom_text(data = table, mapping = aes_string(label = label_name), 
              check_overlap = TRUE, size = 3)
  }
  if (ellipse == TRUE) {
    plot(table[, c(PCx, PCy)], main = title)
    ord = ordiellipse(table[, c(PCx, PCy)], table$type, kind = "sd", 
                      conf = conf)
    cov_ellipse <- function(cov, center = c(0, 0), scale = 1, 
                            npoints = 100) {
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }
    df_ell <- data.frame(matrix(ncol = 0, nrow = 0))
    for (g in (droplevels(as.factor(table$type)))) {
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(table[table$type == 
                                                               g, ], cov_ellipse(ord[[g]]$cov, ord[[g]]$center, 
                                                                                 ord[[g]]$scale))), type = g))
    }
    pcx.y2 = pcx.y + geom_path(data = df_ell, aes(x = df_ell[, 
                                                             PCx], y = df_ell[, PCy], colour = type), size = 1, 
                               linetype = 1)
    if (saveplot == T) {
      ggsave(paste0(savename, "_", PCx, "_vs_", PCy, savetype), 
             dpi = 300, plot = pcx.y2, width = w, height = h)
    pcx.y2} else {pcx.y2}
  }else {
    if (saveplot == T) {
      ggsave(paste0(savename, "_", PCx, "_vs_", PCy, savetype), 
             dpi = 300, plot = pcx.y, width = w, height = h)
    pcx.y } else { pcx.y }
  }
}
