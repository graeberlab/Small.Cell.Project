#' Plot PLS
#' 
#' Plots PLS from scores file (output of PLSR_from_file)
#' 
#' @param file File containing scores matrix
#' @param info.name Vector of sample names
#' @param info.type Vector of sample types in the same order
#' @param title Title of the plot
#' @param labels default=T
#' @param PCx,PCy PCs to display
#' @param ellipse Construct confidence region based on groups in info.type, default = T
#' @param conf default = 0.95 
#' @param saveplot do you want to save the plot
#' @param savetype default=".pdf", ".png" is also possible
#' @param w  width default=8
#' @param h  height default=6
#' @param h  colpalette can add your own vector of colors to color groups. default is ggplot default colors
#' @param h  legendname legend name
#' 
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' 
#' @export
#' 
plot_pls<-function (file, info.name, info.type, title = "", labels = TRUE, 
          PCx = "comp.1", PCy = "comp.2", ellipse = F, conf = 0.95, 
          saveplot = F, savename = "default", savetype = ".pdf", w = 8, 
          h = 6,legendname="type",colpalette=NULL) 
{
  require(ggplot2)
  require(vegan)
  table <- read.table(file, header = TRUE)
  table$type = info.type[match(table$Score, info.name)]
  exp_var = read.delim(paste0(gsub("scores.txt", "", file), 
                              "pve.txt"), row.names = 1)
  exp_var$pve = round(exp_var[, 1] * 100, digits = 2)
  pcx.y <- ggplot(table, aes_string(x = PCx, y = PCy)) + geom_point(size = I(2),aes(color = factor(type))) + theme(legend.position = "right", 
  plot.title = element_text(size = 30), legend.text = element_text(size = 22),legend.title = element_text(size = 20), axis.title = element_text(size = 30), 
  legend.background = element_rect(), axis.text.x = element_text(margin = margin(b = -2)),axis.text.y = element_text(margin = margin(l = -14))) + 
    guides(color = guide_legend(title = "Type")) + labs(title = title,x = paste0(PCx, " (", exp_var$pve[match(PCx, rownames(exp_var))],"%)"), y = paste0(PCy, " (", exp_var$pve[match(PCy, 
    rownames(exp_var))], "%)")) + theme_bw() + 
    if (labels == TRUE) {
      geom_text(data = table, mapping = aes(label = Score), check_overlap = TRUE, size = 3)
    }
  if (!is.null(colpalette)) {
    pcx.y <- pcx.y + scale_color_manual(legendname, 
                                              values = colpalette)
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
      pcx.y2
    }
    else {
      pcx.y2
    }
  }
  else {
    if (saveplot == T) {
      ggsave(paste0(savename, "_", PCx, "_vs_", PCy, savetype), 
             dpi = 300, plot = pcx.y, width = w, height = h)
      pcx.y
    }
    else {
      pcx.y
    }
  }
}
