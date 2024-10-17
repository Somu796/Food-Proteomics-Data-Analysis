changingVennFont <- function(p, font){
  
  grep_grob <- function(gt, lab){
    which(sapply(gt, function(x) grepl(lab, x$name)))
  }
  
  p2 <- ggplot_gtable(ggplot_build(p))
  mygrobs <- p2$grobs
  panel_grob <- mygrobs[[grep_grob(mygrobs, "panel")]]
  venn_grob <- panel_grob$children[[grep_grob(panel_grob$children, "venn")]]
  text_grob <- venn_grob$children[grep_grob(venn_grob$children, "text")]
  text_grob <- do.call(grid::gList, 
                       lapply(text_grob, 
                              function(x) {x$gp$fontfamily <- font; 
                              x}))
  venn_grob$children[grep_grob(venn_grob$children, "text")] <- text_grob
  panel_grob$children[[grep_grob(panel_grob$children, "venn")]] <- venn_grob
  mygrobs[[grep_grob(mygrobs, "panel")]] <- panel_grob
  p2$grobs <- mygrobs
  grid::grid.newpage()
  grid::grid.draw(p2)
  
}