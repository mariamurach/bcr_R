vis_genus <- function(x, title="",adj =c(-0.25, 0.5), facing = "clockwise" , font_size = 0.8,
                      clrs = vivid_colors,
                      ...) {
  b1a_summa <- x %>% group_by( V, J) %>% summarize(counts = sum(counts))
  b1a_ko <- b1a_summa %>% filter(V != ".", J != ".")
  mat <- matrix(nrow = length(unique(b1a_ko$J)), ncol = length(unique(b1a_ko$V)))
  colnames(mat) <- unique(b1a_ko$V)
  rownames(mat) <- unique(b1a_ko$J)
  ord <-  union( rownames(mat), colnames(mat))
  genes <- c(unique(b1a_ko$J), unique(b1a_ko$V))
  ord <- genes[order(match(genes,levels(genes)))]
  genes <-  genes[order(match(genes,levels(genes)))]
  clrs <- vivid_colors[match(genes, names(vivid_colors))] 

  for (j in b1a_ko$J) {
    for (v in b1a_ko$V) {
      sub <- b1a_ko %>% filter(V == v, J == j)
      if(nrow(sub) == 0)
      {mat[j,v] = 0}else{
        mat[j,v] <- sub$counts}
    }
  }
  # convert to proportions / fraction of total
  mat <- mat/sum(mat)
  
  # reorder by V and J decreasing order (left to right)
  vmax <- colSums(mat)
  jmax <- rowSums(mat) 
  mat2 <- mat[order(jmax, decreasing = FALSE), order(vmax, decreasing = TRUE), drop = F]
  zeros <- sum(mat2 == 0)
  # use circlize
  grid.col  <- clrs[1:length(ord)]
  #names(grid.col) <- ord
  chordDiagram(mat2, annotationTrack = c("grid"), big.gap = 3, small.gap = 1.5, 
               order = ord, grid.col = clrs)
  
  # add legends
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, 
                CELL_META$ylim[1], 
                CELL_META$sector.index, 
                facing = facing, 
                niceFacing = T,
                adj = adj, 
                cex = font_size)
  }, bg.border = NA)
  
  # add title
  #title(title, cex.main = 0.1)
}

do_circos <- function(filename, x,lims,  ...){
  pdf( filename, w = 4, h = 4)
  circos.par("canvas.xlim" = lims, "canvas.ylim" = lims)
    vis_genus(x, ...)
  dev.off()
  circos.clear()
}
