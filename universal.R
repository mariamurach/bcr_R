createDir <- function(libdir){
  if (!dir.exists(libdir)){
    dir.create(libdir, recursive = T)
  }else{
    print("dir exists")
  }}
setRLib <- function(){
  rver <- getRversion()
  libdir <- paste0("/home/mm5jy/R/", rver)
  createDir(libdir)
}

theme_mm <- function(text_sizes = 8) {

  theme_bw(base_size = text_sizes) %+replace%    #replace elements we want to change
    theme(
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(margin = margin(rep(2,4))),
      axis.line = element_blank(),
      legend.key = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      ),
      strip.background = element_rect(fill = "white"),
      complete = TRUE,
      axis.title.x = element_blank(),
      
      legend.position="right"
    )
}

theme_mm_lines  <- function(text_sizes = 8) {

  theme_bw(base_size = text_sizes) %+replace%    #replace elements we want to change
    theme(
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.grid = element_line(colour = "grey95"),
      plot.title = element_text(margin = margin(rep(2,4))),
      axis.line = element_blank(),
      legend.key = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      ),
      complete = TRUE,
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "white",
                                      colour = "grey20"), 
      legend.position="right"
    )
}

safe_colorblind_palette <- c("#88CCEE",  "#DDCC77", "#CC6677","#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", 
                             "#E0E678", "#E3B476",  "#FFED61", "#8B8BD9", "#C1F5ED",  "#E6C3F7",
                             "#BF1515")

vivid_colors <- c('#4363d8', '#f58231','#e6194b', '#3cb44b', '#ffe119',  '#911eb4',
          '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
          '#9a6324', '#fffac8', '#800000', '#aaffc3',
          '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')