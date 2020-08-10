# Jeff Jones
# SoCal Bioinformatics Inc. 2019
# 
# plotting helper functions

# set background to white
theme_set(theme_bw())

# ggplot short-cuts
axis_x_rotate <- function(x=90){
  hjust <- 1
  if(x == 0)
    hjust <- 0.5
  
  theme(axis.text.x = element_text(angle = x, hjust = hjust))
}

theme_noleg <- theme(legend.position = "none")

theme_palette <- c("#44709d", "#d97828", "#83992a", "#995d81", 
                   "#a63d40", "#8ea4d2", "#deb340", "#704e2e",  
                   "#69d1c5", "#f7b267", "#496f5d", "#e4572e", 
                   "#48beff", "#8b8c89", "#e56399", "#5c2751"
)

theme_blank_x <- theme(
  axis.text.x=element_blank(), axis.ticks.x=element_blank())

theme_blank_y <- theme(
  axis.text.y=element_blank(), axis.ticks.y=element_blank())

theme_color <- scale_color_manual(values=theme_palette)
theme_fill <- scale_fill_manual(values=theme_palette)

# save a plot
save_plot <- function(image, pth = "../plots/new_plot", type = "pdf", w = 6, h = 7){
  
  pth <- paste(pth, ".", type, sep = "")
  
  if(type == "pdf"){
    pdf(file = pth, pointsize = 8, width = w, height = h)
    plot(image)
    dev.off()
  } else {
    h <- 640 * h/6
    w <- 640 * w/6
    png(filename = pth, pointsize = 10, width = w, height = h)
    plot(image)
    dev.off()
  }
  echo <- paste("saving image ... ", pth, "\n", sep="")
  cat(echo)
}