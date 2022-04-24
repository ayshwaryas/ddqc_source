library(ggplot2)
library(cowplot)
library(dplyr)

no_bkg <- theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank()) 

ggsave1 <- function(filename, plot, n.clusters=30, type="h") {  # custom ggsave function
  if (type == "h") {
    height = 10
    width = 14 / 30 * max(n.clusters, 30)
  }
  if (type == "v") {
    height = 14 / 30 * max(n.clusters, 30)
    width = 14
  }
  if (type == "u") {
    height = 10
    width = 10 + 2 * ceiling(n.clusters / 13)
  }
  
  ggsave(filename = filename, plot = plot + no_bkg, width = width, height = height) #saves plot with custom dimensions 
}

# path <- "Z:\\revisions\\revisions2\\adaptive_mad_threshold\\"
path <- "C:\\Users\\misha\\Downloads\\adaptive_mad_threshold\\"
tissue <- "heart"

df <- read.csv(paste0(path, tissue, "_filtered_cells_count.csv"), row.names = 1)
plot <- ggplot(df, aes(threshold, cells_filtered)) + geom_line()
plot <- plot + facet_wrap(~cluster, ncol = 4, scales="free_y")
ggsave1(paste0(path, tissue, "_threshold_vs_cells.pdf"), plot)

df_sum <- df %>% group_by(threshold) %>% summarise(cells_filtered = sum(cells_filtered))
plot2 <- ggplot(df_sum, aes(threshold, cells_filtered)) + geom_line()
ggsave1(paste0(path, tissue, "_threshold_vs_cells_sum.pdf"), plot2)

plot3 <- ggplot(df, aes(threshold, cells_filtered_pct)) + geom_line()
plot3 <- plot3 + facet_wrap(~cluster, ncol = 4, scales="free_y")
ggsave1(paste0(path, tissue, "_threshold_vs_cells_pct.pdf"), plot3)
