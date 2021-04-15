source("figures_code/config.R")


fout <- paste0(PATH, "percent_comparison.csv")
write("project,tissue,n_cells after cutoff,n_cells after ddqc,percent difference", fout)

data <- read.csv(paste0(PATH, "figure3_table.csv"))
data$Directory <- gsub("/ahg/regevdata/projects/scqc/output_pg/", "", data$Directory)
data$project <- sapply(strsplit(data$Directory, "/"), "[[", 1)
data$tissue <- sapply(strsplit(data$Directory, "/"), "[[", 2)


for (row in 1:nrow(data)) {
  n.cells.ddqc <- data[row, "n_cells.after.ddqc"]
  n.cells.co <- data[row, "n_cells.after.cutoff"]
  pd <- round((n.cells.ddqc - n.cells.co) / n.cells.co, 3)
  write(paste(data[row, "project"], data[row, "tissue"], n.cells.co,  n.cells.ddqc, pd, sep=","), fout, append = TRUE)
}