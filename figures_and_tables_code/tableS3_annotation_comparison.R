old.dir.prefix <- "/Volumes/scqc/output_copies/output_pg04-20-21/"
dir.prefix <- "/Volumes/scqc/output_pg/"
organism <- c("human", "mouse")[2]
method <- c("1.4-none-0", "1.4-joint_clustering_old")[2]



fout <- paste0("~/Downloads/compare_annotations_", organism, "_", method, ".csv")
fout_all <- paste0("~/Downloads/compare_annotations_table.csv")
write("project,tissue,method,cluster,annotation,new.cell.type,new.markers", fout)

cnt <- 0
cnt.unknown <- 0
cnt.50 <- 0
cnt.75 <- 0
cnt.incorrect <- 0



for (project in c("tabula_senis_30m")){#list.dirs(dir.prefix, full.names = FALSE, recursive = FALSE)) {
  for (tissue in list.dirs(paste0(dir.prefix, project), full.names = FALSE, recursive = FALSE)) {
    if (organism == "human") {
      if (project == "human_other") {
        if (tissue != "adipose" && tissue != "Heart_Circulation" && tissue != "krasnow_lung") {
          next 
        }
      } else if (project != "human_tissue_atlas") {
        next
      }
    } else {
      if (project != "tabula_muris" && project != "tabula_muris_smartseq2" && project != "tabula_senis_24m" && project != "tabula_senis_30m") {
        next
      }
    }
    
  
    #old.output.dir <- paste0(old.dir.prefix, project, "/", tissue, "/", method, "/")
    output.dir <- paste0(dir.prefix, project, "/", tissue, "/", method, "/")
    #old.clusters <- read.csv(paste0(old.output.dir, "!clusters.csv"))
    #rownames(old.clusters) <- old.clusters$cluster
    clusters <- read.csv(paste0(output.dir, "!clusters.csv"))
    rownames(clusters) <- clusters$cluster
    
    for (cl in rownames(clusters)) {
      cnt <- cnt + 1
      mrk <- strsplit(as.character(clusters[cl,]$markers), split=";")[[1]]
      
      ann <- tolower(as.character(clusters[cl,]$annotation)) 
      ct <- tolower(as.character(clusters[cl,]$cell_type))
      
      if (ann == "unknown") {
        cnt.unknown <- cnt.unknown + 1
        next
      }
      
      if (as.numeric(clusters[cl,]$X.annotation) > 0.50 && ann != "unknown") {
        cnt.50 <- cnt.50 + 1
      }
      
      if (as.numeric(clusters[cl,]$X.annotation) > 0.75 && ann != "unknown") {
        cnt.75 <- cnt.75 + 1
      }
       
      if (as.numeric(clusters[cl,]$X.annotation) < 0.75) {
        next
      }
      
      
      if (ann != ct && ann != "unknown") {
        if (ann == substring(ct, 1, nchar(ct) - 1)) {
          next
        }
        
        
        if (ann == "endocardial cell" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "cardiac muscle cell" && ct == "cardiomyocytes") {
          next
        }
        
        if (ann == "kidney tubule cell" && ct == "distal tubule cells") {
          next
        }
        
        if (ann == "kidney tubule cell" && ct == "proximal tubule cells") {
          next
        }
        
        if (ann == "t cell" && ct == "t memory cells") {
          next
        }
        
        if (ann == "stromal cell" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "natural killer cell" && ct == "nk cells") {
          next
        }
        
        if (ann == "stromal cell" && ct == "pericytes") {
          next
        }
        
        if (ann == "luminal cell of lactiferous duct" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "keratinocyte" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "basal cell of epidermis" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "basal cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "chondroblast" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "mesenchymal stem cell" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "mesenchymal cell" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "type ii pneumocyte" && ct == "pulmonary alveolar type ii cells") {
          next
        }
        
        if (ann == "fraction a pre-pro b cell" && ct == "b cells") {
          next
        }
        
        if (ann == "t cell" && ct == "gamma delta t cells") {
          next
        }
        
        if (ann == "mesenchymal stem cell of adipose" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "mature natural killer cell" && ct == "nk cells") {
          next
        }
        
        if (ann == "naive b cell" && ct == "b cells") {
          next
        }
        
        if (ann == "pro-b cell" && ct == "b cells") {
          next
        }
        
        if (ann == "immature b cell" && ct == "b cells") {
          next
        }
        
        if (ann == "late pro-b cell" && ct == "b cells") {
          next
        }
        
        if (ann == "hematopoietic precursor cell" && ct == "hematopoietic stem cells") {
          next
        }
        
        if (ann == "immature t cell" && ct == "t cells") {
          next
        }
        
        if (ann == "immature t cell" && ct == "t memory cells") {
          next
        }
        
        if (ann == "immature t cell" && ct == "gamma delta t cells") {
          next
        }
        
        if (ann == "keratinocyte stem cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "pancreatic acinar cell" && ct == "acinar cells") {
          next
        }
        
        if (ann == "lung endothelial cell" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "endothelial cell of hepatic sinusoid" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "skeletal muscle satellite cell" && ct == "satellite cells") {
          next
        }
        
        if (ann == "epithelial cell of proximal tubule" && ct == "proximal tubule cells") {
          next
        }
        
        if (ann == "microglial cell" && ct == "microglia") {
          next
        }
        
        if (ann == "oligodendrocyte precursor cell" && ct == "oligodendrocyte progenitor cells") {
          next
        }
        
        if (ann == "astrocyte of the cerebral cortex" && ct == "astrocytes") {
          next
        }
        
        if (ann == "epithelial cell of large intestine" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "enterocyte of epithelium of large intestine" && ct == "enterocytes") {
          next
        }
        
        if (ann == "large intestine goblet cell" && ct == "goblet cells") {
          next
        }
        
        if (ann == "skeletal muscle satellite stem cell" && ct == "satellite cells") {
          next
        }
        
        if (ann == "uminal epithelial cell of mammary gland" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "fibroblast of cardiac tissue" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "granulocyte" && ct == "neutrophils") {
          next
        }
        
        if (ann == "keratinocyte stem cell" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "luminal epithelial cell of mammary gland" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "fenestrated cell" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "macrophage dendritic cell progenitor" && ct == "macrophages") {
          next
        }
        
        if (ann == "megakaryocyte-erythroid progenitor cell" && ct == "erythroid-like and erythroid precursor cells") {
          next
        }
        
        if (ann == "leukocyte" && ct == "macrophages") {
          next
        }
        
        if (ann == "b cell (plasmocyte)" && ct == "plasma cells") {
          next
        }
        
        if (ann == "endothelial cell (apc)" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "stratified epithelial cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "erythroid cell" && ct == "erythroid-like and erythroid precursor cells") {
          next
        }
        
        if (ann == "b cell (plasmocyte)" && ct == "b cells") {
          next
        }
        
        if (ann == "pancreas exocrine cell" && ct == "acinar cells") {
          next
        }
        
        if (ann == "ec" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "pericyte/ec" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "fenestrated_endothelial" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "myeloid(mono+/dc/mac)" && ct == "dendritic cells") {
          next
        }
        
        if (ann == "smooth muscle cell" && ct == "myocytes") {
          next
        }
        
        if (ann == "smooth muscle cell" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "enterocyte" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "epithelial cell" && ct == "proximal tubule cells") {
          next
        }
        
        if (ann == "loop of henle" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "at2 cell" && ct == "pulmonary alveolar type ii cells") {
          next
        }
        
        if (ann == "endothelial cell (endothelial to mesenchymal transition)" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "natural killer" && ct == "nk cells") {
          next
        }
        
        if (ann == "alveolar epithelial type 2" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "alveolar epithelial type 1" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "b" && ct == "b cells") {
          next
        }
        
        if (ann == "thyroid follicular cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "pre-adipocyte" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "adipocyte_progenitor(diff)" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "capillary" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "capillary aerocyte" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "artery" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "lymphatic" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "atrial myocyte" && ct == "cardiomyocytes") {
          next
        }
        
        if (ann == "bladder cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "bladder cell" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "bladder urothelial cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "brain pericyte" && ct == "smooth muscle cells") {
          next
        }
        
        if (ann == "brush cell of epithelium proper of large intestine" && ct == "tuft cells") {
          next
        }
        
        if (ann == "cardiac neuron" && ct == "schwann cells") {
          next
        }
        
        if (ann == "dn1 thymic pro-t cell" && ct == "nk cells") {
          next
        }
        
        if (ann == "dn4 thymocyte" && ct == "t memory cells") {
          next
        }
        
        if (ann == "endothelial cell of coronary artery" && ct == "endothelial cells") {
          next
        }
        
        if (ann == "enterocyte of epithelium of large intestine" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "epidermal cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "epithelial cell of large intestine" && ct == "enterocytes") {
          next
        }
        
        if (ann == "epithelial cell of lung" && ct == "pulmonary alveolar type ii cells") {
          next
        }
        
        if (ann == "hematopoietic stem cell" && ct == "erythroid-like and erythroid precursor cells") {
          next
        }
        
        if (ann == "kidney collecting duct cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "kupffer cell" && ct == "macrophages") {
          next
        }
        
        if (ann == "large intestine goblet cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "leukocyte" && ct == "t cells") {
          next
        }
        
        if (ann == "lymphocyte" && ct == "b cells") {
          next
        }
        
        if (ann == "mesenchymal cell" && ct == "smooth muscle cells") {
          next
        }
        
        if (ann == "monocyte" && ct == "macrophages") {
          next
        }
        
        if (ann == "myeloid cell" && ct == "macrophages") {
          next
        }
        
        if (ann == "myeloid cell" && ct == "neutrophils") {
          next
        }
        
        if (ann == "myeloid cell" && ct == "dendritic cells") {
          next
        }
        
        if (ann == "naôøωôøωve b cell" && ct == "b cells") {
          next
        }
        
        if (ann == "pancreatic a cell" && ct == "alpha cells") {
          next
        }
        
        if (ann == "pancreatic ductal cell" && ct == "epithelial cells") {
          next
        }
        
        if (ann == "pancreatic stellate cell" && ct == "fibroblasts") {
          next
        }
        
        if (ann == "proerythroblast" && ct == "erythroid-like and erythroid precursor cells") {
          next
        }
        
        if (ann == "proerythroblast + erythroblast" && ct == "erythroid-like and erythroid precursor cells") {
          next
        }
        
        if (ann == "professional antigen presenting cell" && ct == "dendritic cells") {
          next
        }
        
        if (ann == "professional antigen presenting cell" && ct == "b cells") {
          next
        }
        
        if (ann == "skeletal muscle satellite stem cell" && ct == "myoblasts") {
          next
        }
        
        if (ann == "t cell" && ct == "nk cells") {
          next
        }
        
        if (ann == "type b pancreatic cell" && ct == "beta cells") {
          next
        }
        

        
        if (ann == "" && ct == "") {
          next
        }
        
        
        
        if (project == "human_tissue_atlas" && tissue == "Rectum" && cl == 13 && method == "1.4-none-0") {
          next
        }
        
        if (project == "human_tissue_atlas" && tissue == "Muscle" && cl == 13 && method == "1.4-none-0") {
          next
        }
        
        if (project == "human_tissue_atlas" && tissue == "Pleura" && cl == 15 && method == "1.4-joint_clustering_old") {
          next
        }        
        
        if (project == "human_tissue_atlas" && tissue == "Pancreas" && cl == 2 && method == "1.4-joint_clustering_old") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Trachea" && cl == 10 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Trachea" && cl == 14 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Trachea" && cl == 17 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Trachea" && cl == 20 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Adipose" && cl == 25 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Fat" && cl == 20 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Trachea" && cl == 21 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris" && tissue == "Heart_and_Aorta" && cl == 17 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris" && tissue == "Kidney" && cl == 26 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris" && tissue == "Marrow" && cl == 26 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Kidney" && cl == 6 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Kidney" && cl == 8 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_senis_24m" && tissue == "Kidney" && cl == 14 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris" && tissue == "Lung" && cl == 42 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Bone_Marrow" && cl == 16 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Heart_and_Aorta" && cl == 12 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_muris_smartseq2" && tissue == "Liver" && cl == 12 && method == "1.4-none-0") {
          next
        }
        
        if (project == "tabula_senis_24m" && tissue == "Spleen" && cl == 17 && method == "1.4-none-0") {
          next
        }
        
        
        
        write(paste(project, tissue, method, cl, as.character(clusters[cl,]$annotation), 
                    as.character(clusters[cl,]$cell_type),
                    paste(mrk[1:min(25, length(mrk))], collapse=";", sep=""), 
                    sep=","), file = fout, append = TRUE)
        cnt.incorrect <- cnt.incorrect + 1
      }
    }
  }
}
write(paste(organism, project, method, cnt, cnt.50, cnt.75, cnt.incorrect, round((cnt.75 - cnt.incorrect) / cnt.75, 3), sep = ","), file=fout_all, append = TRUE)
