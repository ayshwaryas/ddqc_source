import pandas as pd

annotations = pd.read_csv("/Volumes/scqc/data/misc/pbmc_seurat/annotations.csv")
cells = pd.read_csv("/Users/michaelalperovich/Dropbox/Proga/Primes/primes_storage/output_pg/misc/pbmc_seurat/1.4-mad-2/!cells.csv")
annotations.index = annotations.barcodekey
cells.index = cells.barcodekey

unique_cutoff = annotations.loc[annotations.index.difference(cells.index)]
unique_cutoff.annotations.value_counts()

unique_ddqc = cells.loc[cells.index.difference(annotations.index)]
unique_ddqc.cluster_labels.value_counts()


# ddqc_filtered_counts.loc[ddqc_filtered_counts.index.intersection(unique_cutoff.index)]