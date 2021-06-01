genes_orig = open("genes_orig.tsv", "r")
genes_decorated = open("genes_dec_mgi.tsv", "r")

genes_dict = dict()
for line in genes_decorated.readlines():
    line = line.strip().split()
    if len(line) < 2:
        genes_dict[line[0]] = line[0]
    else:
        genes_dict[line[0]] = line[1]

with open("genes.tsv", "w") as fout:
    for line in genes_orig.readlines():
        gene = line.strip().split()[0]
        fout.write("{}\t{}\n".format(gene, genes_dict.get(gene, gene)))
        