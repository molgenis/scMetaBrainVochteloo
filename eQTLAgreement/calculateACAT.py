import sys
import Shared
from py_qvalue import qvalue
import numpy as np


if len(sys.argv) < 2:
    print("Usage: indir")
    exit()

indir = sys.argv[1]
if not indir.endswith("/"):
    indir += "/"
file = indir+"top_qtl_results_all_MultiTest.txt"
print(file)
fh = Shared.fopen(file,'r')
genes = set()
topfx = {}

header = Shared.index( fh.readline().rstrip().split("\t") )
for line in fh:
    elems = line.strip().split("\t")
    snp = elems[header["snp_id"]]
    gene = elems[header["ensembl_id"]]
    genes.add(gene)
    topfx[gene] = snp
fh.close()

print(f"{len(topfx)} top effects.")

geneacatp={}
fho = open(indir+"acatpvals.txt",'w')
fho.write("gene\tacatp\n")
for chr in range(1,23):

    pvalspergene = {}
    file = indir+f"chr{chr}_qtl_results_all_sortedbysnp.bgz.txt.gz"
    print(file)
    fh = Shared.fopen(file,'r')
    header = Shared.index( fh.readline().rstrip().split("\t") )
    pcol = header["p_value"]
    snpcol = header["snp_id"]
    genecol = header["ensembl_id"]
    lctr = 0
    for line in fh:
        elems = line.strip().split("\t")
        snp = elems[snpcol]
        gene = elems[genecol]
        if gene in genes:
            pfloat = float(elems[pcol])
            gpvals = pvalspergene.get(gene)
            if gpvals is None:
                gpvals = []
            gpvals.append(pfloat)
            pvalspergene[gene] = gpvals
        lctr += 1
        if lctr % 100000 == 0:
            print(f"{lctr} - {len(pvalspergene)} genes",end='\r')
    print(f"{lctr} - {len(pvalspergene)} genes",end='\n')
    fh.close()

    print(f"Calculating ACAT for {len(pvalspergene)} genes")
    for gene in pvalspergene:
        data = pvalspergene.get(gene)
        acatp = Shared.acat_test(data)
        geneacatp[gene] = acatp
        fho.write(f"{gene}\t{acatp}\n")
    print(f"{len(geneacatp)} acat pvalues total")
fho.close()

# perform qvalue correction
genes = []
pvals = []
for gene in geneacatp.keys():
    genes.append(gene)
    pvals.append(geneacatp.get(gene))

qvalue_results = qvalue(np.array(pvals))
q_values = qvalue_results['qvalues']

geneidx = Shared.index(genes)

fh = Shared.fopen(indir+"top_qtl_results_all_MultiTest.txt",'r')
fho = Shared.fopen(indir+"top_qtl_results_all_MultiTest_wACATp.txt",'w')
header = fh.readline().rstrip()

fho.write(header+"\tACATp\tACATp-qvalue\n")
header = Shared.index( header.split("\t") )
for line in fh:
    elems = line.strip().split("\t")
    snp = elems[header["snp_id"]]
    gene = elems[header["ensembl_id"]]
    gidx = geneidx[gene]
    qv = q_values[gidx]
    ap = pvals[gidx]
    outln = line.rstrip()+f"\t{ap}\t{qv}\n"
    fho.write(outln)
fh.close()
fho.close()