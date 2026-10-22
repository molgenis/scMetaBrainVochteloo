import gzip
import Rb
import os
import sys
import Shared

    
# this analysis is performed in one direction; from scMetaBrain to CIMA
def run(ctA, ctB, outdir, genelimitfile):
    print(f"----------- {ctA} vs {ctB} -----------")
    genelimit = None
    if genelimitfile is not None:
        genelimit = set()
        fh = Shared.fopen(genelimitfile,'r')
        for line in fh:
            elems = line.strip().split("\t")
            genelimit.add(elems[0])
        fh.close()

    cimadir = "/path/to/2026-CIMA-OMIX009786/cis_eQTL_Summary_Statics/cis_eQTL/"
    cimafile = "CELLTYPE.txt.gz"

    scMetaBrainDir = "/path/to/output/"
    scMetaBrainQtlFileTpl = "chrCHR_qtl_results_all.txt.gz"
    scMetaBrainSignQtlFileTpl = "top_qtl_results_all_MultiTest.txt"
    scMetaBrainAlleleFile = "/path/to/allDatasets.fillTags.MAF_HWE_CR_filtered.R2_filtered.unique.samples.bgzip_alleles.txt.gz"
    os.makedirs(outdir,exist_ok=True)

    signfile = scMetaBrainDir+"/"+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA)
    ctASignGenes, ctASignPairs = Shared.getLimixEffects(signfile, True, None, genelimit, variantEncoding="CIMA", alleleFile = scMetaBrainAlleleFile)
    
    overlapCtB = Shared.findQTLsInscCIMAFile(cimadir+cimafile.replace("CELLTYPE",ctB), ctASignPairs)

    # compare
    print("Done loading files. Making comparisons.")
    theta = 0
    Shared.makeComparisons(outdir, ctA, ctB, ctASignPairs, overlapCtB, theta, onlyCompareEffectAlleles=False)


ctA = sys.argv[1]
ctB = sys.argv[2]
outdir = sys.argv[3]
genelimit = None
if len(sys.argv) > 4:
    genelimit = sys.argv[4]
run(ctA,ctB,outdir, genelimit)