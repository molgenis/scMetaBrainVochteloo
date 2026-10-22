import gzip
import Rb
import os
import sys
import Shared

# compare scMetaBrain wo DiverseCohorts to DiverseCohorts
def run(ctA, ctB, outdir, genelimitfile, ctAisDiverseCohorts):
    print(f"----------- {ctA} vs {ctB} -----------")
    genelimit = None
    if genelimitfile is not None:
        genelimit = set()
        fh = Shared.fopen(genelimitfile,'r')
        for line in fh:
            elems = line.strip().split("\t")
            genelimit.add(elems[0])
        fh.close()

    diverseCohortsDir = "/path/to/DiverseCohorts/LIMIX/"
    diverseCohortsQtlFileTpl = "chrCHR_qtl_results_all_sortedbysnp.bgz.txt.gz"
    diverseCohortsSignQtlFileTpl = "top_qtl_results_all_MultiTest.txt"

    woDiverseCohortsDir = "/path/to/DiverseCohortsExcluded/LIMIX/"
    woDiverseCohortsQtlFileTpl = "chrCHR_qtl_results_all_sortedbysnp.bgz.txt.gz"
    woDiverseCohortsSignQtlFileTpl = "top_qtl_results_all_MultiTest.txt"

    scMetaBrainAlleleFile = "/path/to/allDatasets.fillTags.MAF_HWE_CR_filtered.R2_filtered.unique.samples.bgzip_alleles.txt.gz"
    os.makedirs(outdir,exist_ok=True)

    ctASignPairs = {}
    overlapCtB = {}

    if ctAisDiverseCohorts:
        signfile = diverseCohortsDir+"/"+ctA+"/"+diverseCohortsSignQtlFileTpl
        ctASignGenes, ctASignPairs = Shared.getLimixEffects(signfile, True, None, genelimit, alleleFile = scMetaBrainAlleleFile, geneEncoding="ENSG", recalculateP=False)
        for chr in range(1,23):
            otherfile = woDiverseCohortsDir+"/"+ctB+"/"+woDiverseCohortsQtlFileTpl.replace("CHR",ctB)
            ctBSignGenes, overlapCtBtmp = Shared.getLimixEffects(otherfile, True, None, genelimit, alleleFile = scMetaBrainAlleleFile, geneEncoding="ENSG", recalculateP=False)
            overlapCtB.update(overlapCtBtmp)
            print(f"{len(overlapCtB)} / {len(ctASignPairs)} pairs found so far")

    else:

        signfile = woDiverseCohortsDir+"/"+ctA+"/"+woDiverseCohortsSignQtlFileTpl
        ctASignGenes, ctASignPairs = Shared.getLimixEffects(signfile, True, None, genelimit, alleleFile = scMetaBrainAlleleFile, geneEncoding="ENSG", recalculateP=False)
        for chr in range(1,23):
            otherfile = diverseCohortsDir+"/"+ctB+"/"+diverseCohortsQtlFileTpl.replace("CHR",ctB)
            ctBSignGenes, overlapCtBtmp = Shared.getLimixEffects(otherfile, True, None, genelimit, alleleFile = scMetaBrainAlleleFile, geneEncoding="ENSG", recalculateP=False)
            overlapCtB.update(overlapCtBtmp)
            print(f"{len(overlapCtB)} / {len(ctASignPairs)} pairs found so far")
    
    print(f"{len(overlapCtB)} / {len(ctASignPairs)} pairs found total")

    

    # compare
    print("Done loading files. Making comparisons.")
    theta = 0.0
    if ctAisDiverseCohorts:
        ctA = "DiverseCohorts-"+ctA
        ctB = "scMbWoDiverseCohorts-"+ctB
    else:
        ctB = "scMbWoDiverseCohorts"+ctB
        ctA = "DiverseCohorts-"+ctA
    Shared.makeComparisons(outdir, ctA, ctB, ctASignPairs, overlapCtB, theta, onlyCompareEffectAlleles=False)


ctA = sys.argv[1]
ctB = sys.argv[2]
outdir = sys.argv[3]
genelimit = None
ctaismb = False
if len(sys.argv) > 4:
    genelimit = sys.argv[4]
if len(sys.argv) > 5:
    data = sys.argv[5]
    if data == "TRUE":
        ctaismb = True

run(ctA,ctB,outdir, None,ctaismb)