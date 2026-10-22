import gzip
import Rb
import os
import sys
import Shared

def run(ctA, ctB, outdir, genelimitfile, ctAisScMetaBrain):
    print(f"----------- {ctA} vs {ctB} -----------")
    genelimit = None
    if genelimitfile is not None:
        genelimit = set()
        fh = Shared.fopen(genelimitfile,'r')
        for line in fh:
            elems = line.strip().split("\t")
            genelimit.add(elems[0])
        fh.close()

    
    
    fujitadir = "/path/to/single_cell_eQTLsumstats/Fujita2022/AnalysisQTL/celltype/"
    fujitafile = "celltype-eqtl-sumstats.CELLTYPE.tsv.gz"
    fujitatopfxfile = "celltype-eqtl-sumstats.CELLTYPEtopfx.txt.gz"


    scMetaBrainDir = "/path/to/Fujita/LIMIX/"
    scMetaBrainQtlFileTpl = "chrCHR_qtl_results_all_sortedbysnp.sortfix.bgz.txt.gz"
    scMetaBrainSignQtlFileTpl = "top_qtl_results_all_acatfix.txt"
    scMetaBrainAlleleFile = "/path/to/allDatasets.fillTags.MAF_HWE_CR_filtered.R2_filtered.unique.samples.bgzip_alleles.txt.gz"
    os.makedirs(outdir,exist_ok=True)

    ctASignPairs = {}
    overlapCtB = {}

    if ctAisScMetaBrain:
        signfile = scMetaBrainDir+"/"+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA)
        ctASignGenes, ctASignPairs = Shared.getLimixEffects(signfile, True, None, genelimit, alleleFile = scMetaBrainAlleleFile, recalculateP=True)
        overlapCtB = Shared.findQTLsInFujita(fujitadir+fujitafile.replace("CELLTYPE",ctB), ctASignPairs)
    else:
        ctASignPairs = Shared.findQTLsInFujita(fujitadir+fujitatopfxfile.replace("CELLTYPE",ctA), None)
        overlapCtB = {}
        for chr in range(1,23):
            ctBfile = scMetaBrainDir + ctB + "/"+scMetaBrainQtlFileTpl.replace("CHR",str(chr))
            ctBSignGenes, ctBtmpPairs = Shared.getLimixEffects(ctBfile, False, ctASignPairs, genelimit, alleleFile = scMetaBrainAlleleFile, recalculateP=True)
            overlapCtB.update(ctBtmpPairs)
        
    # compare
    print("Done loading files. Making comparisons.")
    theta = 1.0

    if ctAisScMetaBrain:
        ctA = "scMbFujita"+ctA
        ctB = "Fujita"+ctB
    else:
        ctB = "scMbFujita"+ctB
        ctA = "Fujita"+ctA

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