import gzip
import Rb
import os
import sys
import Shared

def run(ctA, ctB, outdir, genelimitfile):
    print(f"----------- {ctA} vs {ctB} -----------")
    gtexDir = "/path/to/GTEx_Analysis_v10_eQTL_all_associations/"
    gtexQtlFileTpl = "TISSUE.v10.allpairs.chrCHR.parquet"
    gtexSignQtlFileTpl = "TISSUE.v10.eGenes.txt.gz"

    scMetaBrainAlleleFile = "/path/to/allDatasets.fillTags.MAF_HWE_CR_filtered.R2_filtered.unique.samples.bgzip_alleles.txt.gz"
    scMetaBrainDir = "/path/to/output/"
    scMetaBrainQtlFileTpl = "chrCHR_qtl_results_all.txt.gz"
    scMetaBrainSignQtlFileTpl = "top_qtl_results_all_MultiTest.txt"
       
    os.makedirs(outdir,exist_ok=True)

    ctAMetaBrainFile = scMetaBrainDir+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA)
    ctAGTExFile = gtexDir+gtexSignQtlFileTpl.replace("TISSUE",ctA)
    
    
    ctASignPairs = None
    overlapCtB = {}
    returned = 0
    if os.path.exists(ctAMetaBrainFile):
        print(ctAMetaBrainFile)
        ctASignGenes, ctASignPairs = Shared.getLimixEffects(ctAMetaBrainFile, True, None, genelimit, variantEncoding="GTEX", geneEncoding="ENSG", alleleFile = scMetaBrainAlleleFile)
        for chr in range(1,23):
            chrstr = str(chr)
            qtlfile = gtexDir + gtexQtlFileTpl.replace("TISSUE",ctB).replace("CHR",chrstr)
            overlapChr = Shared.findQTLsInGTExParquet(qtlfile, ctASignPairs)
            overlapCtB.update(overlapChr)
            returned += 1
            print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')
    else:
        print(ctAGTExFile)
        ctASignPairs = Shared.getGTExTopFx(ctAGTExFile)
        for chr in range(1,23):
            chrstr = str(chr)
            qtlfile = scMetaBrainDir +"/"+ctB+"/"+scMetaBrainQtlFileTpl.replace("CELLTYPE",ctB).replace("CHR",chrstr)
            tmpgenes, overlapChr = Shared.getLimixEffects(qtlfile, False, ctASignPairs, genelimit, variantEncoding="GTEX", geneEncoding="ENSG", alleleFile = scMetaBrainAlleleFile)
            overlapCtB.update(overlapChr)
            returned += 1
            print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')

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