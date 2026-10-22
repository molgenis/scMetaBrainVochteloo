import gzip
import Rb
import os
import sys
import Shared


def run(ctA, ctB, outdir, genelimitfile):
    print(f"----------- {ctA} vs {ctB} -----------")
    metaBrainDir = "/path/to/bulkMetaBrain/input/eqtls/"
    metaBrainQtlFileTpl = "chrCHR-AllEffects-sorted.bgz.txt.gz"
    metaBrainSignQtlFileTpl = "merged-withqval.txt.gz"
    
    scMetaBrainAlleleFile = "/path/to/allDatasets.fillTags.MAF_HWE_CR_filtered.R2_filtered.unique.samples.bgzip_alleles.txt.gz"
    scMetaBrainDir = "/path/to/output/"
    scMetaBrainQtlFileTpl = "chrCHR_qtl_results_all.txt.gz"
    scMetaBrainSignQtlFileTpl = "top_qtl_results_all_MultiTest.txt"

    ctASignPairs = None
    overlapCtB = {}

    ctBisMetaBrain = False

    # fopen(sceqtlgenfolder+topfxeqtlgen,'r')
    topfxmetabrain = metaBrainDir + ctA + "/"+metaBrainSignQtlFileTpl
    if topfxmetabrain is not None and os.path.exists(topfxmetabrain):
        # ctASignPairs = findQTLsInscMetaBrainFile(topfxmetabrain,None, None)

        ctASignPairs = Shared.getEffectsMbQTL(topfxmetabrain,None)
        ctBisMetaBrain = True
    elif os.path.exists(scMetaBrainDir+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA)):
        tmpgenes, ctASignPairs = Shared.getLimixEffects(scMetaBrainDir+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA), True, None, genelimit, alleleFile = scMetaBrainAlleleFile)
    else:
        print("Could not find file for ctA: "+ctA)
        sys.exit()

    print(f"{ctA} has {len(ctASignPairs)} significant effects.")
    print(f"{ctB} is scMetaBrain: {ctBisMetaBrain}")

    
    overlapCtB = {}
    returned = 0
    if ctBisMetaBrain:    
        for chr in range(1,23):
            chrstr = str(chr)
            qtlfile = scMetaBrainDir +"/"+ctB+"/"+scMetaBrainQtlFileTpl.replace("CELLTYPE",ctB).replace("CHR",chrstr)
            # overlapChr = findQTLsInscMetaBrainFile(qtlfile, chrstr, ctASignPairs)
            tmpgenes, overlapChr = Shared.getLimixEffects(qtlfile, False, ctASignPairs, genelimit, alleleFile = scMetaBrainAlleleFile)
            overlapCtB.update(overlapChr)
            returned += 1
            print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')
    else:
        for chr in range(1,23):
            chrstr = str(chr)
            qtlfile = metaBrainDir +"/"+ctB+"/"+metaBrainQtlFileTpl.replace("CHR",chrstr)
            # overlapChr = findQTLsInscMetaBrainFile(qtlfile, chrstr, ctASignPairs)
            overlapChr = Shared.getEffectsMbQTL(qtlfile, ctASignPairs)
            overlapCtB.update(overlapChr)
            returned += 1
            print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')
   
    print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')
 
    os.makedirs(outdir,exist_ok=True)

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