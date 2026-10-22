import gzip
import Rb
import os
import sys
import Shared
from py_qvalue import qvalue

def run(ctA, ctB, outdir, useFixedFDR, useFlippedZ, genelimitfile):
    print(f"----------- {ctA} vs {ctB} -----------")
    singleBrainDir = "/path/to/annotation/single_cell_eQTLsumstats/SingleBrain2026/data/"
    singleBrainQtlFileTpl = "TISSUE_eqtl_full_assoc.tsv.gz"
    singleBrainSignQtlFileTpl = "TISSUE_eqtl_top_assoc.tsv.gz-withFixedRandomFxZScoreAndFixedFxQval.tsv.gz"

    scMetaBrainAlleleFile = "/path/to/allDatasets.fillTags.MAF_HWE_CR_filtered.R2_filtered.unique.samples.bgzip_alleles.txt.gz"


    scMetaBrainDir = "/path/to/SingleBrain/LIMIX/"
    scMetaBrainQtlFileTpl = "chrCHR_qtl_results_all_sortedbysnp.bgz.txt.gz"
    scMetaBrainSignQtlFileTpl = "top_qtl_results_all_acatfix.txt"
       
    os.makedirs(outdir,exist_ok=True)

    genelimit = None
    if genelimitfile is not None:
        genelimit = set()
        print("Parsing: "+genelimitfile)
        fh = Shared.fopen(genelimitfile,'r')
        for line in fh:
            elems = line.strip().split(".")
            genelimit.add(elems[0])
            # print(elems[0])
        fh.close()
        print(f"{len(genelimit)} genes to limit to")

    ctAMetaBrainFile = scMetaBrainDir+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA)
    ctAsingleBrainFile = singleBrainDir+singleBrainSignQtlFileTpl.replace("TISSUE",ctA)
        
    ctASignPairs = None
    overlapCtB = {}
    returned = 0
    if os.path.exists(ctAMetaBrainFile):
        print(ctAMetaBrainFile)
        ctASignGenes, ctASignPairs = Shared.getLimixEffects(ctAMetaBrainFile, True, None, genelimit, geneEncoding="ENSG", alleleFile = scMetaBrainAlleleFile)
        
        qtlfile = singleBrainDir + singleBrainQtlFileTpl.replace("TISSUE",ctB)
        overlapChr = Shared.findQTLsInSingleBrain(qtlfile, ctASignPairs, useFixedFDR, genelimit)
        overlapCtB.update(overlapChr)
        returned += 1
        print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')
    else:
        print(ctAsingleBrainFile)
        ctASignPairs = Shared.findQTLsInSingleBrain(ctAsingleBrainFile, None, useFixedFDR, genelimit)
            
        for chr in range(1,23):
            chrstr = str(chr)
            qtlfile = scMetaBrainDir +"/"+ctB+"/"+scMetaBrainQtlFileTpl.replace("CELLTYPE",ctB).replace("CHR",chrstr)
            tmpgenes, overlapChr = Shared.getLimixEffects(qtlfile, False, ctASignPairs, genelimit, geneEncoding="ENSG", alleleFile = scMetaBrainAlleleFile, recalculateP=True)
            overlapCtB.update(overlapChr)
            returned += 1
            print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')

    # compare
    
    print("Done loading files. Making comparisons.")
    theta = 0
    Shared.makeComparisonsSingleBrain(outdir, ctA, ctB, ctASignPairs, overlapCtB, theta, onlyCompareEffectAlleles=False, fixedFxMeta=useFixedFDR, useFlippedRandomZ=useFlippedZ)

ctA = sys.argv[1]
ctB = sys.argv[2]
outdir = sys.argv[3]

useFixedFDR = False
if sys.argv[4] == "YES":
    useFixedFDR = True
useFlippedZ = False
if sys.argv[5] == "YES":
    useFlippedZ = True

genelimit = None
if len(sys.argv) > 6:
    genelimit = sys.argv[6]
run(ctA,ctB,outdir, useFixedFDR, useFlippedZ, genelimit)