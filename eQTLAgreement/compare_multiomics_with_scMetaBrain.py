import gzip
import Rb
import os
import sys
import Shared


def run(ctA, ctB, outdir, genelimitfile):
    print(f"----------- {ctA} vs {ctB} -----------")

    sceqtlgenfolder = "/path/to/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/"
    sceqtlgentopfxfolder = "/path/to/refQTLs/snMultiOme/"
    sceqtlgentopfx = {"B": "B.txt",
                    "DC": "DC.txt",
                    "CD4T": "CD4T.txt",
                    "CD8T": "CD8T.txt",
                    "monocyte": "monocyte.txt",
                    "NK": "NK.txt"}
    sceqtlgenallfx = {"B": "B/qtl_results_all.txt.gz",
                      "DC": "DC/qtl_results_all.txt.gz",
                    "CD4T": "CD4T/qtl_results_all.txt.gz",
                    "CD8T": "CD8T/qtl_results_all.txt.gz",
                    "monocyte": "monocyte/qtl_results_all.txt.gz",
                    "NK": "NK/qtl_results_all.txt.gz"}
                    

    scMetaBrainAlleleFile = "/path/to/allDatasets.fillTags.MAF_HWE_CR_filtered.R2_filtered.unique.samples.bgzip_alleles.txt.gz"
    scMetaBrainDir = "/path/to/output/"
    scMetaBrainQtlFileTpl = "chrCHR_qtl_results_all.txt.gz"
    scMetaBrainSignQtlFileTpl = "top_qtl_results_all_MultiTest.txt"

    ctASignPairs = None

    ctBisMetaBrain = False
    topfxeqtlgen = sceqtlgentopfx.get(ctA)
    # fopen(sceqtlgenfolder+topfxeqtlgen,'r')
    if topfxeqtlgen is not None and os.path.exists(sceqtlgentopfxfolder+topfxeqtlgen):
        ctAGenes, ctASignPairs = Shared.getLimixEffects(sceqtlgentopfxfolder+topfxeqtlgen, True, None, genelimit, variantEncoding = "GTEX", allelesInSnpId=True)
        ctBisMetaBrain = True
    elif os.path.exists(scMetaBrainDir+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA)):
        tmpgenes, ctASignPairs = Shared.getLimixEffects(scMetaBrainDir+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA), True, None, genelimit, variantEncoding="GTEX", alleleFile = scMetaBrainAlleleFile)
    else:
        print("Could not find file for ctA: "+ctA)
        print(os.path.exists(sceqtlgentopfxfolder+topfxeqtlgen))
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
            tmpgenes, overlapChr = Shared.getLimixEffects(qtlfile, False, ctASignPairs, genelimit, variantEncoding = "GTEX", allelesInSnpId=True)
            overlapCtB.update(overlapChr)
            returned += 1
            print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')
    else:
        allfxfile = sceqtlgenallfx.get(ctB)
        tmpgenes, overlapCtB = Shared.getLimixEffects(sceqtlgenfolder+allfxfile, False, ctASignPairs, genelimit, variantEncoding = "GTEX", allelesInSnpId=True)
   
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