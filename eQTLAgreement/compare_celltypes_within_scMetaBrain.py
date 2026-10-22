import gzip
import Rb
import os
import sys
import Shared

def run(ctA,ctB, outdir, genelimitfile):
    print(f"----------- {ctA} vs {ctB} -----------")

    genelimit = None
    if genelimitfile is not None:
        genelimit = set()
        fh = Shared.fopen(genelimitfile,'r')
        for line in fh:
            elems = line.strip().split("\t")
            genelimit.add(elems[0])
        fh.close()

    indir = "/path/to/output/"
    qtlfiletpl = "chrCHR_qtl_results_all.txt.gz"
    signqtlfiletpl = "top_qtl_results_all_MultiTest.txt"
    
    os.makedirs(outdir,exist_ok=True)

    signfile = indir+"/"+ctA+"/"+signqtlfiletpl.replace("CELLTYPE",ctA)
    ctASignGenes, ctASignPairs = Shared.getLimixEffects(signfile, True, None, genelimit)
            
    print(f"{ctA} has {len(ctASignPairs)} significant effects.")

    # compare    
    overlapCtB = {}
    testedGenesCTb = set()
    signfileB = indir+"/"+ctB+"/"+signqtlfiletpl.replace("CELLTYPE",ctB)
    ctBTestedGenes, ctBTopPairs = Shared.getLimixEffects(signfileB, False, None, genelimit)
    
    returned = 0
    for chr in range(1,23):
        chrstr = str(chr)
        qtlfile = indir +"/"+ctB+"/"+qtlfiletpl.replace("CELLTYPE",ctB).replace("CHR",chrstr)
        ctBGenes, ctBSignPairs = Shared.getLimixEffects(qtlfile, False, ctASignPairs.keys(), ctBTestedGenes)

        # overlapChr = findQTLsInFile(qtlfile, chrstr, ctASignGenes, ctASignPairs)
        overlapCtB.update(ctBSignPairs)
        returned += 1
        print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\r')
    
    
    print(f"{returned} processes returned - {len(overlapCtB)} / {len(ctASignPairs)} matches retrieved so far",end='\n')
        
        # compare
    print("Done loading files. Making comparisons.")
    theta = 1.0
    Shared.makeComparisons(outdir, ctA, ctB, ctASignPairs, overlapCtB, theta, onlyCompareEffectAlleles=True)


# celltypes = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
# celltypes.sort()
# print(celltypes)

ctA = sys.argv[1]
ctB = sys.argv[2]
outdir = sys.argv[3]
genelimit = None
if len(sys.argv) > 4:
    genelimit = sys.argv[4]
run(ctA,ctB,outdir, genelimit)
