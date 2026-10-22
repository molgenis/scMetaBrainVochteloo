import gzip
import Rb
import os
import sys
import Shared

    
# this analysis is performed in one direction; from scMetaBrain to CIMA
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

    bryoisdir = "/path/to/single_cell_eQTLsumstats/Bryois2022/"
    bryoisfile = "CELLTYPE.CHROMOSOME.gz"
    bryoisallelefile = "/path/to/single_cell_eQTLsumstats/Bryois2022/snp_pos.txt.gz"
    bryoistopfxfile = "CELLTYPE.top.gz"

    scMetaBrainDir = "/path/to//Bryois/LIMIX/"
    scMetaBrainQtlFileTpl = "chrCHR_qtl_results_all_sortedbysnp.sortfix.bgz.txt.gz"
    scMetaBrainSignQtlFileTpl = "top_qtl_results_all_acatfix.txt"
    scMetaBrainAlleleFile = "/path/to//allDatasets.fillTags.MAF_HWE_CR_filtered.R2_filtered.unique.samples.bgzip_alleles.txt.gz"
    os.makedirs(outdir,exist_ok=True)

    ctASignPairs = {}
    overlapCtB = {}

    if ctAisScMetaBrain:
        signfile = scMetaBrainDir+"/"+ctA+"/"+scMetaBrainSignQtlFileTpl.replace("CELLTYPE",ctA)
        ctASignGenes, ctASignPairs = Shared.getLimixEffects(signfile, True, None, genelimit, alleleFile = scMetaBrainAlleleFile, geneEncoding="ENSG", recalculateP=True)
        for chr in range(1,23):
            overlapCtBtmp = Shared.findQTLsInBryois(bryoisdir+bryoisfile.replace("CELLTYPE",ctB).replace("CHROMOSOME",f"{chr}"), ctASignPairs)
            overlapCtB.update(overlapCtBtmp)
            print(f"{len(overlapCtB)} / {len(ctASignPairs)} pairs found so far")

        # add allele info
        rsidToData = {}
        for pair in overlapCtB.keys():
            data = overlapCtB[pair]
            rsid = data["rsid"]
            rsiddata = rsidToData.get(rsid)
            if rsiddata is None:
                rsiddata = []
            rsiddata.append(data)
            rsidToData[rsid] = rsiddata 

        
        replaced = 0
        fh = Shared.fopen(bryoisallelefile,'r')
        for line in fh:
            if not line.startswith("#"):
                elems = line.strip().split("\t")
                rsid = elems[0]
                rsiddata = rsidToData.get(rsid)
                if rsiddata is not None:
                    alleles = elems[4]+"/"+elems[3]
                    effectallele = elems[3]
                    for dc in rsiddata:
                        # check if eff allele is a1 or a2...
                        dc["alleles"] = alleles
                        dc["effectallele"] = effectallele
                        replaced += 1
        fh.close()
        print(f"{replaced} / {len(overlapCtB)} have replaced alleles")
    else:

        ctASignPairs = Shared.findQTLsInBryoisTopFx(bryoisdir+bryoistopfxfile.replace("CELLTYPE",ctA))
        
        overlapCtB = {}
        for chr in range(1,23):
            ctBfile = scMetaBrainDir + ctB + "/"+scMetaBrainQtlFileTpl.replace("CHR",str(chr))
            ctBSignGenes, ctBtmpPairs = Shared.getLimixEffects(ctBfile, False, ctASignPairs, genelimit, alleleFile = scMetaBrainAlleleFile, geneEncoding="ENSG", recalculateP=True)
            overlapCtB.update(ctBtmpPairs)
    
    print(f"{len(overlapCtB)} / {len(ctASignPairs)} pairs found total")

    

    # compare
    print("Done loading files. Making comparisons.")
    theta = 0.0
    if ctAisScMetaBrain:
        ctA = "scMbBryois"+ctA
        ctB = "Bryois"+ctB
    else:
        ctB = "scMbBryois"+ctB
        ctA = "Bryois"+ctA
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