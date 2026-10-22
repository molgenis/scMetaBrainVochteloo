from shared import *

def countUnique(setA, setB):
    unique = 0
    for item in setA:
        if item not in setB:
            unique += 1
    return unique

def getUniqueGenes(setA, setB):
    unique = []
    for gene in setA:
        if gene not in setB:
            unique.append(gene)
    unique.sort()
    return unique


def getSharedGenes(setA, setB):
    overlap = []
    for gene in setA:
        if gene in setB:
            overlap.append(gene)
    overlap.sort()
    return overlap

def addToMap(di, item):
    ctr = di.get(item)
    if ctr is None:
        ctr = 0
    ctr += 1
    di[item] = ctr

def getWithCountBiggerThanOne(di):
    items = set()
    for gene, ct in di.items():
        if ct > 1:
            items.add(gene)
    items = list(items)
    items.sort()
    return items

##############################################################
# data loading section
##############################################################

celltypes = ["AST","END","EX","IN","MIC","OLI","OPC","PER"]
celltypeset = set()
celltypeset.update(celltypes)

gwascolocfile = "/path/to/coloc_results_all_cell_types_all_traits.wfdr.txt"
phenomap = getPhenoMap(gwascolocfile,ds1col="Cell type",ds2col="GWAS",cs1col="eQTL credible set",cs2col="GWAS credible set",existingphenomap=None)

# per phenotype overview
phenos = list(phenomap.keys())
phenos.sort()
gwasphenos = list()
celltypephenos = list()
for pheno in phenos:
    if pheno not in celltypeset:
        gwasphenos.append(pheno)
    else:
        celltypephenos.append(pheno)

origphenos = set()
origphenos.update(phenos)

# gwasphenos = ["AD_EADB_GCST90704646"]
gwasphenoset = set()
gwasphenoset.update(gwasphenos)

bulkcolocfile = "/path/to/coloc_results_all_bulk_tissues_all_cell_types.wfdr.txt"
sccolocfile = "/path/to/coloc_results_all_datasets.zero_cs_removed.wfdr.txt"
phenomap = getPhenoMap(sccolocfile,ds1col="Dataset 1",ds2col="Dataset 2",cs1col="Dataset 1 credible set",cs2col="Dataset 2 credible set",existingphenomap=phenomap)
phenomap = getPhenoMap(bulkcolocfile,ds1col="Cell type",ds2col="Bulk tissue",cs1col="SC credible set",cs2col="Bulk credible set",existingphenomap=phenomap)

bulkphenos = list()
for pheno in phenomap.keys():
    if pheno not in origphenos:
        bulkphenos.append(pheno)
bulkphenos.sort()

# load in locus information (for GWASes)
sumstatdir = "/path/to/gwas_summary_statistics_v2/"
locusbreakerfiles = {
    "AD_EADB_GCST90704646":sumstatdir+"AD_EADB_2026/GCST90704646/3-locus_breaker/AD_GCST90704646.sortedByPosition.bgzip_loci.txt",
    "ADHD_vanDerLaan_GCST90568441":sumstatdir+"ADHD_vanDerLaan_2025/GCST90568441/4-locus_breaker/ADHD_GCST90568441.hg19tohg38.sortedByPosition.bgzip_loci.txt",
    "ALS_vanRheenen_GCST90027163":sumstatdir+"ALS_Rheenen_2021/4-locus_breaker/ALS_GCST90027163.hg19tohg38.sortedByPosition.bgzip_loci.txt",
    "BD_OConnel_eur":sumstatdir+"BD_OConnel_2025/5-locus_breaker/bip2024_eur_no23andMe.hg19tohg38.withBetas.sortedByPosition.bgzip_loci.txt",
    "EA_Okbay":sumstatdir+"EA_Okbay_2022/4-locus_breaker/EA_Okbay_additive_excl_23andMe.hg19tohg38.sortedByPosition.bgzip_loci.txt",
    "Epilepsy_ILAECoCE_TRANS_all_epilepsy":sumstatdir+"Epilepsy_ILAECoCE_2023/ILAE3_TRANS_all_epilepsy_final/6-locus_breaker/ILAE3_TRANS_all_epilepsy_final.hg19tohg38.sortedByPosition.bgzip.PcolRenamed_loci.txt",
    "Intellegence_Savage_GCST006250":sumstatdir+"Intelligence_Savage_2018/4-locus_breaker/Intelligence_GCST006250.hg19tohg38.sortedByPosition.bgzip_loci.txt",
    "MDD_PGC_div":sumstatdir+"MDD_PGC_2025/pgc-mdd2025_no23andMe_div/4-locus_breaker/ILAE3_TRANS_all_epilepsy_final.hg19tohg38.sortedByPosition.bgzip.tsv.gz_loci.txt",
    "MS_Patsopoulos":sumstatdir+"MS_Patsopoulos_2019/4-locus_breaker/Patsopoulos_MS.SEerrorfixed.MRBasecleaned.hg19tohg38.sortedByPosition.bgzip_loci.txt",
    "PD_Hampton_all":sumstatdir+"PD_Hampton_2025/GP2_ALL_EUR_ALL_DATASET/3-locus_breaker/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.sortedByPosition.bgzip_loci.txt",
    "SCZ_Trubetskoy_primary":sumstatdir+"SCZ_Trubetskoy_2022/PGC3_SCZ_wave3.primary/4-locus_breaker/PGC3_SCZ_wave3.primary.autosome.hg19tohg38.sortedByPosition.bgzip_loci.txt"
}

for gwasname, locusbreakerfile in locusbreakerfiles.items():
    print(f"{gwasname}\t{locusbreakerfile}")
    pheno = phenomap.get(gwasname)
    pheno.addLocusBreaker(locusbreakerfile)



for pheno in gwasphenos:
    gwasphenoset = set()
    gwasPhenoObj = phenomap.get(pheno)
    fho = open("./summary/"+pheno+"_summary.txt",'w')
    gwasphenoset.add(pheno)

    ##############################################################
    # overall summary
    ##############################################################

    nrCredibleSetsGWAS = gwasPhenoObj.countCredibleSets()
    nrTestedGenesGWAS = gwasPhenoObj.countGenes(None, False)
    nrTestedGenesGWASSingleVariant = gwasPhenoObj.countGenes(True, False)
    nrTestedGenesGWASMultiVariant = gwasPhenoObj.countGenes(False, False)
    nrSignificantGWASGenes = gwasPhenoObj.countGenes(None, True)
    nrSignificantGWASGenesSingleVariant = gwasPhenoObj.countGenes(True, True)
    nrSignificantGWASGenesMultiVariant = gwasPhenoObj.countGenes(False, True)
    nrSignificantGWASCredibleSets = gwasPhenoObj.countCredibleSets(True)
    
    fho.write(pheno+"\n")
    fho.write(f"{len(gwasPhenoObj.loci)} loci\n")
    fho.write(f"{ nrCredibleSetsGWAS } credible sets\n")
    fho.write(f"{  nrTestedGenesGWAS} tested genes\n")
    fho.write(f"{  nrTestedGenesGWASSingleVariant} tested genes (single variant coloc)\n")
    fho.write(f"{  nrTestedGenesGWASMultiVariant} tested genes (multi variant coloc)\n")

    fho.write(f"{  nrSignificantGWASGenes} significantly colocalizing genes\n")
    fho.write(f"{  nrSignificantGWASGenesSingleVariant} significantly colocalizing genes (single variant coloc)\n")
    fho.write(f"{  nrSignificantGWASGenesMultiVariant} significantly colocalizing genes (multi variant coloc)\n")
    fho.write(f"{  nrSignificantGWASCredibleSets} significantly colocalizing credible sets\n")
           
    fho.write("##############################################################\n"
            "# gwas: summary per cell type\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()


    # Summary per cell type
    header = ("Cell type"
                "\tNumber of coloc tests"
                "\tNumber of coloc tests (single variant coloc)"
                "\tNumber of coloc tests (multi variant coloc)"

                "\tNumber of tested genes"
                "\tNumber of tested genes (single variant coloc)"
                "\tNumber of tested genes (multi variant coloc)"

                "\tNumber of tested credible sets (eQTL)"
                "\tNumber of tested credible sets (GWAS)"

                # "\tNumber of overall coloc genes"
                "\tNumber of significant coloc tests"
                "\tNumber of significant coloc tests (single variant coloc)"
                "\tNumber of significant coloc tests (multi variant coloc)"

                "\tNumber of significant coloc genes"
                "\tNumber of significant coloc genes (single variant coloc)"
                "\tNumber of significant coloc genes (multi variant coloc)"
                
                "\tNumber of significant credible sets (eQTL)"
                "\tNumber of significant credible sets (GWAS)"

                # "\tPercentage of total significantly colocalizing tests"
                # "\tPercentage of total  significantly colocalizing tests (single variant coloc)"
                # "\tPercentage of total  significantly colocalizing tests (multi variant coloc)"

                "\tPercentage of total significantly colocalizing genes"
                "\tPercentage of total significantly colocalizing genes (single variant coloc)"
                "\tPercentage of total significantly colocalizing genes (multi variant coloc)"
                "\tPercentage of significantly colocalizing credible sets (GWAS)"

    )
    print()
    print()
    print(header)
    fho.write(header+"\n")
    for celltype in celltypes:
        phenoObj = phenomap.get(celltype)
        
        testsoverall = phenoObj.countColocTests(useSingleVariantFilter=None, useSignificanceFilter=False, countGenes=False, filterPhenoList=gwasphenoset)
        testssinglevariant = phenoObj.countColocTests(useSingleVariantFilter=True, useSignificanceFilter=False, countGenes=False, filterPhenoList=gwasphenoset)
        testsmultivariant = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=False, countGenes=False, filterPhenoList=gwasphenoset)

        testedgenesoverall = phenoObj.countColocTests(useSingleVariantFilter=None, useSignificanceFilter=False, countGenes=True, filterPhenoList=gwasphenoset)
        testedgenessinglevariant = phenoObj.countColocTests(useSingleVariantFilter=True, useSignificanceFilter=False, countGenes=True, filterPhenoList=gwasphenoset)
        testedgenesmultivariant = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=False, countGenes=True, filterPhenoList=gwasphenoset)
        testedcrediblesets = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=False, countGenes=False, countCredibleSets = True, filterPhenoList=gwasphenoset) # need to fix this

        signtestsoverall = phenoObj.countColocTests(useSingleVariantFilter=None, useSignificanceFilter=True, countGenes=False, filterPhenoList=gwasphenoset)
        signtestssinglevariant = phenoObj.countColocTests(useSingleVariantFilter=True, useSignificanceFilter=True, countGenes=False, filterPhenoList=gwasphenoset)
        signtestsmultivariant = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=True, countGenes=False, filterPhenoList=gwasphenoset)

        signgenesoverall = phenoObj.countColocTests(useSingleVariantFilter=None, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
        signgenessinglevariant = phenoObj.countColocTests(useSingleVariantFilter=True, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
        signgenesmultivariant = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
        signcrediblesets = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=True, countGenes=False, countCredibleSets = True, filterPhenoList=gwasphenoset) # need to fix this
        

        # propA = round((signtestsoverall / testsoverall) ,2) if testsoverall > 0 else 0
        # propB = round((signtestssinglevariant / testssinglevariant) ,2) if testssinglevariant > 0 else 0
        # propC = round((signtestsmultivariant / testsmultivariant) ,2) if testsmultivariant > 0 else 0
        propD = round((signgenesoverall / nrSignificantGWASGenes) * 100 ,2) if nrSignificantGWASGenes > 0 else 0
        propE = round((signgenessinglevariant / nrSignificantGWASGenesSingleVariant) * 100 ,2) if nrSignificantGWASGenesSingleVariant > 0 else 0
        propF = round((signgenesmultivariant / nrSignificantGWASGenesMultiVariant) * 100 ,2) if nrSignificantGWASGenesMultiVariant > 0 else 0
        propH = round((signcrediblesets[1] / nrSignificantGWASCredibleSets) * 100 ,2) if nrSignificantGWASCredibleSets > 0 else 0
        
        lnout = (f"{celltype}\t{testsoverall}\t{testssinglevariant}\t{testsmultivariant}"
                    f"\t{testedgenesoverall}\t{testedgenessinglevariant}\t{testedgenesmultivariant}"
                    f"\t{testedcrediblesets[0]}\t{testedcrediblesets[1]}"
                    f"\t{signtestsoverall}\t{signtestssinglevariant}\t{signtestsmultivariant}"
                    f"\t{signgenesoverall}\t{signgenessinglevariant}\t{signgenesmultivariant}"
                    f"\t{signcrediblesets[0]}\t{signcrediblesets[1]}"
                    f"\t{propD}\t{propE}\t{propF}\t{propH}"
        )
        fho.write(lnout+"\n")
        print(lnout)

    fho.write("\n")
    fho.write("\n")
    fho.write("##############################################################\n"
            "# between cell type comparisons\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()

    # Colocalizing genes per cell type


# -	AST	END	EX	IN	MIC	OLI	OPC	PER	ColocGenes	
#NrUnique	PercentageUniqueWithinCellType	PercentageUniqueOverAllCellTypes	UniqueGenesThatAreEQTLInOtherCellType
    header = ("Celltype")
    for celltypeA in celltypes:
        header += (
            f"\t{celltype} overlapping significantly colocalizing genes"
            f"\t{celltype} overlapping significantly colocalizing genes (single variant coloc)"
            f"\t{celltype} overlapping significantly colocalizing genes (multi variant coloc)"
            f"\t{celltype} overlapping significantly colocalizing GWAS credible sets"
        )
    header += "\tUnique genes\tUnique genes (single variant coloc)\tUnique genes (multi variant coloc)\tUnique GWAS credible sets"
    header += "\tPercentage unique genes\tPercentage unique genes (single variant coloc)\tPercentage unique genes (multi variant coloc)\tPercentage unique GWAS credible sets"
    header += "\tUnique genes that were also tested for coloc in another cell type"

    print(header)
    fho.write(header+"\n")

    # Comparing colocalizing genes between cell types
    sharedgenesSingleVariant = {}
    for celltypeA in celltypes:
        phenoObjA = phenomap.get(celltypeA)
        ctAgenes = set()
        ctAgenesSingleVariant = set()
        ctAgenesMultiVariant = set()
        ctAcs = set()

        for coloc in phenoObjA.colocresults:
            if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                ctAgenes.add(coloc.gene)
                if coloc.singlevariant:
                    ctAgenesSingleVariant.add(coloc.gene)
                else:
                    ctAgenesMultiVariant.add(coloc.gene)
                    ctAcs.add(coloc.getCredibleSetByPhenoObj(gwasPhenoObj))

        ctBgenes = set()
        ctBgenesSingleVariant = set()
        ctBgenesMultivariant = set()
        ctBcs = set()
        ctBeGenes = set()
        outln = celltypeA
        for celltypeB in celltypes:
            phenoObjB = phenomap.get(celltypeB)
            overlapgenes = set()
            overlapgenesSingleVariant = set()
            overlapgenesMultiVariant = set()
            overlapCs = set()
            for coloc in phenoObjB.colocresults:
                if coloc.hasPhenotypeInSet(gwasphenoset) and celltypeA != celltypeB:
                    ctBeGenes.add(coloc.gene)
                if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                    if coloc.gene in ctAgenes:
                        overlapgenes.add(coloc.gene)
                        if celltypeA != celltypeB:
                            ctBgenes.add(coloc.gene)
                        if coloc.singlevariant:
                            overlapgenesSingleVariant.add(coloc.gene)
                            if celltypeA != celltypeB:
                                ctBgenesSingleVariant.add(coloc.gene)
                        else:
                            gwascs = coloc.getCredibleSetByPhenoObj(gwasPhenoObj)
                            if gwascs in ctAcs:
                                overlapgenesMultiVariant.add(coloc.gene)
                                overlapCs.add(gwascs)
                                if celltypeA != celltypeB:
                                    ctBgenesMultivariant.add(coloc.gene)
                                    ctBcs.add(gwascs)
            outln += ( f"\t{len(overlapgenes)}"
                f"\t{len(overlapgenesSingleVariant)}"
                f"\t{len(overlapgenesMultiVariant)}"
                f"\t{len(overlapCs)}"
            )

        uniqueGenes = countUnique(ctAgenes,ctBgenes)
        uniqueGenesSingleVariant = countUnique(ctAgenesSingleVariant,ctBgenesSingleVariant)
        uniqueGenesMultiVariant = countUnique(ctAgenesMultiVariant,ctBgenesMultivariant)
        uniqueGWASCs = countUnique(ctAcs,ctBcs)
        uniqueGenesThatAreEQTLInOtherCellType = set()
        for gene in ctAgenes:
            if gene not in ctBgenes:
                if gene in ctBeGenes:
                    uniqueGenesThatAreEQTLInOtherCellType.add(gene)

        propA = round((uniqueGenes/len(ctAgenes))*100,2) if len(ctAgenes) > 0 else 0
        propB = round((uniqueGenesSingleVariant/len(ctAgenesSingleVariant))*100,2) if len(ctAgenesSingleVariant) > 0 else 0
        propC = round((uniqueGenesMultiVariant/len(ctAgenesMultiVariant))*100,2) if len(ctAgenesMultiVariant) > 0 else 0
        propD = round((uniqueGWASCs/len(ctAcs))*100,2) if len(ctAcs) > 0 else 0

        outln += (f"\t{uniqueGenes}"
                f"\t{uniqueGenesSingleVariant}"
                f"\t{uniqueGenesMultiVariant}"
                f"\t{uniqueGWASCs}"

                f"\t{propA}"
                f"\t{propB}"
                f"\t{propC}"
                f"\t{propD}"
                f"\t{len(uniqueGenesThatAreEQTLInOtherCellType)}"
        )
        fho.write(outln+"\n")
        print(outln)


    # Colocalizing genes per cell type
    fho.write("\n")
    fho.write("\n")
    fho.write("##############################################################\n"
            "# coloc genes per cell type\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()
    
    header ="Colocalizing genes per cell type\nCell type\tGenes\tGenes (single variant coloc)\tGenes (multi variant coloc)"
    print(header)
    fho.write(header+"\n")
    
    for celltypeA in celltypes:
        phenoObjA = phenomap.get(celltypeA)
        ctAgenes = set()
        ctAgenesSingleVariant = set()
        ctAgenesMultiVariant = set()

        for coloc in phenoObjA.colocresults:
            if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                ctAgenes.add(coloc.gene)
                if coloc.singlevariant:
                    ctAgenesSingleVariant.add(coloc.gene)
                else:
                    ctAgenesMultiVariant.add(coloc.gene)

        ctAgenes = list(ctAgenes)
        ctAgenesSingleVariant = list(ctAgenesSingleVariant)
        ctAgenesMultiVariant = list(ctAgenesMultiVariant)
        ctAgenes.sort()
        ctAgenesSingleVariant.sort()
        ctAgenesMultiVariant.sort()

        lnout = celltypeA +"\t"+";".join(ctAgenes)+"\t"+";".join(ctAgenesSingleVariant)+"\t"+";".join(ctAgenesMultiVariant)
        print(lnout)
        fho.write(lnout+"\n")

    # Unique coloc genes per cell type
    
    fho.write("\n")
    fho.write("\n")
    fho.write("##############################################################\n"
            "# unique genes per cell type\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()
    
    header ="Unique colocalizing genes per cell type\nCell type\tGenes\tGenes (single variant coloc)\tGenes (multi variant coloc)"
    print(header)
    fho.write(header+"\n")
    for celltypeA in celltypes:
        phenoObjA = phenomap.get(celltypeA)
        ctAgenes = set()
        ctAgenesSingleVariant = set()
        ctAgenesMultiVariant = set()

        for coloc in phenoObjA.colocresults:
            if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                ctAgenes.add(coloc.gene)
                if coloc.singlevariant:
                    ctAgenesSingleVariant.add(coloc.gene)
                else:
                    ctAgenesMultiVariant.add(coloc.gene)

        ctBgenes = set()
        ctBgenesSingleVariant = set()
        ctBgenesMultiVariant = set()
        for celltypeB in celltypes:
            if celltypeB != celltypeA:
                phenoObjB = phenomap.get(celltypeB)
                for coloc in phenoObjB.colocresults:
                    if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                        ctBgenes.add(coloc.gene)
                        if coloc.singlevariant:
                            ctBgenesSingleVariant.add(coloc.gene)
                        else:
                            ctBgenesMultiVariant.add(coloc.gene)
        print()
        unique = getUniqueGenes(ctAgenes, ctBgenes)
        uniqueSingleVariant = getUniqueGenes(ctAgenesSingleVariant, ctBgenesSingleVariant)
        uniqueMultiVariant = getUniqueGenes(ctAgenesMultiVariant, ctBgenesMultiVariant)
     
        lnout = (celltypeA 
                +"\t"+";".join(unique)
                +"\t"+";".join(uniqueSingleVariant)
                +"\t"+";".join(uniqueMultiVariant)
        )
        print(lnout)
        fho.write(lnout+"\n")

    fho.write("\n")
    fho.write("\n")
    fho.write("##############################################################\n"
            "# shared genes between pairs of cell types\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()

    # Shared genes per cell type
    # Celltype1	Celltype2	NrGenes	Genes
    header = "Cell type 1\tCell type 2\tNumber overlapping genes\tNumber overlapping genes (single variant coloc)\tNumber overlapping genes (multi variant coloc)\tGenes\tGenes (single variant coloc)\tGenes (multi variant coloc)"
    fho.write(header+"\n")
    print(header)


    for i in range(0,len(celltypes)):
        celltypeA = celltypes[i]      
        phenoObjA = phenomap.get(celltypeA)
        ctAgenes = set()
        ctAgenesSingleVariant = set()
        ctAgenesMultiVariant = set()

        for coloc in phenoObjA.colocresults:
            if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                ctAgenes.add(coloc.gene)
                if coloc.singlevariant:
                    ctAgenesSingleVariant.add(coloc.gene)
                else:
                    ctAgenesMultiVariant.add(coloc.gene)

        ctBgenes = set()
        ctBgenesSingleVariant = set()
        ctBgenesMultiVariant = set()
        for celltypeB in celltypes:
            if celltypeB != celltypeA:
                phenoObjB = phenomap.get(celltypeB)
                for coloc in phenoObjB.colocresults:
                    if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                        ctBgenes.add(coloc.gene)
                        if coloc.singlevariant:
                            ctBgenesSingleVariant.add(coloc.gene)
                        else:
                            ctBgenesMultiVariant.add(coloc.gene)
        shared = getSharedGenes(ctAgenes,ctBgenes)
        sharedSingleVariant = getSharedGenes(ctAgenesMultiVariant,ctBgenesSingleVariant)
        sharedMultiVariant = getSharedGenes(ctAgenesMultiVariant,ctBgenesMultiVariant)
        outln = (celltypeA
                +"\tAll other cell types"
                +"\t"+str(len(shared))
                +"\t"+str(len(sharedSingleVariant))
                +"\t"+str(len(sharedMultiVariant))
                +"\t"+";".join(shared)
                +"\t"+";".join(sharedSingleVariant)
                +"\t"+";".join(sharedMultiVariant))
        fho.write(outln+"\n")
        print(outln)

        for j in range(i+1,len(celltypes)):
            celltypeB = celltypes[j]      
            phenoObjB = phenomap.get(celltypeB)
            ctBgenes = set()
            ctBgenesSingleVariant = set()
            ctBgenesMultiVariant = set()
            
            if celltypeB != celltypeA:
                phenoObjB = phenomap.get(celltypeB)
                for coloc in phenoObjB.colocresults:
                    if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                        ctBgenes.add(coloc.gene)
                        if coloc.singlevariant:
                            ctBgenesSingleVariant.add(coloc.gene)
                        else:
                            ctBgenesMultiVariant.add(coloc.gene)

            shared = getSharedGenes(ctAgenes,ctBgenes)
            sharedSingleVariant = getSharedGenes(ctAgenesSingleVariant,ctBgenesSingleVariant)
            sharedMultiVariant = getSharedGenes(ctAgenesMultiVariant,ctBgenesMultiVariant)
            
            outln = celltypeA+"\t"+celltypeB+"\t"+";".join(shared)+"\t"+";".join(sharedSingleVariant)+"\t"+";".join(sharedSingleVariant)
            outln = (celltypeA
                +"\t"+celltypeB
                +"\t"+str(len(shared))
                +"\t"+str(len(sharedSingleVariant))
                +"\t"+str(len(sharedMultiVariant))
                +"\t"+";".join(shared)
                +"\t"+";".join(sharedSingleVariant)
                +"\t"+";".join(sharedMultiVariant))
            fho.write(outln+"\n")
            print(outln)


    fho.write("\n")
    fho.write("\n")
    fho.write("##############################################################\n"
            "# genes colocalizing in multiple cell types\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()

    geneCtr = {}
    geneCtrSingleVariant = {}
    geneCtrMultiVariant = {}
    fho.write("Genes with colocalization in > 1 cell type")
    for celltype in celltypes:
        phenoObjA = phenomap.get(celltype)
        seen = set()
        seenSingleVariant = set()
        seenMultiVariant = set()
        for coloc in phenoObjA.colocresults:
            if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                if coloc.gene not in seen:
                    addToMap(geneCtr,coloc.gene)
                    seen.add(coloc.gene)
                if coloc.singlevariant:
                    if coloc.gene not in seenSingleVariant:
                        addToMap(geneCtrSingleVariant,coloc.gene)
                        seenSingleVariant.add(coloc.gene)
                else:
                    if coloc.gene not in seenMultiVariant:
                        addToMap(geneCtrMultiVariant,coloc.gene)
                        seenMultiVariant.add(coloc.gene)
    print(geneCtr)
    outln = "Overall:\t"+str(len(getWithCountBiggerThanOne(geneCtr)))+"\t"+";".join(getWithCountBiggerThanOne(geneCtr))
    print(outln)
    fho.write(outln+"\n")
    outln = "Single variant coloc:\t"+str(len(getWithCountBiggerThanOne(geneCtrSingleVariant)))+"\t"+";".join(getWithCountBiggerThanOne(geneCtrSingleVariant))
    print(outln)
    fho.write(outln+"\n")
    outln = "Multi variant coloc:\t"+str(len(getWithCountBiggerThanOne(geneCtrMultiVariant)))+"\t"+";".join(getWithCountBiggerThanOne(geneCtrMultiVariant))
    print(outln)
    fho.write(outln+"\n")
    # Genes that colocalize in >1 cell type:

    # Genes per locusbreaker locus

    fho.write("\n")
    fho.write("\n")
    fho.write("##############################################################\n"
            "# genes per locusbreaker locus\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()

    header = ("Colocalizing genes per locus\n"
                "LocusNr"
                "\tChromosome"
                "\tStart basepair position"
                "\tEnd basepair position"
                "\tNumber of tested genes"
                "\tNumber of tested genes (single variant coloc)"
                "\tNumber of tested genes (multi variant coloc)"
                "\tNumber of colocalizing genes"
                "\tNumber of colocalizing genes (single variant coloc)"
                "\tNumber of colocalizing genes (multi variant coloc)"
                "\tColocalizing genes (single variant coloc)"
                "\tColocalizing genes (multi variant coloc)"
                )
    print(header)
    fho.write(header+"\n")
    locusctr = 1
    wcolocctr = 0
    for locus in gwasPhenoObj.loci:
        outln = f"{locusctr}\t{locus.chr}\t{locus.start}\t{locus.stop}"
        locusctr += 1
        testedgenes = set()
        testedgenesSingleVariant = set()
        testedgenesMultiVariant = set()
        genes = set()
        genesSingleVariant = set()
        genesMultivariant = set()
        hascoloc = False
        for coloc in gwasPhenoObj.colocresults:
            if locus.overlaps(coloc.gwasvarchr,coloc.gwasvarpos):
                testedgenes.add(coloc.gene)
                if coloc.significant:
                    genes.add(coloc.gene)
                    hascoloc = True
                if coloc.singlevariant:
                    testedgenesSingleVariant.add(coloc.gene)
                    if coloc.significant:
                        genesSingleVariant.add(coloc.gene)
                else:
                    testedgenesMultiVariant.add(coloc.gene)
                    if coloc.significant:
                        genesMultivariant.add(coloc.gene)
        genes = list(genes)
        genesSingleVariant = list(genesSingleVariant)
        genesMultivariant = list(genesMultivariant)

        svgenestr = ";".join(genesSingleVariant)
        mvgenestr = ";".join(genesMultivariant)

        outln += (f"\t{len(testedgenes)}"
                    f"\t{len(testedgenesSingleVariant)}"
                    f"\t{len(testedgenesMultiVariant)}"
                    f"\t{len(genes)}"
                    f"\t{len(genesSingleVariant)}"
                    f"\t{len(genesMultivariant)}"
                    f"\t{svgenestr}"
                    f"\t{mvgenestr}"
        )
        print(outln)
        fho.write(outln+"\n")
        if hascoloc:
            wcolocctr += 1
    perc = round ( (wcolocctr/len(gwasPhenoObj.loci)) * 100) if len(gwasPhenoObj.loci) > 0 else 0
    outln = f"{wcolocctr} / {len(gwasPhenoObj.loci)} ({perc}%) have at least one colocalizing gene"
    print(outln)
    fho.write(outln+"\n")

    fho.write("\n")
    fho.write("\n")
    fho.write("##############################################################\n"
            "# Comparisons with previous studies\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()

    # Colocalizing in our study but not colocalizing but an eQTL in another publication
    header = "Cell type\tSignificant colocalizing genes\tOverlapping eGenes in Bryois in cell type\tOverlapping eGenes in Fujita in cell type\tOverlapping eGenes in SingleBrain in cell type\tTotal Shared Genes\tUnique genes to this scMetaBrain\tColocalizing genes that are not significant in previous studies\tPercentage"
    print(header)
    fho.write(header+"\n")
    signGenesBryois, eGenesBryois, signGenesFujita, eGenesFujita, signGenesSingleBrain, eGenesSingleBrain = getEGeneListsFromOtherPubs(celltypes)
    for celltype in celltypes:
        phenoObj = phenomap.get(celltype)
        colocGenes = set()
        for coloc in phenoObj.colocresults:
            if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                colocGenes.add(coloc.gene)
        
        otherpubs = set()
        egenesBryoisForCellType = signGenesBryois.get(celltype)
        if egenesBryoisForCellType is None:
            egenesBryoisForCellType = set()
        otherpubs.update(egenesBryoisForCellType)
        egenesFujitaForCellType = signGenesFujita.get(celltype)
        if egenesFujitaForCellType is None:
            egenesFujitaForCellType = set()
        otherpubs.update(egenesFujitaForCellType)
        egenesSingleBrainForCellType = signGenesSingleBrain.get(celltype)
        if egenesSingleBrainForCellType is None:
            egenesSingleBrainForCellType = set()
        otherpubs.update(egenesSingleBrainForCellType)


        overlapBryois = len(getSharedGenes(colocGenes,egenesBryoisForCellType))
        overlapFujita = len(getSharedGenes(colocGenes,egenesFujitaForCellType))
        overlapSingleBrain = len(getSharedGenes(colocGenes,egenesSingleBrainForCellType))
        overallOverlap = len(getSharedGenes(colocGenes,otherpubs))
        unique = getUniqueGenes(colocGenes,otherpubs)
        uniquestr = ";".join(unique)

        Percentage = round((len(unique) / len(colocGenes)) * 100,2) if len(colocGenes) > 0 else 0

        outln = (
                f"{celltype}"
                f"\t{len(colocGenes)}"
                f"\t{overlapBryois}"
                f"\t{overlapFujita}"
                f"\t{overlapSingleBrain}"
                f"\t{overallOverlap}"
                f"\t{len(unique)}"
                f"\t{uniquestr}"
                f"\t{Percentage}"

        )
        print(outln)
        fho.write(outln+"\n")
    

    fho.write("\n")
    fho.write("\n")
    fho.write("##############################################################\n"
            "# Comparisons with bulk\n"
            "##############################################################\n\n")
    fho.write("\n")
    fho.write("\n")
    print()
    print()

    header = ("Celltype"
                "\tNumber of colocalizing genes"
                "\tNumber of colocalizing genes (single variant coloc)"
                "\tNumber of colocalizing genes (multi variant coloc)"
                "\tNumber of colocalizing eQTL credible sets"
    )
    for bulkpheno in bulkphenos:
        header += ( f"\t{bulkpheno} - Significant overlapping genes"
                f"\t{bulkpheno} - Significant overlapping genes (single variant coloc)"
                f"\t{bulkpheno} - Significant overlapping genes (multi variant coloc)"
                f"\t{bulkpheno} - Significant overlapping credible sets (multi variant coloc)"
                )
    header += ("\tPercentage unique genes"
                "\tPercentage unique genes (single variant coloc)"
                "\tPercentage unique genes (multi variant coloc)"
                "\tPercentage unique credible sets"
                "\tGenes not colocalizing in bulk (single variant coloc)"
                "\tGenes not colocalizing in bulk (multi variant coloc)"
    )

    print(header)
    fho.write(header+"\n")
        
    for celltype in celltypes:
        phenoObj = phenomap.get(celltype)
        significantcrediblesets = set()
        genes = set()
        genessinglevariant = set()
        genesmultivariant = set()
        signcrediblesets = set()
        for coloc in phenoObj.colocresults:
            if coloc.significant and coloc.hasPhenotypeInSet(gwasphenoset):
                significantcrediblesets.add(coloc.getCredibleSetByPhenoObj(phenoObj))
                genes.add(coloc.gene)
                if coloc.singlevariant:
                    genessinglevariant.add(coloc.gene)
                else:
                    genesmultivariant.add(coloc.gene)
                    signcrediblesets.add(coloc.getCredibleSetByPhenoObj(phenoObj))

        outln = f"{celltype}\t{len(genes)}\t{len(genessinglevariant)}\t{len(genesmultivariant)}\t{len(signcrediblesets)}"
        celltypeset = set()
        celltypeset.add(celltype)

        alloverlappinggenes = set()
        alloverlappinggenessinglevariant = set()
        alloverlappinggenesmultivariant = set()
        alloverlappingcrediblesets = set()

        for bulkpheno in bulkphenos:
            bulkPhenoObj = phenomap.get(bulkpheno)
            overlapgenes = set()
            overlapgenessinglevariant = set()
            overlapgenesmultivariant = set()
            overlapcrediblesets = set()
            for coloc in bulkPhenoObj.colocresults:
                if coloc.significant and coloc.hasPhenotypeInSet(celltypeset):
                    celltypecs = coloc.getCredibleSetByPhenoObj(phenoObj)
                    if celltypecs in significantcrediblesets:
                        overlapgenes.add(coloc.gene)
                        alloverlappinggenes.add(coloc.gene)
                        if coloc.singlevariant:
                            overlapgenessinglevariant.add(coloc.gene)
                            alloverlappinggenessinglevariant.add(coloc.gene)
                        else:
                            overlapgenesmultivariant.add(coloc.gene) 
                            alloverlappinggenesmultivariant.add(coloc.gene)
                            overlapcrediblesets.add(celltypecs)
                            alloverlappingcrediblesets.add(celltypecs)
          
            outln += (f"\t{len(overlapgenes)}"
                        f"\t{len(overlapgenessinglevariant)}"
                        f"\t{len(overlapgenesmultivariant)}"
                        f"\t{len(overlapcrediblesets)}"
            )

        unique = getUniqueGenes(genes,alloverlappinggenes)
        uniquesinglevariant = getUniqueGenes(genessinglevariant,alloverlappinggenessinglevariant)
        uniquemultivariant = getUniqueGenes(genesmultivariant,alloverlappinggenesmultivariant)
        uniquecs = getUniqueGenes(signcrediblesets,alloverlappingcrediblesets)



        propA = round((len(unique) / len(genes)) * 100,2) if len(genes) > 0 else 0
        propB = round((len(uniquesinglevariant) / len(genessinglevariant)) * 100,2) if len(genessinglevariant) > 0 else 0
        propC = round((len(uniquemultivariant) / len(genesmultivariant)) * 100,2) if len(genesmultivariant) > 0 else 0
        propD = round((len(uniquecs) / len(signcrediblesets)) * 100,2) if len(signcrediblesets) > 0 else 0

        uniqueStr1 = ";".join(uniquesinglevariant)
        uniqueStr2 = ";".join(uniquemultivariant)
        outln += f"\t{propA}\t{propB}\t{propC}\t{propD}\t{uniqueStr1}\t{uniqueStr2}"
        print(outln)
        fho.write(f"{outln}\n")
    fho.close()



    






