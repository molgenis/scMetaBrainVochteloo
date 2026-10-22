from shared import *


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


# ##############################################################
# # gwas overview
# ##############################################################



# # 
# print()
# print()
# print()
fho = open('summary.txt','w')
fho.write("##############################################################\n"
          "# gwas overview\n"
          "##############################################################\n\n")
header = ("Phenotype"
        
        "\tNumber of locusbreaker loci"
        "\tNumber of locusbreaker loci with at least one colocalizing gene"
        "\tNumber of locusbreaker loci with at least one colocalizing gene (single variant coloc)"
        "\tNumber of locusbreaker loci with at least one colocalizing gene (multi variant coloc)"
        "\tPercentage of locusbreaker loci with at least one colocalizing gene"
        "\tPercentage of locusbreaker loci with at least one colocalizing gene (single variant coloc)"
        "\tPercentage of locusbreaker loci with at least one colocalizing gene (multi variant coloc)"
        
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

        )

fho.write(header+"\n")
print(header)

# per GWAS overview
testedGenesOverall = set()
colocGenesOverall = set()
testedGenesSingleVariant = set()
colocGenesSingleVariant = set()
testedGenesMultiVariant = set()
colocGenesMultiVariant = set()

celltypegenetraitcombos = set()
celltypegenetraitcombosSingleVariant = set()
celltypegenetraitcombosMultiVariant = set()

# gwasphenos = [""AD_EADB_GCST90704646""]
# phenoObj = phenomap.get("AD_EADB_GCST90704646")
# print(f"{len(phenoObj.crediblesets)}")
# count = set()
# for cs in phenoObj.crediblesets:
#     if not cs.singlevariant:
#         count.add(cs)
#         print(cs.toStr())
# print(f"{len(count)}")
# exit()


for pheno in gwasphenos:
    
    phenoObj = phenomap.get(pheno)

    lociWithSignificantColocGene = set()
    lociWithSignificantColocGeneSingleVariant = set()
    lociWithSignificantColocGeneMultiVariant = set()
    for coloc in phenoObj.colocresults:
        if coloc.significant:
            for locus in phenoObj.loci:
                #print(f"{locus.chr}:{locus.start} - {coloc.gwasvarchr}:{coloc.gwasvarpos}")
                if locus.overlaps(coloc.gwasvarchr,coloc.gwasvarpos):
                    lociWithSignificantColocGene.add(locus)
                    if coloc.singlevariant:
                        lociWithSignificantColocGeneSingleVariant.add(locus)
                    else:
                        lociWithSignificantColocGeneMultiVariant.add(locus)

    
    # nrgenes = phenoObj.countGenes(None, False)
    # nrgenesSingleVariant = phenoObj.countGenes(True, False)
    # nrgenesMultiVariant = phenoObj.countGenes(False, False)
    # nrcrediblesets = phenoObj.countCredibleSets()   

    propA = round((len(lociWithSignificantColocGene)/len(phenoObj.loci))*100,2) if len(phenoObj.loci) > 0 else 0
    propB = round((len(lociWithSignificantColocGeneSingleVariant)/len(phenoObj.loci))*100,2) if len(phenoObj.loci) > 0 else 0
    propC = round((len(lociWithSignificantColocGeneMultiVariant)/len(phenoObj.loci))*100,2) if len(phenoObj.loci) > 0 else 0

    lnout = (f"{pheno}"
        f"\t{len(phenoObj.loci)}\t{len(lociWithSignificantColocGene)}\t{len(lociWithSignificantColocGeneSingleVariant)}\t{len(lociWithSignificantColocGeneMultiVariant)}"
        f"\t{propA}\t{propB}\t{propC}" )

    for coloc in phenoObj.colocresults:
        testedGenesOverall.add(coloc.gene)
        if coloc.singlevariant:
            testedGenesSingleVariant.add(coloc.gene)
            if coloc.significant:
                colocGenesSingleVariant.add(coloc.gene)
                colocGenesOverall.add(coloc.gene)
        else:
            testedGenesMultiVariant.add(coloc.gene)
            if coloc.significant:
                colocGenesMultiVariant.add(coloc.gene)
                colocGenesOverall.add(coloc.gene)
    
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

    lnout += (f"\t{testsoverall}\t{testssinglevariant}\t{testsmultivariant}"
                f"\t{testedgenesoverall}\t{testedgenessinglevariant}\t{testedgenesmultivariant}"
                f"\t{testedcrediblesets[1]}\t{testedcrediblesets[0]}"

                f"\t{signtestsoverall}\t{signtestssinglevariant}\t{signtestsmultivariant}"
                f"\t{signgenesoverall}\t{signgenessinglevariant}\t{signgenesmultivariant}"
                f"\t{signcrediblesets[1]}\t{signcrediblesets[0]}"
    )
    fho.write(lnout+"\n")
    print(lnout)
# exit()
fho.write("\n")
fho.write("\n")
print()
print()

propA = round((len(colocGenesOverall) / len(testedGenesOverall)) * 100,2)
propB = round((len(colocGenesSingleVariant) / len(testedGenesSingleVariant)) * 100,2)
propC = round((len(colocGenesMultiVariant) / len(testedGenesMultiVariant)) * 100,2)
outln = f"Nr of colocalizing / tested genes = (%) (overall): {len(colocGenesOverall)} / {len(testedGenesOverall)} = {propA}%"
print(outln)
fho.write(outln+"\n")
outln =f"Nr of colocalizing / tested genes = (%) (single variant coloc): {len(colocGenesSingleVariant)} / {len(testedGenesSingleVariant)} = {propB}%"
print(outln)
fho.write(outln+"\n")
outln =f"Nr of colocalizing / tested genes = (%) (multi variant coloc): {len(colocGenesMultiVariant)} / {len(testedGenesMultiVariant)} = {propC}%"
print(outln)
fho.write(outln+"\n")






##############################################################
# gwas: determining the single cell specific colocalizations
##############################################################
fho.write("##############################################################\n"
          "# gwas: determining the single cell specific colocalizations\n"
          "##############################################################\n\n")
fho.write("\n")
fho.write("\n")
print()
print()

fho.write("\n")
header=("GWAS"
        "\tNumber of significant colocalizing genes"
        "\tNumber of significant colocalizing genes (single variant coloc)"
        "\tNumber of significant colocalizing genes (multi variant coloc)"
        "\tNumber of significant colocalizing credible sets"

        "\tNumber of significant colocalizing genes in a single cell type"
        "\tNumber of significant colocalizing genes in a single cell type (single variant coloc)"
        "\tNumber of significant colocalizing genes in a single cell type (multi variant coloc)"
        "\tNumber of significant colocalizing credible sets in a single cell type (GWAS)"

        "\tPercentage of significant colocalizing genes in a single cell type "
        "\tPercentage of significant colocalizing genes in a single cell type (single variant coloc)"
        "\tPercentage of significant colocalizing genes in a single cell type (multi variant coloc)"
        "\tPercentage of significant colocalizing credible sets in a single cell type (GWAS)"
            )
print(header)
fho.write(header+"\n")

celltypespecgenesOverall = set()
celltypespecgenesSingleVariant = set()
celltypespecgenesMultiVariant = set()

for pheno in gwasphenos:
    phenoObj = phenomap.get(pheno)

    # credible sets are unique per gene; shared between cell types.
    nrSignGenesOverall = phenoObj.countColocTests(useSingleVariantFilter=None, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
    nrSignGenesSingleVariant = phenoObj.countColocTests(useSingleVariantFilter=True, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
    nrSignGenesMultiVariant = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
    nrSignCredibleSets = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=True, countGenes=False, countCredibleSets = True, filterPhenoList=gwasphenoset)
    nrSignCredibleSets = nrSignCredibleSets[0]

    cellTypeSpecificGenesOverall = phenoObj.determineCellTypeSpecificColocalizations(useSingleVariantFilter=None, countGenes=True, countCredibleSets=False)
    cellTypeSpecificGenesSingleVariant = phenoObj.determineCellTypeSpecificColocalizations(useSingleVariantFilter=True, countGenes=True, countCredibleSets=False)
    cellTypeSpecificGenesMultiVariant  = phenoObj.determineCellTypeSpecificColocalizations(useSingleVariantFilter=False, countGenes=True, countCredibleSets=False)
    cellTypeSpecificCredibleSetsMultiVariant  = phenoObj.determineCellTypeSpecificColocalizations(useSingleVariantFilter=False, countGenes=False, countCredibleSets=True)

    celltypespecgenesOverall.update(cellTypeSpecificGenesOverall)
    celltypespecgenesSingleVariant.update(cellTypeSpecificGenesSingleVariant)
    celltypespecgenesMultiVariant.update(cellTypeSpecificGenesMultiVariant)

    nrCellTypeSpecificGenesOverall = len(cellTypeSpecificGenesOverall)
    nrCellTypeSpecificGenesSingleVariant = len(cellTypeSpecificGenesSingleVariant) 
    nrCellTypeSpecificGenesMultiVariant = len(cellTypeSpecificGenesMultiVariant)
    nrCellTypeSpecificCredibleSetsMultiVariant = len(cellTypeSpecificCredibleSetsMultiVariant)

    propA = round((nrCellTypeSpecificGenesOverall/nrSignGenesOverall) *100 ,2) if nrSignGenesOverall > 0 else 0
    propB = round((nrCellTypeSpecificGenesSingleVariant/nrSignGenesSingleVariant) *100 ,2) if nrSignGenesSingleVariant > 0 else 0
    propC = round((nrCellTypeSpecificGenesMultiVariant/nrSignGenesMultiVariant) *100 ,2) if nrSignGenesMultiVariant > 0 else 0
    propD = round((nrCellTypeSpecificCredibleSetsMultiVariant/nrSignCredibleSets) *100 ,2) if nrSignCredibleSets > 0 else 0  

    outln =     (f"{pheno}\t{nrSignGenesOverall}\t{nrSignGenesSingleVariant}\t{nrSignGenesMultiVariant}\t{nrSignCredibleSets}"
                f"\t{nrCellTypeSpecificGenesOverall}\t{nrCellTypeSpecificGenesSingleVariant}\t{nrCellTypeSpecificGenesMultiVariant}\t{nrCellTypeSpecificCredibleSetsMultiVariant}"
                f"\t{propA}\t{propB}\t{propC}\t{propD}"
             )
    print(outln)
    fho.write(outln+"\n")
fho.write("\n")
fho.write("\n")
print()
print()

##############################################################
# gwas: determining the single cell specific colocalizations
##############################################################
fho.write("##############################################################\n"
          "# gwas: determining the single cell specific colocalizations\n"
          "##############################################################\n\n")
# check whether the uniquely colocalizing genes 
# are not a gwas coloc result gene for another GWAS
celltypesPerGeneOverall = {}
celltypesPerGeneSingleVariant = {}
celltypesPerGeneMultiVariant = {}
for celltype in celltypes:
    phenoObj = phenomap.get(celltype)
    for coloc in phenoObj.colocresults:
        gene = coloc.gene
        # if coloc.significant:
        #     print(f"{celltype} - {coloc.gene} - {coloc.gene in celltypespecgenesOverall} - {coloc.hasPhenotypeInSet(gwasphenos)}")
        if coloc.significant and coloc.gene in celltypespecgenesOverall and coloc.hasPhenotypeInSet(gwasphenos):
            
            other = celltypesPerGeneOverall.get(gene)
            if other is None:
                other = set()
                celltypesPerGeneOverall[gene] = other
            other.add(celltype)
        if coloc.singlevariant and coloc.significant and coloc.gene in celltypespecgenesSingleVariant and coloc.hasPhenotypeInSet(gwasphenos):            
            other = celltypesPerGeneSingleVariant.get(gene)
            if other is None:
                other = set()
                celltypesPerGeneSingleVariant[gene] = other
            other.add(celltype)
        if not coloc.singlevariant and coloc.significant and coloc.gene in celltypespecgenesMultiVariant and coloc.hasPhenotypeInSet(gwasphenos):            
            other = celltypesPerGeneMultiVariant.get(gene)
            if other is None:
                other = set()
                celltypesPerGeneMultiVariant[gene] = other
            other.add(celltype)
    # exit()

# double check whether the credible sets do not colocalize in another cell type for another GWAS
ctrOverall = 0
for gene, coloccelltypes in celltypesPerGeneOverall.items():
    # print(f"{gene} - {coloccelltypes}")
    if len(coloccelltypes) == 1:
        ctrOverall += 1
ctrSingleVariant = 0
for gene, coloccelltypes in celltypesPerGeneSingleVariant.items():
    if len(coloccelltypes) == 1:
        ctrSingleVariant += 1
ctrMultiVariant = 0
for gene, coloccelltypes in celltypesPerGeneMultiVariant.items():
    if len(coloccelltypes) == 1:
        ctrMultiVariant += 1

propA = round( (ctrOverall / len(colocGenesOverall)) * 100 ,2 )
propB = round( (ctrSingleVariant / len(colocGenesSingleVariant)) * 100, 2 )
propC = round( (ctrMultiVariant / len(colocGenesMultiVariant)) * 100, 2 )

outln = f"Nr of (overall) genes colocalizing in single cell type / colocalizing = (%): {ctrOverall} / {len(colocGenesOverall)} = {propA}%"
print(outln)
fho.write(outln+"\n")
outln =f"Nr of (single variant coloc) genes colocalizing in single cell type / colocalizing = (%): {ctrSingleVariant} / {len(colocGenesSingleVariant)} = {propB}%"
print(outln)
fho.write(outln+"\n")
outln =f"Nr of (multi variant coloc) genes colocalizing in single cell type / colocalizing = (%) : {ctrMultiVariant} / {len(colocGenesMultiVariant)} = {propC}%"
print(outln)
fho.write(outln+"\n")


fho.write("##############################################################\n"
          "# pairwise comparisons; gwas to cell types and vice versa\n"
          "##############################################################\n\n")

# pairwise comparisons
for pair in [ [gwasphenos, celltypephenos], [celltypephenos, gwasphenos] ]:
    ref = pair[0]
    alt = pair[1]
    for singlevariantFilter in [False, True]:
        significantFilter = True
        
        print()
        print()

        header = "Pheno\tTotal number of colocalization tests\tTotal number of genes"
        if not singlevariantFilter:
            header+="\tTotal number of credible sets"
        header += "\tTotal number of significant colocalization tests\tTotal number of significant genes"
        if not singlevariantFilter:
            header+="\tTotal number of significant credible sets"
        header += "\t"

        for pheno in alt:
            if singlevariantFilter:
                header+="\t"+pheno+" colocalizing genes"
            else:
                header+="\t"+pheno+" colocalizing genes\t"+pheno+" colocalizing tests\t"+pheno+" colocalizing credible sets"
        
        print(header)
        fho.write("\n")
        fho.write(f"\nSingleVariant: {singlevariantFilter}, significant: {significantFilter}\n")
        fho.write(header+"\n")
        significant = True
        for pheno1 in ref:
            phenoObj1 = phenomap.get(pheno1)
            nrColocTests = phenoObj1.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=False, countGenes=False, filterPhenoList=gwasphenoset)
            nrGenes = phenoObj1.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=False, countGenes=True, filterPhenoList=gwasphenoset)
            nrSignColocTests = phenoObj1.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes=False, filterPhenoList=gwasphenoset)
            nrSignGenes = phenoObj1.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes=True, filterPhenoList=gwasphenoset)

            outln = f"{pheno1}\t{nrColocTests}\t{nrGenes}"

            if not singlevariantFilter:
                outln+= f"\t{phenoObj1.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=False, countGenes=False, countCredibleSets=True, filterPhenoList=gwasphenoset)[0]}"

            outln += f"\t{nrSignColocTests}\t{nrSignGenes}"
            if not singlevariantFilter:
                outln+= f"\t{phenoObj1.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes=False, countCredibleSets=True, filterPhenoList=gwasphenoset)[0]}"

            outln += "\t"
            for pheno2 in alt:
                phenoObj2 = phenomap.get(pheno2)
                overlapGenes = phenoObj1.determineOverlap(phenoObj2, useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes = True)
            
                if singlevariantFilter:
                    outln += f"\t{len(overlapGenes)}"
                else:
                    overlapTests  = phenoObj1.determineOverlap(phenoObj2, useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes = False)
                    overlapCredibleSets  = phenoObj1.determineOverlap(phenoObj2, useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes = False, countCredibleSets = True)
                    outln += f"\t{len(overlapGenes)}\t{len(overlapTests)}\t{len(overlapCredibleSets)}"
            print(outln)
            fho.write(outln+"\n")

fho.write("\n")
fho.write("\n")
print()
print()


fho.write("##############################################################\n"
          "# overlaps within each cell type / gwas\n"
          "##############################################################\n\n")

for pair in [ [gwasphenos, celltypephenos], [celltypephenos, gwasphenos] ]:
    ref = pair[0]
    alt = pair[1]

    for singlevariantFilter in [False, True]:
        if singlevariantFilter:
            header = "Cell type\tNumber of colocalizing genes"
            for pheno in ref:
                header += "\t"+pheno
            fho.write(header+"\n")
            print(header)
            for phenoA in ref:
                phenoObj1 = phenomap.get(phenoA)       
                allSignificantGenes = set()
                outln = f"{phenoA}"
                for coloc in phenoObj1.colocresults:
                    if coloc.singlevariant and coloc.significant and coloc.hasPhenotypeInSet(alt):
                        allSignificantGenes.add(coloc.gene)
                outln += f"\t{len(allSignificantGenes)}"
                for phenoB in ref:
                    phenoObj2 = phenomap.get(phenoB)
                    overlap = set()
                    for coloc in phenoObj2.colocresults:
                        if coloc.singlevariant and coloc.significant and coloc.hasPhenotypeInSet(alt) and coloc.gene in allSignificantGenes:
                            overlap.add(coloc.gene)
                    outln += f"\t{len(overlap)}"
                print(outln)
                
                fho.write(outln+"\n")
            fho.write("\n")
            fho.write("\n")
            print()
            print()
        else:
            header = "Cell type\tNumber of colocalizing genes\tNumber of colocalizing credible sets"
            for phenoA in ref:
                header += "\tOverlapping genes "+phenoA+"\tOverlapping credible sets "+phenoA
            fho.write(header+"\n")
            print(header)

            for phenoA in ref:
                phenoObj1 = phenomap.get(phenoA)
                allSignificantGenes = set()
                allSignificantCredibleSets = set()

                outln = f"{phenoA}"
                for coloc in phenoObj1.colocresults:
                    if not coloc.singlevariant and coloc.significant and coloc.hasPhenotypeInSet(alt):
                        allSignificantGenes.add(coloc.gene)
                        allSignificantCredibleSets.add(coloc.getOtherCredibleSetByStr(phenoA))

                outln += f"\t{len(allSignificantGenes)}\t{len(allSignificantCredibleSets)}"
                for phenoB in ref:
                    phenoObj2 = phenomap.get(phenoB)
                    overlapGWASCredibleSets = set()
                    overlapGenes = set()
                    for coloc in phenoObj2.colocresults:
                        if not coloc.singlevariant and coloc.significant and coloc.hasPhenotypeInSet(alt) and coloc.gene in allSignificantGenes:
                            gwasCs = coloc.getOtherCredibleSetByStr(phenoB)
                            if gwasCs in allSignificantCredibleSets:
                                overlapGWASCredibleSets.add(gwasCs)
                            overlapGenes.add(coloc.gene)
                    outln += f"\t{len(overlapGenes)}\t{len(overlapGWASCredibleSets)}"
                print(outln)
                fho.write(outln+"\n")
            fho.write("\n")
            fho.write("\n")
            print()
            print()

fho.write("\n")
fho.write("\n")
print()
print()



fho.write("##############################################################\n"
          "# comparisons with bulk MetaBrain, per GWAS\n"
          "##############################################################\n\n")

# gwasphenos = ["AD_EADB_GCST90704646"]
header = (
    "GWAS"
    "\tNumber of significantly colocaling genes"
    "\tNumber of significantly colocaling genes (single variant coloc)"
    "\tNumber of significantly colocaling genes (multi variant coloc)"
    "\tNumber of significantly colocaling credible sets (multi variant coloc)"

    "\tNumber of tested genes in bulk"
    "\tNumber of tested genes in bulk (single variant coloc)"
    "\tNumber of tested genes in bulk (multi variant coloc)"
    "\tNumber of tested credible sets in bulk (multi variant coloc)"

    "\tNumber of significantly colocalizing genes in bulk"
    "\tNumber of significantly colocalizing genes in bulk (single variant coloc)"
    "\tNumber of significantly colocalizing genes in bulk (multi variant coloc)"
    "\tNumber of significantly colocalizing credible sets in bulk (multi variant coloc)"
    "\t"
)
for bulkpheno in bulkphenos:
    header += ( f"\t{bulkpheno} - Significant overlapping genes"
                f"\t{bulkpheno} - Significant overlapping genes (single variant coloc)"
                f"\t{bulkpheno} - Significant overlapping genes (multi variant coloc)"
                f"\t{bulkpheno} - Significant overlapping credible sets (multi variant coloc)"
                )

print(header)
fho.write(header+"\n")
for pheno in gwasphenos:
    phenoObj = phenomap.get(pheno)

    # credible sets are unique per gene; shared between cell types.
    nrSignGenesOverall = phenoObj.countColocTests(useSingleVariantFilter=None, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
    nrSignGenesSingleVariant = phenoObj.countColocTests(useSingleVariantFilter=True, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
    nrSignGenesMultiVariant = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=True, countGenes=True, filterPhenoList=gwasphenoset)
    nrSignCredibleSets = phenoObj.countColocTests(useSingleVariantFilter=False, useSignificanceFilter=True, countGenes=False, countCredibleSets = True, filterPhenoList=gwasphenoset)
    nrSignCredibleSets = nrSignCredibleSets[0]
    
    testedGenes = set()
    testedGenesSingleVariant = set()
    testedGenesMultiVariant = set()
    testedCredibleSets = set()

    significantGenes = set()
    significantGenesSingleVariant = set()
    significantGenesMultiVariant = set()
    significantCredibleSets = set()

    overlapPerBulkPheno = {}

    for colocresult in phenoObj.colocresults:
        if colocresult.significant:
            # get cell type eQTL credible set
            phenoCs = colocresult.getCredibleSetByStr(pheno)
            celltypeCs = colocresult.getOtherCredibleSetByStr(pheno)
            
            
            # determine if celltype has matching coloc result for this gene in bulk
            celltypeColocResults = celltypeCs.phenotype.colocresults
            for colocresultB in celltypeColocResults:
                # if colocresultB.hasCredibleSet(celltypeCs):
                #     print(f"{celltypeCs.toStr()} - {colocresultB.toStr()} - celltypecs: {colocresultB.hasCredibleSet(celltypeCs)} - bulk: {colocresultB.hasPhenotypeInSet(bulkphenos)} ")

                # determine if it is for the same cell type credible set   
                if colocresultB.hasCredibleSet(celltypeCs) and colocresultB.hasPhenotypeInSet(bulkphenos):
                    
                    gene = colocresultB.gene
                    testedGenes.add(gene)

                    if colocresultB.significant:
                        colocresultB.gwascredibleset = phenoCs
                        significantGenes.add(gene)
                        bulkCs = colocresultB.getOtherCredibleSetByCsObj(celltypeCs)
                        other = overlapPerBulkPheno.get(bulkCs.phenotype.name)
                        if other is None:
                            other = set()
                            overlapPerBulkPheno[bulkCs.phenotype.name] = other
                        other.add(colocresultB)

                    if colocresultB.singlevariant:
                        testedGenesSingleVariant.add(gene)
                        if colocresultB.significant:
                            significantGenesSingleVariant.add(gene)
                    else:
                        testedGenesMultiVariant.add(gene)
                        testedCredibleSets.add(phenoCs)
                        if colocresultB.significant:
                            significantGenesMultiVariant.add(gene)
                            significantCredibleSets.add(phenoCs)
    
    bulkstr = ""
    for bulkpheno in bulkphenos:
        overlap = overlapPerBulkPheno.get(bulkpheno)
        
        genes = set()
        genessinglevariant = set()
        genesmultivariant = set()
        crediblesets = set()
        if overlap is not None:
            for coloc in overlap:
                genes.add(coloc.gene)
                if coloc.singlevariant:
                    genessinglevariant.add(coloc.gene)
                else:
                    genesmultivariant.add(coloc.gene)
                    crediblesets.add(coloc.gwascredibleset)
        bulkstr += f"\t{len(genes)}\t{len(genessinglevariant)}\t{len(genesmultivariant)}\t{len(crediblesets)}"

    outln = (f"{pheno}"
            f"\t{nrSignGenesOverall}\t{nrSignGenesSingleVariant}\t{nrSignGenesMultiVariant}\t{nrSignCredibleSets}"
            f"\t{len(testedGenes)}\t{len(testedGenesSingleVariant)}\t{len(testedGenesMultiVariant)}\t{len(testedCredibleSets)}"
            f"\t{len(significantGenes)}\t{len(significantGenesSingleVariant)}\t{len(significantGenesMultiVariant)}\t{len(significantCredibleSets)}"
            f"\t{bulkstr}"
            )
    print(outln)
    fho.write(outln+"\n")
    

fho.write("\n")
fho.write("\n")
print()
print()


fho.write("##############################################################\n"
          "# comparisons with other, publishes datasets\n"
          "##############################################################\n\n")

signGenesBryois, eGenesBryois, signGenesFujita, eGenesFujita, signGenesSingleBrain, eGenesSingleBrain = getEGeneListsFromOtherPubs(celltypes)

header = ("Pheno"
        "\tSignificantly colocalizing genes"
        "\tBryois overlapping significant eGenes"
        "\tFujita overlapping significant eGenes"
        "\tSingleBrain overlapping significant eGenes"
        "\tSignificantly colocalizing genes not overlapping previously published eGenes"
)

fho.write(header+"\n")
print(header)
for pheno in gwasphenos:
    phenoObj = phenomap.get(pheno)
    allgenes = set()
    for coloc in phenoObj.colocresults:
        if coloc.significant:
            scCs = coloc.getOtherCredibleSetByStr(pheno)
            allgenes.add(coloc.gene)

    unique = 0
    overlapBryois = 0
    overlapFujita = 0
    overlapSingleBrain = 0
    for gene in allgenes:
        if gene not in eGenesFujita and gene not in eGenesBryois and gene not in eGenesSingleBrain:
            unique += 1
        if gene in eGenesFujita:
            overlapFujita += 1
        if gene in eGenesSingleBrain:
            overlapSingleBrain+= 1
        if gene in eGenesBryois:
            overlapBryois+= 1
    outln = f"{pheno}\t{len(allgenes)}\t{overlapBryois}\t{overlapFujita}\t{overlapSingleBrain}\t{unique}"
    fho.write(outln+"\n")
    print(outln)

fho.write("\n")
fho.write("\n")
print()
print()

header = ("Cell type"
        "\tSignificantly colocalizing genes"
        "\tBryois overlapping significant eGenes"
        "\tFujita overlapping significant eGenes"
        "\tSingleBrain overlapping significant eGenes"
        "\tSignificantly colocalizing genes not overlapping previously published eGenes"
)
fho.write(header+"\n")
print(header)
for celltype in celltypes:
    phenoObj = phenomap.get(celltype)

    allgenes = set()
    for coloc in phenoObj.colocresults:
        if coloc.significant and coloc.hasPhenotypeInSet(gwasphenos):
            allgenes.add(coloc.gene)
    unique = 0
    overlapBryois = 0
    overlapFujita = 0
    overlapSingleBrain = 0
    othereGenesBryois = signGenesBryois.get(celltype)
    othereGenesFujita = signGenesFujita.get(celltype)
    othereGenesSingleBrain = signGenesSingleBrain.get(celltype)
    totaloverlapping = 0
    for gene in allgenes:
        isOverlapping = False
        if othereGenesBryois is not None and gene in othereGenesBryois:
            isOverlapping = True
            overlapBryois += 1
        if othereGenesFujita is not None and gene in othereGenesFujita:
            isOverlapping = True
            overlapFujita += 1
        if othereGenesSingleBrain is not None and gene in othereGenesSingleBrain:
            isOverlapping = True
            overlapSingleBrain += 1
        if isOverlapping:
            totaloverlapping += 1
    unique = len(allgenes) - totaloverlapping
    outln = f"{celltype}\t{len(allgenes)}\t{overlapBryois}\t{overlapFujita}\t{overlapSingleBrain}\t{unique}"
    fho.write(outln+"\n")    
    print(outln)
    

fho.write("\n")
fho.write("\n")
print()
print()



fho.write("\n")
fho.write("\n")
print()
print()
        

fho.write("##############################################################\n"
          "# list overlapping genes per gwas pair\n"
          "##############################################################\n\n")
# Coloc overlap between GWASes
# Disease1	Disease2	Nr Shared	Overlapping Coloc Genes
header = "GWAS1\tGWAS2\tNr shared single variant\tNr shared multi variant\tShared genes single variant\tShared genes multi variant"
fho.write(header+"\n")
print(header)

multiGWASgenes = set()
for i in range(0,len(gwasphenos)):
    phenoAname = gwasphenos[i]
    phenoA = phenomap.get(phenoAname)

    genesSingleVariant = set()
    genesMultiVariant = set()
    for coloc in phenoA.colocresults:
        if coloc.significant:
            if coloc.singlevariant:
                genesSingleVariant.add(coloc.gene)
            else:
                genesMultiVariant.add(coloc.gene)

    for j in range(i+1,len(gwasphenos)): 
        phenoBname = gwasphenos[j] 
        phenoB = phenomap.get(phenoBname)
        overlapSingleVariant = set()
        overlapMultiVariant = set()
        for coloc in phenoB.colocresults:
            if coloc.significant:
                if coloc.singlevariant and coloc.gene in genesSingleVariant:
                    overlapSingleVariant.add(coloc.gene)
                    multiGWASgenes.add(coloc.gene)
                elif not coloc.singlevariant and coloc.gene in genesMultiVariant:
                    overlapMultiVariant.add(coloc.gene)
                    multiGWASgenes.add(coloc.gene)
        listA = list(overlapSingleVariant)
        listA.sort()
        listAStr = ";".join(listA)
        listB = list(overlapMultiVariant)
        listB.sort()
        listBStr = ";".join(listB)
        if len(overlapSingleVariant) > 0 or len(overlapMultiVariant) > 0:
            
            outln = f"{phenoAname}\t{phenoBname}\t{len(overlapSingleVariant)}\t{len(overlapMultiVariant)}\t{listAStr}\t{listBStr}"
            print(outln)
            fho.write(outln+"\n")

fho.close()
print(f"{len(multiGWASgenes)} genes colocalize in > 1 GWAS")

# todo:
# - intersect with bulk
# - disease specific summaries (how?)