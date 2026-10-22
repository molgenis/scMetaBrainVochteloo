

class ColocResult:
    def __init__(self, credibleSetA, credibleSetB, method, fdr, significant):
        self.credibleSetA = credibleSetA
        self.credibleSetB = credibleSetB
        self.singlevariant = (method == "coloc.abf()")
        self.fdr = fdr
        self.significant = (significant == "True")
    
    def toStr(self):
        print(f"{self.credibleSetA.gene} - {self.credibleSetA.number} - {self.singlevariant}")

class CredibleSet:
    def __init__(self, gene, phenotype, number):
        self.phenotype = phenotype
        self.gene = gene
        self.number = number
        self.singlevariant = (number == "-")
        self.coloc = set()

    def addColoc(self, colocresult):
        self.coloc.add(colocresult)

    def toStr(self):
        return self.phenotype.name+":"+self.gene+":"+self.number

class Phenotype:
    def __init__(self, name):
        self.name = name
        self.colocresults = []
        self.crediblesets = set()

    def countGenes(self, useSingleVariantFilter, useSignificanceFilter):
        genes = set()
        for coloc in self.colocresults:
            if not useSignificanceFilter or (useSignificanceFilter and coloc.significant):
                if useSingleVariantFilter is None:
                    genes.add(coloc.credibleSetA.gene)
                elif useSingleVariantFilter and coloc.singlevariant:
                    genes.add(coloc.credibleSetA.gene)
                elif not useSingleVariantFilter and not coloc.singlevariant:
                    genes.add(coloc.credibleSetA.gene)
        return len(genes)
    
    def countCredibleSets(self):
        csset = set()
        for cs in self.crediblesets:
            if not cs.singlevariant:
                csset.add(cs)
        return len(csset)

    def countColocTests(self, useSingleVariantFilter, useSignificanceFilter, countGenes, countCredibleSets=False):
        coloctests = set()
        for coloc in self.colocresults:
            if ((not useSignificanceFilter or (useSignificanceFilter and coloc.significant)) and 
                    ( ( useSingleVariantFilter is None)
                        or (useSingleVariantFilter and coloc.singlevariant)
                        or (not useSingleVariantFilter and not coloc.singlevariant)
                    )):
                if countGenes:
                    coloctests.add(coloc.credibleSetA.gene)
                elif countCredibleSets:
                    if coloc.credibleSetA.phenotype == self:
                        coloctests.add(coloc.credibleSetA)
                    else:
                        coloctests.add(coloc.credibleSetB)
                else:
                    coloctests.add(coloc)
        return len(coloctests)

    def countNonColoc(self, countGenes=False, useSingleVariantFilter=False):
        sign = set()
        other = set()
        

        
        for coloc in self.colocresults:  
            obj = coloc.credibleSetA
            if (coloc.credibleSetB.phenotype == self):
                obj = coloc.credibleSetB

            if ( ( useSingleVariantFilter is None)
                    or (useSingleVariantFilter and coloc.singlevariant)
                    or (not useSingleVariantFilter and not coloc.singlevariant)
                ):
                if countGenes:
                    if coloc.significant:
                        sign.add(obj.gene)
                    other.add(obj.gene)
                else:
                    
                    if coloc.significant:
                        sign.add(obj)

        unique = set()
        if countGenes:
            for gene in other:
                if gene not in sign:
                    unique.add(gene)
        else:
            for cs in self.crediblesets:
                if cs not in sign:
                    unique.add(cs)

        return len(unique)        

    def determineOverlap(self, phenoObj2, useSingleVariantFilter, useSignificanceFilter, countGenes, countCredibleSets=False):
        overlap = set()
        for coloc in self.colocresults:  
            if (coloc.credibleSetA.phenotype == self and coloc.credibleSetB.phenotype == phenoObj2) or (coloc.credibleSetB.phenotype == self and coloc.credibleSetA.phenotype == phenoObj2):
                if ((not useSignificanceFilter or (useSignificanceFilter and coloc.significant)) and 
                    ( ( useSingleVariantFilter is None)
                        or (useSingleVariantFilter and coloc.singlevariant)
                        or (not useSingleVariantFilter and not coloc.singlevariant)
                    )):
                    if countGenes:
                        overlap.add(coloc.credibleSetA.gene)
                    elif countCredibleSets:
                        if coloc.credibleSetA.phenotype == self:
                            overlap.add(coloc.credibleSetA)
                        else:
                            overlap.add(coloc.credibleSetB)
                    else:
                        overlap.add(coloc)
        return overlap

def index(elems):
    out = {}
    for i in range(len(elems)):
        out[elems[i]] = i
    return out

file = "coloc_results_all_datasets.zero_cs_removed.wfdr.txt"
file = "/path/to/coloc_results_all_bulk_tissues_all_cell_types.wfdr.txt"

# Gene
# Dataset 1
# Dataset 2
# Number of dataset 1 variants
# Number of dataset 2 variants
# Number of overlapping variants
# Dataset 1 variant
# Dataset 1 variant chromosome
# Dataset 1 variant position
# Dataset 1 variant p-value
# Dataset 1 variant beta
# Dataset 1 variant standard error
# Dataset 2 variant
# Dataset 2 variant chromosome
# Dataset 2 variant position
# Dataset 2 variant p-value
# Dataset 2 variant beta
# Dataset 2 variant standard error
# Distance between variants
# PPH0
# PPH1
# PPH2
# PPH3
# PPH4
# Dataset 1 credible set
# Dataset 2 credible set
# Dataset 1 lambda
# Dataset 2 lambda
# Dataset 1 variant is significant
# Dataset 2 variant is significant
# Analysis type
# FDR
# Significant

phenomap = {}
csmap = {}
colocresults = []

fh = open(file,'r')
header = fh.readline().rstrip().split("\t")
idx = index(header)
for k,v in idx.items():
    print(f"{k}\t{v}")
print()

# parse file
for line in fh:
    elems = line.strip().split("\t")
    gene = elems[idx['Gene']]
    pheno1 = elems[idx['Cell type']]
    pheno2 = elems[idx['Bulk tissue']]
    cs1nr = elems[idx['SC credible set']]
    cs2nr = elems[idx['Bulk credible set']]
    method = elems[idx['Analysis type']]
    
    fdr = float(elems[idx['FDR']])
    sign = elems[idx['Significant']]

    phenObj1 = phenomap.get(pheno1)
    phenObj2 = phenomap.get(pheno2)
    if phenObj1 is None:
        phenObj1 = Phenotype(pheno1)
        phenomap[pheno1] = phenObj1
        
    if phenObj2 is None:
        phenObj2 = Phenotype(pheno2)
        phenomap[pheno2] = phenObj2
    
    csObj1 = csmap.get(pheno1+":"+gene+":"+cs1nr)
    if csObj1 is None:
        csObj1 = CredibleSet(gene, phenObj1, cs1nr)
        csmap[pheno1+":"+gene+":"+cs1nr] = csObj1
        phenObj1.crediblesets.add(csObj1)

    csObj2 = csmap.get(pheno2+":"+gene+":"+cs2nr)
    if csObj2 is None:
        csObj2 = CredibleSet(gene, phenObj2, cs2nr)
        csmap[pheno2+":"+gene+":"+cs2nr] = csObj2
        phenObj2.crediblesets.add(csObj2)
    
    

    colocObj = ColocResult(csObj1, csObj2, method, fdr, sign)
    # coloc is bidirectional
    phenObj1.colocresults.append(colocObj)
    phenObj2.colocresults.append(colocObj)
    csObj1.addColoc(colocObj)
    csObj2.addColoc(colocObj)
    colocresults.append(colocObj)

print(f"{len(colocresults)} for {len(phenomap)} phenotypes.")

# count several things
phenos = list(phenomap.keys())
phenos.sort()

fho = open('summary.txt','w')
header = ("Phenotype"
        "\tNumber of genes"
        "\tNumber of genes (single variant coloc)"
        "\tNumber of genes (multi variant coloc)"
        "\tNumber of credible sets"

        "\tNumber of overall coloc tests"
        "\tNumber of overall coloc genes"
        "\tNumber of overall significant coloc tests"
        "\tNumber of overall significant coloc genes"
        "\tNumber of multivariant coloc tests"
        "\tNumber of multivariant coloc genes"
        "\tNumber of multivariant coloc credible sets"
        "\tNumber of significant multivariant coloc tests"
        "\tNumber of significant multivariant coloc genes"
        "\tNumber of significant multivariant coloc credible sets"
        "\tNumber of singlevariant coloc tests"
        "\tNumber of singlevariant coloc genes"
        "\tNumber of significant singlevariant coloc tests"
        "\tNumber of significant singlevariant coloc genes"
        )

        # "\t"
        # "\tUnique genes"
        # "\tUnique genes (single variant)"
        # "\tUnique genes (multi variant)"
        # "\tUnique credible sets"
fho.write(header+"\n")
print(header)

celltypes = ["AST","END","EX","IN","MIC","OLI","OPC","PER"]
celltypeset = set()
celltypeset.update(celltypes)
celltypes = celltypeset

# per phenotype overview
for pheno in phenos:
        phenoObj = phenomap.get(pheno)

        nrgenes = phenoObj.countGenes(None, False)
        nrgenesSingleVariant = phenoObj.countGenes(True, False)
        nrgenesMultiVariant = phenoObj.countGenes(False, False)
        nrcrediblesets = phenoObj.countCredibleSets()   

        lnout = (f"{pheno}"
                f"\t{nrgenes}\t{nrgenesSingleVariant}\t{nrgenesMultiVariant}\t{nrcrediblesets}")

        for singlevariantFilter in [None, False, True]:
            for significantFilter in [False, True]:
                countTests = phenoObj.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes=False)
                countGenes = phenoObj.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes=True)
                lnout += f"\t{countTests}\t{countGenes}"
                if singlevariantFilter is not None and not singlevariantFilter:
                    countCredSets = phenoObj.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes=False, countCredibleSets = True)
                    lnout += f"\t{countCredSets}"
        
        # lnout += "\t"
        # # determine number of non-colocalizing genes
        # for singlevariantFilter in [None, False, True]:
        #     uniquegenes = phenoObj.countNonColoc(countGenes=False, useSingleVariantFilter=singlevariantFilter)
        #     lnout += f"\t{uniquegenes}"
        # uniquecrediblesets = phenoObj.countNonColoc(countGenes=True, useSingleVariantFilter=False)
        # lnout += f"\t{uniquecrediblesets}"
        # # determine number of non-colocalizing credible sets


        fho.write(lnout+"\n")
        print(lnout)


# pairwise comparisons  
for singlevariantFilter in [False, True]:
    for significantFilter in [False, True]:
    
        print()
        print()

        header = "Pheno\tTotal number of colocalization tests\tTotal number of genes\t"
        if not singlevariantFilter:
            header+="Total number of credible sets\t"
        for pheno in phenos:
            if pheno not in celltypes:
                if singlevariantFilter:
                    header+="\t"+pheno+" colocalizing genes"
                else:
                    header+="\t"+pheno+" colocalizing genes\t"+pheno+" colocalizing tests\t"+pheno+" colocalizing credible sets"
        print(header)
        fho.write("\n")
        fho.write(f"\nSingleVariant: {singlevariantFilter}, significant: {significantFilter}\n")
        fho.write(header+"\n")
        significant = True
        for pheno1 in phenos:
            if pheno1 in celltypes:
                phenoObj1 = phenomap.get(pheno1)
                nrColocTests = phenoObj1.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes=False)
                nrGenes = phenoObj1.countColocTests(useSingleVariantFilter=singlevariantFilter, useSignificanceFilter=significantFilter, countGenes=True)
                outln = f"{pheno1}\t{nrColocTests}\t{nrGenes}\t"
                if not singlevariantFilter:
                    outln+= f"{phenoObj1.countCredibleSets()}\t"

                for pheno2 in phenos:
                    if pheno2 not in celltypes:
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

fho.close()

