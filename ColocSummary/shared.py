import gzip

def fopen(file, mode):
	if file.endswith(".gz"):
		if mode == "r":
			return gzip.open(file,'rt')
		else:
			return gzip.open(file,'wt')
	return open(file, mode)

class Locus:
    def __init__(self, chr, start, stop):
        self.chr = chr
        self.start = start
        self.stop = stop

    def overlaps(self, chr, pos):
        if self.chr != chr:
            return False
        if pos >= self.start and pos <= self.stop:
            return True
        return False

def index(elems):
    out = {}
    for i in range(len(elems)):
        out[elems[i]] = i
    return out

class Gene:
	symbol = None
	ensg = None
	chr = None
	start = None
	stop = None	

class ColocResult:
    def __init__(self, credibleSetA, credibleSetB, method, fdr, significant):
        self.credibleSetA = credibleSetA
        self.credibleSetB = credibleSetB
        self.gene = credibleSetA.gene
        self.singlevariant = (method == "coloc.abf()")
        self.fdr = fdr
        self.significant = (significant == "True")
    
    def toStr(self):
        print(f"{self.credibleSetA.gene} - {self.credibleSetA.number} - {self.singlevariant}")

    def hasCredibleSet(self, credibleset):
        return (self.credibleSetA == credibleset or self.credibleSetB == credibleset)

    def hasPhenotypeInSet(self, phenotypeSet):
        return (self.credibleSetA.phenotype.name in phenotypeSet or self.credibleSetB.phenotype.name in phenotypeSet)

    def getOtherCredibleSetByStr(self, phenoname):
        if self.credibleSetA.phenotype.name == phenoname:
            return self.credibleSetB
        elif self.credibleSetB.phenotype.name == phenoname:
            return self.credibleSetA
        else:
            return None

    def getOtherCredibleSetByCsObj(self, credibleset):
        if self.credibleSetA == credibleset:
            return self.credibleSetB
        elif self.credibleSetB == credibleset:
            return self.credibleSetA
        else:
            return None

    def getCredibleSetByPhenoObj(self, phenoObj):
        if self.credibleSetA.phenotype == phenoObj:
            return self.credibleSetA
        elif self.credibleSetB.phenotype == phenoObj:
            return self.credibleSetB
        else:
            return None

    def getCredibleSetByStr(self, phenoname):
        if self.credibleSetA.phenotype.name == phenoname:
            return self.credibleSetA
        elif self.credibleSetB.phenotype.name == phenoname:
            return self.credibleSetB
        else:
            return None
    
    def toStr(self):
        return f"{self.credibleSetA.toStr()} - {self.credibleSetB.toStr()} - sign: {self.significant}"

class CredibleSet:
    def __init__(self, gene, phenotype, number):
        self.phenotype = phenotype
        self.gene = gene
        self.number = number
        self.singlevariant = (number == "-")
        self.colocresults = set()

    def addColoc(self, colocresult):
        self.colocresults.add(colocresult)

    def hasSignificantColoc(self):
        for coloc in self.colocresults:
            if coloc.significant:
                return True
        return False

    def toStr(self):
        return self.phenotype.name+":"+self.gene+":"+self.number

    def __lt__(self, other):
        return self.toStr() < other.toStr()

class Phenotype:
    def __init__(self, name):
        self.name = name
        self.colocresults = []
        self.crediblesets = set()
        self.loci = []

    def addLocusBreaker(self, locusbreakerfile):
        fh = open(locusbreakerfile,'r')
        header = fh.readline().rstrip().split("\t")
        idx = index(header)
        loci = []
        for line in fh:
            elems = line.strip().split("\t")
            chr = int(elems[idx['chr']])
            start = int(elems[idx['start']])
            stop = int(elems[idx['end']])
            locus = Locus(chr,start,stop)
            self.loci.append(locus)

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
    
    def countCredibleSets(self, onlySignificant=False, filterPhenoList=None):
        csset = set()
        for cs in self.crediblesets:
            if not cs.singlevariant:
                if not onlySignificant or cs.hasSignificantColoc():
                    csset.add(cs)
        return len(csset)

    def countColocTests(self, useSingleVariantFilter, useSignificanceFilter, countGenes, countCredibleSets=False, filterPhenoList=None):
        selfCs = set()
        otherCs = set()
        coloctests = set()
        for coloc in self.colocresults:
            if ((not useSignificanceFilter or (useSignificanceFilter and coloc.significant)) and 
                    ( ( useSingleVariantFilter is None)
                        or (useSingleVariantFilter and coloc.singlevariant)
                        or (not useSingleVariantFilter and not coloc.singlevariant)
                    )):
                if filterPhenoList is None or coloc.hasPhenotypeInSet(filterPhenoList):
                    if countGenes:
                        coloctests.add(coloc.gene)
                    elif countCredibleSets:
                        if coloc.credibleSetA.phenotype == self: # in case of GWAS comparisons this is the number of GWAS credible sets
                            selfCs.add(coloc.credibleSetA)
                            otherCs.add(coloc.credibleSetB)
                        else:
                            selfCs.add(coloc.credibleSetB)
                            otherCs.add(coloc.credibleSetA)
                    else:
                        coloctests.add(coloc)
        if countCredibleSets:
            return [len(selfCs), len(otherCs)] 
        return len(coloctests)

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
                        overlap.add(coloc.gene)
                    elif countCredibleSets:
                        if coloc.credibleSetA.phenotype == self:
                            overlap.add(coloc.credibleSetA)
                        else:
                            overlap.add(coloc.credibleSetB)
                    else:
                        overlap.add(coloc)
        return overlap

    def determineCellTypeSpecificColocalizations(self, useSingleVariantFilter, countGenes, countCredibleSets=False):
        celltypespecific = set()

        # determine the coloc results per gene
        geneToColocResult = {}
        for coloc in self.colocresults:
            if coloc.significant:
                if (useSingleVariantFilter is None) or (useSingleVariantFilter and coloc.singlevariant) or (not useSingleVariantFilter and not coloc.singlevariant):
                    other = geneToColocResult.get(coloc.credibleSetA.gene)
                    if other is None:
                        other = list()
                        geneToColocResult[coloc.credibleSetA.gene] = other
                    other.append(coloc)

        # iterate the genes
        for gene, genecoloclist in geneToColocResult.items():  
            if useSingleVariantFilter or countGenes:
                # count the cell types in which this gene is colocalizing
                coloccelltypes = set()
                for coloc in genecoloclist:
                    ctCs = coloc.credibleSetA
                    if ctCs.phenotype == self: # self in this case is a GWAS pheno object; the other must be a cell type
                        ctCs = coloc.credibleSetB
                    coloccelltypes.add(ctCs.phenotype)
                if len(coloccelltypes) == 1:
                    celltypespecific.add(gene)
            else:
                # group by GWAS credible set
                phenoCsToColocResult = {}
                for coloc in genecoloclist:
                    ctCs = coloc.credibleSetA
                    phenoCs = coloc.credibleSetB
                    if ctCs.phenotype == self: # self in this case is a GWAS pheno object; the other must be a cell type
                        ctCs = coloc.credibleSetB
                        phenoCs = coloc.credibleSetA
                    other = phenoCsToColocResult.get(phenoCs)
                    if other is None:
                        other = []
                        phenoCsToColocResult[phenoCs] = other
                    other.append(coloc)

                for phenoCs, coloclist in phenoCsToColocResult.items():
                    if len(coloclist) == 1:
                        celltypespecific.add(phenoCs)
                
        return celltypespecific



def getPhenoMap(file,ds1col,ds2col,cs1col,cs2col,existingphenomap):
# Gene
# Number of GWAS variants
# Number of eQTL variants
# Number of overlapping variants
# GWAS variant
# GWAS variant chromosome
# GWAS variant position
# GWAS variant p-value
# GWAS variant beta
# GWAS variant standard error
# eQTL variant
# eQTL variant chromosome
# eQTL variant position
# eQTL variant p-value
# eQTL variant beta
# eQTL variant standard error
# Distance between variants
# PPH0
# PPH1
# PPH2
# PPH3
# PPH4
# GWAS credible set
# eQTL credible set
# GWAS lambda
# eQTL lambda
# Cell type
# Nearest gene to GWAS variant
# Distance between GWAS variant and nearest gene
# Nearest gene is coloc gene
# GWAS variant is significant
# eQTL variant is significant
# Analysis type
# GWAS
# FDR
# Significant
    phenomap = None
    csmap = None
    if existingphenomap is not None:
        phenomap = existingphenomap
        # create cs map; not so efficient but whatever
        csmap = {}
        for phenoname, phenoObj in phenomap.items():
            for cs in phenoObj.crediblesets:
                csmap[f"{phenoname}:{cs.gene}:{cs.number}"] = cs
    else:
        phenomap = {}
        csmap = {}
    colocresults = []

    fh = open(file,'r')
    header = fh.readline().rstrip().split("\t")
    idx = index(header)
    # for k,v in idx.items():
    #     print(f"{k}\t{v}")
    print()

    # parse file
    for line in fh:
        elems = line.strip().split("\t")
        gene = elems[idx['Gene']]
        pheno1 = elems[idx[ds1col]]
        pheno2 = elems[idx[ds2col]]
        cs1nr = elems[idx[cs1col]]
        cs2nr = elems[idx[cs2col]]
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
        
        # GWAS variant chromosome
        # GWAS variant position
        gwasvarchridx = idx.get("GWAS variant chromosome")
        gwasvarposidx = idx.get("GWAS variant position")
        if gwasvarchridx is not None:
            colocObj.gwasvarchr = int(elems[gwasvarchridx])
            colocObj.gwasvarpos = int(elems[gwasvarposidx])

        # coloc is bidirectional
        phenObj1.colocresults.append(colocObj)
        phenObj2.colocresults.append(colocObj)
        csObj1.addColoc(colocObj)
        csObj2.addColoc(colocObj)
        colocresults.append(colocObj)

    print(f"{len(colocresults)} for {len(phenomap)} phenotypes.")
    return phenomap

def getEGeneListsFromOtherPubs(celltypes):
    

    matches = {
        "EX": "Excitatory.neurons",
        "IN": "Inhibitory.neurons",
        "AST": "Astrocytes",
        "PER": "Pericytes",
        "END": "Endothelial.cells",
        "OLI": "Oligodendrocytes",
        "OPC": "OPCs...COPs",
        "MIC": "Microglia"
    }
    
    signGenesBryois = {}
    eGenesBryois = set()
    for celltype in celltypes:
        matchingcelltype = matches.get(celltype)
        if matchingcelltype is not None:
            signGenesBryois[celltype] = set()
            topfxfile = f"/path/to/single_cell_eQTLsumstats/Bryois2022/{matchingcelltype}.top.gz"
            fh = gzip.open(topfxfile,'rt')
            header = fh.readline().strip().split("\t")
            idx = index(header)
            for line in fh:
                elems = line.strip().split("\t")
                gene = elems[idx['symbol_ensembl']].split("_")[0]
                adjp = float(elems[idx['adj_p']])
                if adjp < 0.05:
                    signGenesBryois[celltype].add(gene)
                    eGenesBryois.add(gene)
            fh.close()


    matches = {
        "EX": "Exc",
        "IN": "Inh",
        "AST": "Ast",
        "PER": "End",
        "END": "End",
        "OLI": "Oli",
        "OPC": "OPC",
        "MIC": "Mic"
    }

    signGenesFujita = {}
    eGenesFujita = set()
    for celltype in celltypes:
        matchingcelltype = matches.get(celltype)
        if matchingcelltype is not None:
            signGenesFujita[celltype] = set()
            topfxfile = f"/path/to/single_cell_eQTLsumstats/Fujita2022/AnalysisQTL/celltype/celltype-eqtl-sumstats.{matchingcelltype}topfx.txt.gz"
            fh = gzip.open(topfxfile,'rt')
            header = fh.readline().strip().split("\t")
            idx = index(header)
            for line in fh:
                elems = line.strip().split("\t")
                qval = float(elems[idx['qval']])    
                if qval < 0.05:
                    gene = elems[idx['gene_symbol']]
                    signGenesFujita[celltype].add(gene)
                    eGenesFujita.add(gene)
            fh.close()

    matches = {
        "OPC": "SingleBrain_OPC",
        "EX": "SingleBrain_Ext",
        "IN":"SingleBrain_IN",
        "AST":"SingleBrain_Ast",
        "OLI":"SingleBrain_OD",
        "MIC":"SingleBrain_MG",
        "END":"SingleBrain_End"
    }

    # load ensg to gene map
    ensgsymbolmap = {}
    fh = gzip.open("/path/to/single_cell_eQTLsumstats/SingleBrain2026/data/ensgtosymbol.txt.gz",'rt')
    for line in fh:
        elems = line.strip().split("\t")
        ensgsymbolmap[elems[0]] = elems[1]
    fh.close()

    signGenesSingleBrain = {}
    eGenesSingleBrain = set()
    for celltype in celltypes:
        matchingcelltype = matches.get(celltype)
        if matchingcelltype is not None:
            signGenesSingleBrain[celltype] = set()
            topfxfile = f"/path/to/single_cell_eQTLsumstats/SingleBrain2026/data/{matchingcelltype}_eqtl_top_assoc.tsv.gz-withFixedRandomFxZScoreAndFixedFxQval.tsv.gz"
            fh = gzip.open(topfxfile,'rt')
            header = fh.readline().strip().split("\t")
            idx = index(header)
            for line in fh:
                elems = line.strip().split("\t")
                qval = float(elems[idx['fixedfxqval']])    
                if qval < 0.05:
                    gene = elems[idx['feature']]
                    # look up symbol somehow
                    symbol = ensgsymbolmap.get(gene)
                    if symbol is not None:
                        signGenesSingleBrain[celltype].add(symbol)
                        eGenesSingleBrain.add(symbol)
            fh.close()
    return [signGenesBryois, eGenesBryois, signGenesFujita, eGenesFujita, signGenesSingleBrain, eGenesSingleBrain]