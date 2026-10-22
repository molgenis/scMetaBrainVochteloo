
from scipy.stats.distributions import norm
from py_qvalue import qvalue
import numpy as np
from scipy import stats
import Rb
import gzip
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
import traceback
import math as m
import scipy.stats as sps

def fopen(file,mode):
    if file.endswith(".gz"):
        return gzip.open(file,mode+"t")
    else:
        return open(file,mode)

def index(arr):
    
    idx = {}
    for i in range(0,len(arr)):
        idx[arr[i]] = i
    
    return idx
    
def pValToBetaAndSE(p, cor, maf, n):
    if p < 1e-200: # prevent clippinf og norm.ppf
        p = 1e-200
    z = norm.ppf(p)
    
    if cor > 0:
        z *=-1
    # print(f"{p}\t{cor}\t{maf}\{n}")
    # p2 = norm.sf(z)
    # print(f"p {p}\tz {z}\tp2 {1-p2}\tcorr {cor}")
    return zscoreToBetaAndSE(z, maf, n)

### ACAT function
def acat_test(pvalues,weights=None):
    '''acat_test()
    Aggregated Cauchy Assocaition Test
    A p-value combination method using the Cauchy distribution.

    Inspired by: https://github.com/yaowuliu/ACAT/blob/master/R/ACAT.R

    Author: Ryan Neff

    Inputs:
        pvalues: <list or numpy array>
            The p-values you want to combine.
        weights: <list or numpy array>, default=None
            The weights for each of the p-values. If None, equal weights are used.

    Returns:
        pval: <float>
            The ACAT combined p-value.
    '''
    if any(np.isnan(pvalues)):
        raise Exception("Cannot have NAs in the p-values.")
    if any([(i>1)|(i<0) for i in pvalues]):
        raise Exception("P-values must be between 0 and 1.")
    if any([i==1 for i in pvalues])&any([i==0 for i in pvalues]):
        raise Exception("Cannot have both 0 and 1 p-values.")
    if any([i==0 for i in pvalues]):
        print("Warn: p-values are exactly 0.")
        return 0
    if any([i==1 for i in pvalues]):
        print("Warn: p-values are exactly 1.")
        new_p = 1 - (1/len(pvalues)) # recommended by R version
        new_pvals = []
        for pval in pvalues:
                if pval > new_p:
                        new_pvals.append(new_p)
                else:
                        new_pvals.append(pval)
        pvalues = new_pvals
        #return 1
    if weights==None:
        weights = [1/len(pvalues) for i in pvalues]
    elif len(weights)!=len(pvalues):
        raise Exception("Length of weights and p-values differs.")
    elif any([i<0 for i in weights]):
        raise Exception("All weights must be positive.")
    else:
        weights = [i/len(weights) for i in weights]

    pvalues = np.array(pvalues)
    weights = np.array(weights)

    if any([i<1e-16 for i in pvalues])==False:
        cct_stat = sum(weights*np.tan((0.5-pvalues)*np.pi))
    else:
        is_small = [i<(1e-16) for i in pvalues]
        is_large = [i>=(1e-16) for i in pvalues]
        cct_stat = sum((weights[is_small]/pvalues[is_small])/np.pi)
        cct_stat += sum(weights[is_large]*np.tan((0.5-pvalues[is_large])*np.pi))

    if cct_stat>1e15:
        pval = (1/cct_stat)/np.pi
    else:
        pval = 1 - sps.cauchy.cdf(cct_stat)

    return pval

# See Zhang et al, Nature Genetics volume 48, pages481–487(2016) for details
# Formula taken from Supplementary Information
def zscoreToBetaAndSE(z, maf, n):
    if maf <= 0:
        return [0,0]
    if z == 0:
        return [0,0]
    
    chisq = z * z
    a = 2 * maf * (1-maf) * (n+chisq)
    # print(f"{z}\t{maf}\t{a}\t{m.sqrt(a)}\n")
    beta = z / m.sqrt(a)
    se = 1 / m.sqrt(a)
    return [beta, se]



def checkAlleles(alleles1,effectA1,alleles2,effectA2):

    refa = alleles1
    alta = effectA1
    refb = alleles2
    altb = effectA2
    # print()
    # print(f"{refa} - {alta}  /// {refb} - {altb}")

    if "/" in refa:
        refasplit = refa.split("/")
        refa = refasplit[0]
        if refa == alta:
            refa = refasplit[1]

    if "/" in refb:
        refbsplit = refb.split("/")
        refb = refbsplit[0]
        if refb == altb:
            refb = refbsplit[1]

    nrshared = 0
    ala = [refa,alta]
    alb = [refb,altb]
    for a in ala:
        for b in alb:
            if a == b:
                nrshared += 1
    # check whether there are at least two equal alleles
    if nrshared < 2:
        print("Not enough matching alleles")
        return None
    # print(f"{alta==altb}")
    return alta == altb


def compareDatasets(ct1,ct2,nrsignct1,b1,se1,b2,se2,theta,p2, z1,z2, outfilename):
    fh = fopen(outfilename,'w')

    fh.write("CT1\tCT2\tNrSignCT1\tNrOverlapCT2\tAC_NoFDRFilter\tRb\tRb_se\tRb_p\tCorrelZ\tNrSignificant\tPropSignificant\tAC_Significant\tPi1\tNote\n")
    
    if len(b1) > 0:
        # calculate Rb over all shared effects
        rb, rb_se, rb_p = Rb.calcu_cor_true(b1, se1, b2, se2, theta, skip_jacknife=False)
        note = "OK"
        # determine FDR over p2
        pi1 = 0
        q_values = None
        try:
            qvalue_results = qvalue(np.array(p2))
            # print(p2)
            q_values = qvalue_results['qvalues']
            # print(q_values)
            pi1 = 1-qvalue_results['pi0']
            # print(pi1)
        except:
            # for p in p2:
            #     print(p)
            # calculate FDR in stead
            print(traceback.format_exc())
            note = "QVALUEFAIL"
            q_values = stats.false_discovery_control(p2)
            pass
        
        ac = 0
        ac_sign = 0
        
        nrsign = 0
        propsign = 0
        for i in range(0,len(b1)):
            if b1[i] * b2[i] >= 0:
                ac+=1
                if q_values[i] < 0.05:
                    ac_sign +=1
            if q_values[i] < 0.05:        
                nrsign += 1
        ac /= len(b1)
        
        if nrsign > 0:
            ac_sign /= nrsign
        else:
            ac_sign = 0

        propsign = nrsign / len(b1)
        if ac > 1:
            ac = 1
        if ac < 0:
            ac = 0

        ac = round(ac,3)
        rb = round(rb,3)
        rb_se = round(rb_se,3)
        propsign = round(propsign,3)
        ac_sign = round(ac_sign,3)
        if pi1 > 1:
            pi1 = 1
        if pi1 < 0:
            pi1 = 0
        pi1 = round(pi1,3)
        
        if len(b1) < 250:
            if note != "OK":
                note += "LIMITEDOBSERVATIONS"
            else:
                note = "LIMITEDOBSERVATIONS"
       
       
        correlz = stats.pearsonr(z1,z2)
        correlz = correlz.statistic

        print(f"{ct1}\t{ct2}")
        print(f"Nr SignA: {nrsignct1}")
        print(f"Overlap: {len(b1)}")
        print(f"AC: {ac}")
        print(f"Rb: {rb}")
        print(f"Rb se: {rb_se}")
        print(f"Correl_Z: {correlz}")
        print(f"NrSign: {nrsign}")
        print(f"AC_sign: {ac_sign}")
        print(f"Pi1: {pi1}")
        print(f"Note: {note}")

        # 
        fh.write(f"{ct1}\t{ct2}\t{nrsignct1}\t{len(b1)}\t{ac}\t{rb}\t{rb_se}\t{rb_p}\t{correlz}\t{nrsign}\t{propsign}\t{ac_sign}\t{pi1}\t{note}\n")
    else:
        note = "NOOBSERVATIONS"
        ac = ""
        rb = ""
        rb_p = ""
        rb_se = ""
        nrsign = ""
        pi1 = ""
        fh.write(f"{ct1}\t{ct2}\t{nrsignct1}\t{len(b1)}\t{ac}\t{rb}\t{rb_se}\t{rb_p}\t{nrsign}\t{propsign}\t{ac_sign}\t{pi1}\t{note}\n")
    fh.close()

def compareDatasetsSingleBrain(ct1,ct2,nrsignct1,b1,se1,b2,se2,theta,p2, z1,z2, outfilename):
    fh = fopen(outfilename,'w')

    fh.write("CT1\tCT2\tNrSignCT1\tNrOverlapCT2\tAC_NoFDRFilter\tRb\tRb_se\tRb_p\tCorrelZ\tNrSignificant\tPropSignificant\tAC_Significant\tPi1\tNote\n")
    
    if len(b1) > 0:
        # calculate Rb over all shared effects
        rb, rb_se, rb_p = Rb.calcu_cor_true(b1, se1, b2, se2, theta, skip_jacknife=False)
        note = "OK"
        # determine FDR over p2
        pi1 = 0
        q_values = None
        print(p2)
        try:
            qvalue_results = qvalue(np.array(p2))
            # print(p2)
            q_values = qvalue_results['qvalues']
            print("")
            print("qvals:")
            print(qvalue_results)
            # print(q_values)
            pi1 = 1 - qvalue_results['pi0']
            # print(pi1)
        except:
            # for p in p2:
            #     print(p)
            # calculate FDR in stead
            print(traceback.format_exc())
            note = "QVALUEFAIL"
            q_values = stats.false_discovery_control(p2)
            pass
        
        print(f"pi1: {pi1}")
        ac = 0
        ac_sign = 0
        
        nrsign = 0
        propsign = 0
        for i in range(0,len(z1)):
            if z1[i] * z2[i] >= 0:
                ac+=1
                if q_values[i] < 0.05:
                    ac_sign +=1
            if q_values[i] < 0.05:        
                nrsign += 1
        ac /= len(b1)
        
        if nrsign > 0:
            ac_sign /= nrsign
        else:
            ac_sign = 0

        propsign = nrsign / len(b1)
        if ac > 1:
            ac = 1
        if ac < 0:
            ac = 0

        ac = round(ac,3)
        rb = round(rb,3)
        rb_se = round(rb_se,3)
        propsign = round(propsign,3)
        ac_sign = round(ac_sign,3)

        if pi1 > 1:
            pi1 = 1
        if pi1 < 0:
            pi1 = 0
        pi1 = round(pi1,3)
        
        if len(b1) < 250:
            if note != "OK":
                note += "LIMITEDOBSERVATIONS"
            else:
                note = "LIMITEDOBSERVATIONS"


        correlz = stats.pearsonr(z1,z2)
        correlz = correlz.statistic

        print(f"{ct1}\t{ct2}")
        print(f"Nr SignA: {nrsignct1}")
        print(f"Overlap: {len(b1)}")
        print(f"AC: {ac}")
        print(f"Rb: {rb}")
        print(f"Rb se: {rb_se}")
        print(f"Correl_Z: {correlz}")
        print(f"NrSign: {nrsign}")
        print(f"AC_sign: {ac_sign}")
        print(f"Pi1: {pi1}")
        print(f"Note: {note}")

        
        # 
        fh.write(f"{ct1}\t{ct2}\t{nrsignct1}\t{len(b1)}\t{ac}\t{rb}\t{rb_se}\t{rb_p}\t{correlz}\t{nrsign}\t{propsign}\t{ac_sign}\t{pi1}\t{note}\n")
    else:
        note = "NOOBSERVATIONS"
        ac = ""
        rb = ""
        rb_p = ""
        rb_se = ""
        nrsign = ""
        pi1 = ""
        fh.write(f"{ct1}\t{ct2}\t{nrsignct1}\t{len(b1)}\t{ac}\t{rb}\t{rb_se}\t{rb_p}\t{nrsign}\t{propsign}\t{ac_sign}\t{pi1}\t{note}\n")
    fh.close()

def makeComparisonsSingleBrain(outdir, ctA, ctB, ctASignPairs, overlapCtB, thetaVal, onlyCompareEffectAlleles=False, fixedFxMeta = False, useFlippedRandomZ=False):
    # determine allellic concordance, rb, etc

    fhoCompOut = fopen(outdir+ctA+"-vs-"+ctB+".txt",'w')
    compoutheader = "Gene\tSNP\talleles1\teffectallele1\talleles2\teffectallele2\tb1\tse1\tz1\tp1\tb2\tse2\tz2\tp2\tAllelesEqual"
    compoutheader += "\trandomz\tfixedz\trandomp\tfixedp"

    fhoCompOut.write(compoutheader+"\n")
    b1 = []
    se1 = []
    z1 = []
    b2 = []
    se2 = []
    z2 = []

    theta = []
    p1 = None
    p2 = []
    samedir = 0
    for pair in ctASignPairs.keys():
        ctAPair = ctASignPairs.get(pair)
        ctBPair = overlapCtB.get(pair)
        if ctBPair is not None:
            allelesOK = None
            if onlyCompareEffectAlleles:
                allelesOK = ctAPair["effectallele"] == ctBPair["effectallele"]
            else:
                allelesOK = checkAlleles(ctAPair["alleles"],ctAPair["effectallele"],ctBPair["alleles"],ctBPair["effectallele"])
            # print(f"{pair}")
            # print(f"{ctAPair}")
            # print(f"{ctBPair}")
            # print()
            if allelesOK is not None:
                b1est = ctAPair["beta"]
                b2est = ctBPair["beta"]
                se1est = ctAPair["se"]
                se2est = ctBPair["se"]

                # determine which dataset is singlebrain
                # "randomz": randomz,
                # "fixedz": fixedz,
                # "randomp": randomp,
                # "fixedp": fixedp
                

                include = True

                if b1est is None or b2est is None or se1est is None or se2est is None:
                    include = False
                elif m.isnan(b1est) or m.isnan(b2est) or m.isnan(se1est) or m.isnan(se2est):
                    include = False
                elif m.isinf(b1est) or m.isinf(b2est) or m.isinf(se1est) or m.isinf(se2est):
                    include = False
                
                if include:
                    z1est = None
                    z2est = None
                    p1est = None
                    p2est = None

                    if (ctAPair.get('fixedz') is not None):
                        # this dataset is SingleBrain
                        if fixedFxMeta:
                            z1est = ctAPair['fixedz']
                            p1est = ctAPair['fixedp']
                        else:
                            z1est = ctAPair['randomz']
                            if useFlippedRandomZ:
                                z1est = ctAPair['randomzflipped']
                            p1est = ctAPair['randomp']
                        p2est = ctBPair['p']
                        z2est = ctBPair['beta'] / ctBPair['se']
                    else:
                        if fixedFxMeta:
                            z2est = ctBPair['fixedz']
                            p2est = ctBPair['fixedp']
                        else:
                            z2est = ctBPair['randomz']
                            if useFlippedRandomZ:
                                z2est = ctBPair['randomzflipped']
                            p2est = ctBPair['randomp']
                        p1est = ctAPair['p']
                        z1est = ctAPair['beta'] / ctAPair['se']

                    b1.append(ctAPair["beta"])
                    se1.append(ctAPair["se"])
                    z1.append(z1est)  

                    if not allelesOK:
                        b2est = -b2est
                        z2est = -z2est
                    b2.append(b2est)
                    se2.append(se2est)
                    z2.append(z2est)  
                    p2.append(p2est)

                    alleles1 = ctAPair["alleles"]
                    effectallele1 = ctAPair["effectallele"]
                    alleles2 = ctBPair["alleles"]
                    effectallele2 = ctBPair["effectallele"]

                    # n.append(min(ctAPair[-1],ctBPair[-1]))
                    theta.append(0) # sample overlap between datasets
                    snp = ctAPair["snp"]
                    gene = ctAPair["gene"]
                    ol = f"{gene}\t{snp}\t{alleles1}\t{effectallele1}\t{alleles2}\t{effectallele2}\t{b1est}\t{se1est}\t{z1est}\t{[p1est]}\t{b2est}\t{se2est}\t{z2est}\t{p2est}\t{allelesOK}\n"
                    
                    randomz = ctAPair.get("randomz")
                    fixedz = ctAPair.get("fixedz")
                    randomp = ctAPair.get("randomp")
                    fixedp = ctAPair.get("fixedp")
                    if fixedz is None:
                        randomz = ctBPair.get("randomz")
                        fixedz = ctBPair.get("fixedz")
                        randomp = ctBPair.get("randomp")
                        fixedp = ctBPair.get("fixedp")
                        
                        ol += f"\t{randomz}\t{fixedz}\t{randomp}\t{fixedp}\n"
                    fhoCompOut.write(ol)                                
    fhoCompOut.close()

    # Shared.compareDatasets(ctA,ctB, b1,se1,b2,se2,n,p2, )
    outfilename = outdir+ctA+"-vs-"+ctB+"-summary.txt"
    # print(p2)
    compareDatasetsSingleBrain(ctA,ctB,len(ctASignPairs),b1,se1,b2,se2,theta,p2, z1,z2, outfilename)

def makeComparisons(outdir, ctA, ctB, ctASignPairs, overlapCtB, thetaVal, onlyCompareEffectAlleles=False, addSingleBraincolumns=False):
    # determine allellic concordance, rb, etc

    fhoCompOut = fopen(outdir+ctA+"-vs-"+ctB+".txt",'w')
    compoutheader = "Gene\tSNP\talleles1\teffectallele1\talleles2\teffectallele2\tb1\tse1\tp1\tb2\tse2\tp2\tz1\tz2\tAllelesEqual"
    if addSingleBraincolumns:
        compoutheader += "\trandomz\tfixedz\trandomp\tfixedp"
    fhoCompOut.write(compoutheader+"\n")
    b1 = []
    se1 = []
    b2 = []
    se2 = []
    theta = []
    z1 = []
    z2 = []
    p1 = None
    p2 = []
    samedir = 0
    for pair in ctASignPairs.keys():
        ctAPair = ctASignPairs.get(pair)
        ctBPair = overlapCtB.get(pair)
        if ctBPair is not None:
            allelesOK = None
            if onlyCompareEffectAlleles:
                allelesOK = ctAPair["effectallele"] == ctBPair["effectallele"]
            else:
                allelesOK = checkAlleles(ctAPair["alleles"],ctAPair["effectallele"],ctBPair["alleles"],ctBPair["effectallele"])
            # print(f"{pair}")
            # print(f"{ctAPair}")
            # print(f"{ctBPair}")
            # print()
            if allelesOK is not None:
                b1est = ctAPair["beta"]
                b2est = ctBPair["beta"]
                se1est = ctAPair["se"]
                se2est = ctBPair["se"]
                include = True

                if b1est is None or b2est is None or se1est is None or se2est is None:
                    include = False
                elif m.isnan(b1est) or m.isnan(b2est) or m.isnan(se1est) or m.isnan(se2est):
                    include = False
                elif m.isinf(b1est) or m.isinf(b2est) or m.isinf(se1est) or m.isinf(se2est):
                    include = False
                
                if include:
                    b1.append(b1est)
                    se1.append(se1est)
                    z1est = b1est / se1est
                    if not allelesOK:
                        b2est = -b2est
                    z2est = b2est / se2est
                    z1.append(z1est)
                    z2.append(z2est)
                    b2.append(b2est)
                    se2.append(se2est)
                    p2.append(ctBPair["p"])
                    alleles1 = ctAPair["alleles"]
                    effectallele1 = ctAPair["effectallele"]
                    alleles2 = ctBPair["alleles"]
                    effectallele2 = ctBPair["effectallele"]

                    # n.append(min(ctAPair[-1],ctBPair[-1]))
                    theta.append(thetaVal) # sample overlap between datasets
                    snp = ctAPair["snp"]
                    gene = ctAPair["gene"]
                    ol = f"{gene}\t{snp}\t{alleles1}\t{effectallele1}\t{alleles2}\t{effectallele2}\t{b1est}\t{se1est}\t{ctAPair["p"]}\t{b2est}\t{se2est}\t{ctBPair["p"]}\t{z1est}\t{z2est}\t{allelesOK}\n"
                    if addSingleBraincolumns:
                        randomz = ctAPair.get("randomz")
                        fixedz = ctAPair.get("fixedz")
                        randomp = ctAPair.get("randomp")
                        fixedp = ctAPair.get("fixedp")
                        if fixedz is None:
                            randomz = ctBPair.get("randomz")
                            fixedz = ctBPair.get("fixedz")
                            randomp = ctBPair.get("randomp")
                            fixedp = ctBPair.get("fixedp")
                            
                        ol += f"\t{randomz}\t{fixedz}\t{randomp}\t{fixedp}\n"
                    fhoCompOut.write(ol)                                
    fhoCompOut.close()

    # Shared.compareDatasets(ctA,ctB, b1,se1,b2,se2,n,p2, )
    outfilename = outdir+ctA+"-vs-"+ctB+"-summary.txt"
    compareDatasets(ctA,ctB,len(ctASignPairs),b1,se1,b2,se2,theta,p2, z1,z2, outfilename)

def getLimixEffects(file, filtersignificance, variantgenepairs, genelimitset, variantEncoding=None, geneEncoding=None, printErrors=False, alleleFile = None, allelesInSnpId=False, recalculateP=False):
    print("Parsing "+file)
    fh = fopen(file,'r')
    header = fh.readline().rstrip().split("\t")

    header = index(header)
    snpcol = header["snp_id"]
    snpposcol = header["snp_position"]
    snpchrcol = header["snp_chromosome"]
    genecol = header["feature_id"]
    betacol = header["beta"]
    secol = header["beta_se"]
    qvalcol = header.get("qval")
    ncol = header["n_samples"]
    effalcol = header["assessed_allele"]
    pcol = header["p_value"]
    ensgcol = header.get("ensembl_id")
    if ensgcol is None:
        ensgcol = header.get("ENSG")

    output = {}

    if filtersignificance and not type(qvalcol) == int:
        print("Error in limix file: "+file)
        print("qvalue column is not found")
        exit()
    else:
        print(f"Filtering for significance: {filtersignificance}")
        print(f"{qvalcol} qval col")
    
    testedgenes = set()
    genes = set()
    lctr = 0
    wctr = 0

    rsids = set()
    for line in fh:
        try:
            elems = line.rstrip().split("\t")
            snp = elems[snpcol]
            if variantEncoding is not None:
                if variantEncoding == "CIMA" or variantEncoding == "GTEX":
                    snpchr = elems[snpchrcol]
                    snppos = elems[snpposcol]
                    snp = "chr"+snpchr+"_"+snppos

            gene = elems[genecol]
            ensg = elems[ensgcol]
            pair = gene+"_"+snp
            if geneEncoding == "ENSG":
                pair = ensg+"_"+snp
                gene = ensg
            if genelimitset is None or gene in genelimitset:
                rsids.add(elems[snpcol])
  
                included = True
                
                # print(line)
                if filtersignificance:
                    q = float(elems[qvalcol])
                    if q > 0.05:
                        included = False
                    # print(f"{included} - {included}")
                elif variantgenepairs is not None:
                    if pair not in variantgenepairs:
                        included = False
                # print(included)
                if included:    
                    beta = float(elems[betacol])
                    se = float(elems[secol])
                    p = float(elems[pcol])
                    alleles = "ref"
                    if allelesInSnpId:
                        idelems = elems[snpcol].split(":")
                        alleles = idelems[2]+"/"+idelems[3]
                    effallele = elems[effalcol]
                    samplesize = int(float(elems[ncol]))

                    # limix gives a p-value of 0.9999 to variants that are not significant in the odds-ratio test (i.e. the genotypes are collinear with any covariates or the random effects matrix)
                    # this throws off the pi1 / q-value calculations because it causes a spike of tests
                    # with a p-value near 1. The code below recalculates the p-values for such variants using the 
                    # beta and sd that were extracted from the full model with a term for the variant
                    # the recalculated p-value will likely be very insignificant, and this fixes the issue with q-values.
                    # this is mainly an issue in smaller datasets, lower expressed genes and cell with low cell type counts
                    if recalculateP:
                        if p > 0.9:
                            z = abs(beta / se)
                            p = sps.norm.sf(abs(z))*2

                    output[pair] = {
                        "snp": snp,
                        "rsid": elems[snpcol],
                        "gene": gene,
                        "beta": beta,
                        "se": se,
                        "alleles": alleles,
                        "effectallele": effallele,
                        "samplesize": samplesize,
                        "p": p
                    }
                    genes.add(gene)
                    wctr += 1
                       
                lctr += 1
                if lctr % 100000 == 0:
                    print(f"{lctr} lines parsed, {len(output)} kept",end='\r')
                    
        except Exception as e:
            if printErrors:
                print(Exception, e)
                print("Error with line:")
                print(line)

            # print(f"{snpcol}	{snpposcol}	{snpchrcol}	{genecol}	{betacol}	{secol}	{qvalcol}	{alcol}	{ncol}	{effalcol}	{pcol}	{ensgcol}")
            # print(f"{len(elems)} elements")
            # print()
            # print(traceback.format_exc())
            # exit()
    # load in allele info if available
    print(f"{lctr} lines parsed, {len(output)} kept",end='\n')

    if alleleFile is not None:
        print(f"Loading allele information from: "+alleleFile)

        rsidToData = {}
        for pair in output.keys():
            data = output[pair]
            rsid = data["rsid"]
            rsiddata = rsidToData.get(rsid)
            if rsiddata is None:
                rsiddata = []
            rsiddata.append(data)
            rsidToData[rsid] = rsiddata

        replaced = 0
        fh = fopen(alleleFile,'r')
        for line in fh:
            if not line.startswith("#"):
                elems = line.strip().split("\t")
                rsid = elems[2]
                rsiddata = rsidToData.get(rsid)
                if rsiddata is not None:
                    alleles = elems[3]+"/"+elems[4]
                    for dc in rsiddata:
                        # check if eff allele is a1 or a2...
                        if dc["effectallele"] in alleles:
                            dc["alleles"] = alleles
                            replaced += 1
        fh.close()
        print(f"{replaced} / {len(output)} have replaced alleles")
   
    return [genes, output]

def findQTLsInscCIMAFile(cimafile,pairs):

    print("Parsing: "+cimafile)
    fh = fopen(cimafile,'r')
    header = fh.readline().rstrip().split("\t")
    
    
    header = index(header)
    variantcol = header["variant_id"]
    slopecol = header["slope"]
    secol = header["slope_se"]
    refalcol = header["A2(REF)"]
    effalcol = header["A1(ALT/effect allele)"]
    samplesizecol = header["sample_sizes"]
    pvalcol = header["pval_nominal"]    
    genecol = header["phenotype_id"]
    ctcol = header.get("celltype")
    analysiscol = header.get("analysis")

    output = {}
    lnctr = 0
    for line in fh:
        elems = line.rstrip().split("\t")
        try:
            gene = elems[genecol]
            snp = elems[variantcol]
            beta = float(elems[slopecol])
            se = float(elems[secol])
            refal = elems[refalcol]
            effallele = elems[effalcol]
            alleles = refal+"/"+effallele
            pair = gene +"_"+snp
            samplesize = int(elems[samplesizecol])
            p = float(elems[pvalcol])
            if pairs is None or pair in pairs:
                output[pair] = {
                    "snp": snp,
                    "gene": gene,
                    "beta": beta,
                    "se": se,
                    "alleles": alleles,
                    "effectallele": effallele,
                    "samplesize": samplesize,
                    "p": p
                }

            lnctr += 1
        except Exception as e:
            
            print(traceback.format_exc())
            print(Exception, e)
            print("Error with line:")
            print(line)
            # exit()
        if lnctr % 100000 == 0:
            print(f"{lnctr} lines parsed, {len(output)} / {len(pairs)} found",end='\r')
    print(f"{lnctr} lines parsed, {len(output)} / {len(pairs)} found",end='\n')
    return output

def findQTLsInGTExParquet(parquetfile, pairs):
    print("Reading: "+parquetfile)

    table = pq.read_table(parquetfile)
    # print(table.column_names)
    output = {}
    reader = table.to_reader()
    rctr = 0
    mctr = 0
    for recordbatch in reader:
        for d in recordbatch.to_pylist():
            ensg = d['gene_id'].split(".")[0]
            variant = d['variant_id'].split("_") # chr22_15603514_C_T_b38
            snp = variant[0]+"_"+variant[1]
            ref = variant[2]
            alt = variant[3]
            p = d['pval_nominal']
            beta = d['slope']
            se = d['slope_se']
            pair = ensg+"_"+snp
            if pairs is None or pair in pairs:
                alleles = ref+"/"+alt
                output[pair] = {
                    "snp": snp,
                    "gene": ensg,
                    "beta": beta,
                    "se": se,
                    "alleles": alleles,
                    "effectallele": alt,
                    "samplesize": 0,
                    "p": p
                }

            rctr +=1
            if rctr % 100000 == 0:
                print(f"{rctr} lines read, {len(output)} / {len(pairs)} matched",end='\r')
    print(f"{rctr} lines read, {len(output)} / {len(pairs)} matched",end='\n')
    return output

def readGeneAnnotation(file):
    output = {}
    fh = fopen(file,'r')
    for line in fh:
        elems = line.strip().split("\t")
        ensg = elems[1].split(".")[0]
        symb = elems[2]
        output[symb] = ensg
    return output

def getGTExTopFx(qtlfile):
    fh = fopen(qtlfile,'r')
    header = fh.readline().rstrip().split("\t")
    snpcol = None
    snpposcol = None
    snpchrcol = None
    genecol = None
    betacol = None
    secol= None
    ncol = None
    qvalcol = None
    alcol = None
    effalcol = None

    for i in range(0,len(header)):
        if header[i] == "gene_id":
            genecol = i
        if header[i] == "variant_id":
            snpcol = i
        if header[i] == "slope":
            betacol = i
        if header[i] == "slope_se":
            secol = i
        if header[i] == "qval":
            qvalcol = i
        
    output = {}
    for line in fh:
        elems = line.rstrip().split("\t")
        ensg = elems[genecol].split(".")[0]
        variant = elems[snpcol].split("_") # chr22_15603514_C_T_b38
        snp = variant[0]+"_"+variant[1]
        ref = variant[2]
        alt = variant[3]
        beta = float(elems[betacol])
        se = float(elems[secol])
        pair = ensg+"_"+snp
        alleles = ref+"/"+alt
        qval = float(elems[qvalcol])
        if qval <= 0.05:
            output[pair] = {
                    "snp": snp,
                    "gene": ensg,
                    "beta": beta,
                    "se": se,
                    "alleles": alleles,
                    "effectallele": alt,
                    "samplesize": 0,
                    "p": 1
                }
    print(f"{len(output)} combos in {qtlfile}")
    return output

def getEffectsMbQTL(file, pairs):
    fh = fopen(file,'r')
    header = fh.readline().rstrip().split("\t")
    header = index(header)

    snpcol = header["SNP"]
    snpposcol = header["SNPPos"]
    snpchrcol = header["SNPChr"]
    genecol = header["GeneSymbol"]
    betacol = header["MetaBeta"]
    secol = header["MetaSE"]
    ncol = header["MetaPN"]
    qvalcol = header.get("qval")
    if qvalcol is None:
        qvalcol = header.get("Sign") 
    alcol = header["SNPAlleles"]
    effalcol = header["SNPEffectAllele"]
    pcol = header["MetaP"]

    # Gene    GeneChr GenePos GeneStrand      GeneSymbol      SNP     SNPChr  SNPPos  SNPAlleles      SNPEffectAllele SNPEffectAlleleFreq     QTLType MetaP   MetaPN  MetaPZ  MetaBeta        MetaSE  MetaI2  NrDatasets      DatasetCorrelationCoefficients(Dataset)    DatasetZScores(Dataset) DatasetSampleSizes(Dataset)

    output = {}
    lctr = 0
    mctr = 0
    fname = file.split("/")[-1]
    for line in fh:
        elems = line.rstrip().split("\t")
        gene = elems[genecol]
        snp = elems[snpcol]
        snpelems = snp.split(":")
        if len(snpelems) > 2:
            snp = snpelems[2]
        snpchr = elems[snpchrcol]
        snppos = elems[snpposcol]
        alleles = elems[alcol]
        alleleElems = alleles.split("/")
        # snp = snpchr+":"+snppos+""
        pair = gene+"_"+snp
        # print(pair)
        if pairs is None or pair in pairs:
            beta = float(elems[betacol])
            se = float(elems[secol])
            samplesize = int(elems[ncol])
            alt = elems[effalcol]
            p = float(elems[pcol])
            if qvalcol != None and pairs is None:
                qval = float(elems[qvalcol])
                if qval<=0.05:
                    output[pair] = {
                    "snp": snp,
                    "gene": gene,
                    "beta": beta,
                    "se": se,
                    "alleles": alleles,
                    "effectallele": alt,
                    "samplesize": samplesize,
                    "p": p
                }
                mctr += 1
            else:
                output[pair] = {
                    "snp": snp,
                    "gene": gene,
                    "beta": beta,
                    "se": se,
                    "alleles": alleles,
                    "effectallele": alt,
                    "samplesize": samplesize,
                    "p": p
                }
                mctr += 1
        lctr += 1
        if lctr % 100000 == 0:
           print(f"{lctr} lines parsed {mctr} matches - {fname}",end='\r')
    print(f"{lctr} lines parsed {mctr} matches - {fname}",end='\n') 
    return output

def findQTLsInOneK1KFile(file, chrstr, pairs):

    fh = fopen(file,'r')
    header = fh.readline().rstrip().split("\t")
    header = index(header)
    snpcol = header["RSID"]
    snpposcol = header["POS"]
    snpchrcol = header["CHR"]
    genecol = header["GENE"]
    corcol = header["SPEARMANS_RHO"]
    pvalcol = header["P_VALUE"]
    afcol = header["A2_FREQ_ONEK1K"]
    qvalcol = header["Q_VALUE"]
    alcol = header["A1"]
    effalcol = header["A2"]

# CELL_ID CELL_TYPE       RSID    SNPID   GENE    GENE_ID CHR     POS     A1      A2      A2_FREQ_ONEK1K  A2_FREQ_HRC     SPEARMANS_RHO   S_STATISTICS    P_VALUE Q_VALUE FDR     RSQUARE GENOTYPED       ROUND
    output = {}
    lctr = 0
    mctr = 0
    for line in fh:
        elems = line.rstrip().split("\t")
        gene = elems[genecol]
        snp = elems[snpcol]
        snpchr = elems[snpchrcol]
        snppos = elems[snpposcol]

        pair = gene+"_"+snp
        if pairs is None or pair in pairs:
            corr = float(elems[corcol])
            p = float(elems[pvalcol])
            n = 982
            maf = float(elems[afcol])
            conv = Rb.pValToBetaAndSE(p,corr,maf,n) # convert p to z to beta and SE!
            if conv is not None:
                beta = conv[0]
                se = conv[1]
                alleles = elems[alcol]
                effallele = elems[effalcol]
                alleles = alleles+"/"+effallele
                
                output[pair] = {
                    "snp": snp,
                    "gene": gene,
                    "beta": beta,
                    "se": se,
                    "alleles": alleles,
                    "effectallele": effallele,
                    "samplesize": n,
                    "p": p
                }

                mctr += 1
        lctr += 1
        if lctr % 100000 == 0:
            print(f"{lctr} lines parsed {mctr} matches - {file}",end='\r')
    print(f"{lctr} lines parsed {mctr} matches  - {file}",end='\n') 
    return output

def findQTLsInSingleBrain(file, pairs, useFixedFxFDR=False, genelimitset=None):
# feature
# variant_id
# chr
# pos
# ref
# alt
# Allele
# beta_tissue_0
# sd_tissue_0
# z_tissue_0
# beta_tissue_1
# sd_tissue_1
# z_tissue_1
# beta_tissue_2
# sd_tissue_2
# z_tissue_2
# beta_tissue_3
# sd_tissue_3
# z_tissue_3
# fixed_beta
# fixed_sd
# fixed_z
# Random_Z
# Fixed_P
# Random_P
# Fixed_bonf
# Random_bonf
# Fixed_FDR
# Random_FDR
# qval
# Random_Z_flipped
# fixedfxqval
    print("Parsing: "+file)
    fh = fopen(file,'r')
    header = fh.readline().rstrip().split("\t")
    header = index(header)

    genecol = header["feature"]
    variantcol = header["variant_id"]
    slopecol = header["fixed_beta"]
    secol = header["fixed_sd"]
    refalcol = header["ref"]
    effalcol = header["Allele"]
    
    randomzcol = header["Random_Z"]
    fixedzcol = header["fixed_z"]

    randompcol = header["Random_P"]
    fixedpcol = header["Fixed_P"]

    randomqcol = header.get("qval")
    fixedqcol = header.get("fixedfxqval")

    output = {}
    lnctr = 0
    for line in fh:
        elems = line.rstrip().split("\t")
        try:
            gene = elems[genecol].split(".")[0]
            if genelimitset is None or gene in genelimitset:
                snp = elems[variantcol]
                beta = float(elems[slopecol])
                if beta != 0.0: # skip weird results with exact zero estimate
                    se = float(elems[secol])
                    refal = elems[refalcol]
                    effallele = elems[effalcol]
                    alleles = refal+"/"+effallele
                    pair = gene +"_"+snp

                    randomp = None
                    randomz = None
                    randomzflipped = None

                    fixedz = float(elems[fixedzcol])
                    fixedp = float(elems[fixedpcol])

                    try:
                        randomp = float(elems[randompcol])
                        randomz = float(elems[randomzcol])
                        randomzflipped = randomz
                        if randomz * fixedz < 0:
                            randomzflipped *= -1
                    except:
                        pass

                    if pairs is None or pair in pairs:
                        include = False
                        if pairs is None:
                            randomqval = float(elems[randomqcol])
                            fixedqval = float(elems[fixedqcol])

                            if (useFixedFxFDR and fixedqval <= 0.05) or (not useFixedFxFDR and randomqval <= 0.05):
                                include = True
                        elif pair in pairs:
                            include = True

                        if include:

                            if randomp is not None and fixedp is not None:
                                if not np.isnan(randomp) and not np.isnan(randomp):
                                    output[pair] = {
                                        "snp": snp,
                                        "gene": gene,
                                        "beta": beta,
                                        "se": se,
                                        "alleles": alleles,
                                        "effectallele": effallele,
                                        
                                        "randomz": randomz,
                                        "randomzflipped": randomzflipped,
                                        "fixedz": fixedz,
                                        "randomp": randomp,
                                        "fixedp": fixedp,
                                    }

            lnctr += 1
        except Exception as e:
            
            print(traceback.format_exc())
            print(Exception, e)
            print("Error with line:")
            print(line)
            # exit()
        if lnctr % 100000 == 0:
            if pairs is not None:
                print(f"{lnctr} lines parsed, {len(output)} / {len(pairs)} found",end='\r')
            else:
                print(f"{lnctr} lines parsed, {len(output)} pairs",end='\r')
    if pairs is not None:
        print(f"{lnctr} lines parsed, {len(output)} / {len(pairs)} found",end='\n')
    else:
        print(f"{lnctr} lines parsed, {len(output)} pairs.",end='\n')
    print(f"{len(output)} loaded from singlebrain" )
    return output    

def findQTLsInFujita(fujitafile,pairs):

    print("Parsing: "+fujitafile)
    fh = fopen(fujitafile,'r')
    header = fh.readline().rstrip().split("\t")
    
    
    header = index(header)
    variantcol = header["snps"]
    slopecol = header["beta"]
    secol = header["se"]
    refalcol = header["REF"]
    effalcol = header["ALT"]
    
    pvalcol = header["pvalue"]    
    genecol = header["gene_symbol"]
    significancecol = header['significant_by_2step_FDR']
    
    output = {}
    lnctr = 0
    for line in fh:
        elems = line.rstrip().split("\t")
        try:
            gene = elems[genecol]
            snp = elems[variantcol]
            beta = float(elems[slopecol])
            se = float(elems[secol])
            refal = elems[refalcol]
            effallele = elems[effalcol]
            alleles = refal+"/"+effallele
            pair = gene +"_"+snp


            significant = elems[significancecol]
            if significant == "No":
                significant = False
            elif significant == "Yes":
                significant = True
            else:
                significant = False

            
            p = float(elems[pvalcol])
            if (pairs is None and significant) or (pairs is not None and pair in pairs):
                output[pair] = {
                    "snp": snp,
                    "gene": gene,
                    "beta": beta,
                    "se": se,
                    "alleles": alleles,
                    "effectallele": effallele,
                    "samplesize": None,
                    "p": p
                }

            lnctr += 1
        except Exception as e:
            
            print(traceback.format_exc())
            print(Exception, e)
            print("Error with line:")
            print(line)
            # exit()
        if lnctr % 100000 == 0:
            print(f"{lnctr} lines parsed, {len(output)} kept",end='\r')
    print(f"{lnctr} lines parsed, {len(output)} kept",end='\n')
    return output

def findQTLsInBryois(bryoisfile,pairs):

    print("Parsing: "+bryoisfile)
    fh = fopen(bryoisfile,'r')
    
    variantcol = 1
    slopecol = 4
    pvalcol = 3
    genecol = 0

    output = {}
    lnctr = 0
    for line in fh:
        elems = line.rstrip().split(" ")
        try:
            gene = elems[genecol].split("_")[-1]
            snp = elems[variantcol]
            beta = float(elems[slopecol])
            
            pair = gene +"_"+snp
            # print(pair)
            
            p = float(elems[pvalcol])

            
            if pairs is None or pair in pairs:
                # convert to z-score
                if p < 1e-200: # prevent clippinf og norm.ppf
                    p = 1e-200
                z = norm.ppf(p)
                se = abs(beta) / abs(z)

                output[pair] = {
                    "snp": snp,
                    "rsid": snp,
                    "gene": gene,
                    "beta": beta,
                    "se": se,
                    "samplesize": None,
                    "p": p
                }

            lnctr += 1
        except Exception as e:
            
            print(traceback.format_exc())
            print(Exception, e)
            print("Error with line:")
            print(line)
            # exit()
        if lnctr % 100000 == 0:
            print(f"{lnctr} lines parsed, {len(output)} / {len(pairs)} found",end='\r')
    print(f"{lnctr} lines parsed, {len(output)} / {len(pairs)} found",end='\n')
    return output

def findQTLsInBryoisTopFx(file):
    print("Parsing: "+file)
    fh = fopen(file,'r')
    output = {}

    header = fh.readline().split("\t")
    # symbol_ensembl	SNP			dist_TSS	beta	pval	bpval	adj_p	beta_metabrain	p_metabrain	Replication
    hidx = index(header)
    genecol = hidx['symbol_ensembl']
    snpcol = hidx['SNP']
    effalcol = hidx['effect_allele']
    refalcol = hidx['other_allele']
    betacol = hidx['beta']
    
    pcol = hidx['pval']
    adjpcol = hidx['adj_p']


    hidx = index(header)
    for line in fh:
        elems = line.strip().split("\t")

        adjp = float(elems[adjpcol])
        if adjp < 0.05:
            gene = elems[genecol]
            snp = elems[snpcol]
            beta = float(elems[betacol])

            refallele = elems[refalcol]
            effallele = elems[effalcol]
            alleles = refallele+"/"+effallele
            p = float(elems[pcol])

            # convert to z-score
            if p < 1e-200: # prevent clippinf og norm.ppf
                p = 1e-200
            z = norm.ppf(p)
            se = abs(beta) / abs(z)

            gene = gene.split("_")[-1]
            pair = gene +"_"+snp

            print(pair)
            
            output[pair] = {
                    "snp": snp,
                    "gene": gene,
                    "beta": beta,
                    "se": se,
                    "alleles": alleles,
                    "effectallele": effallele,
                    "samplesize": None,
                    "p": p
                }
    fh.close()
    print(f"{len(output)} significant")
    return output
