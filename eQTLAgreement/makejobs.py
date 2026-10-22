import os


curdir="/path/to/output/"

fh = open("./tpl.tpl",'r')
tpl = fh.readlines()
fh.close()

########### scmetabrain
logdir=curdir+"/logs/scmetabrain/"
jobdir=curdir+"/jobs/scmetabrain/"
outdir=curdir+"/output/scmetabrain/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

celltypes = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]

for ctA in celltypes:
    for ctB in celltypes:
        fho = open(jobdir+ctA+"-"+ctB+".sh",'w')
        for line in tpl:
            ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
            ol = ol.replace("SCRIPT","compare_celltypes_within_scMetaBrain.py")
            ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
            ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
            ol = ol.replace("CMD",f"{ctA} {ctB}")
            ol = ol.replace("OUTDIR",outdir)
            fho.write(ol)
        fho.close()

logdir=curdir+"/logs/scmetabrain-subset/"
jobdir=curdir+"/jobs/scmetabrain-subset/"
outdir=curdir+"/output/scmetabrain-subset/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

celltypes = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
genesubset = "/path/to/meanexp-merged-intersect-selection0.1.txt"
for ctA in celltypes:
    for ctB in celltypes:
        fho = open(jobdir+ctA+"-"+ctB+".sh",'w')
        for line in tpl:
            ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
            ol = ol.replace("SCRIPT","compare_celltypes_within_scMetaBrain.py")
            ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
            ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
            ol = ol.replace("CMD",f"{ctA} {ctB}")
            ol = ol.replace("OUTDIR",outdir+" "+genesubset)
            fho.write(ol)
        fho.close()

# subsampled
logdir=curdir+"/logs/scmetabrain-downsampled/"
jobdir=curdir+"/jobs/scmetabrain-downsampled/"
outdir=curdir+"/output/scmetabrain-downsampled/"

os.makedirs(jobdir, exist_ok=True)
os.makedirs(logdir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)

celltypes = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]

for ctA in celltypes:
    for ctB in celltypes:
        fho = open(jobdir+ctA+"-"+ctB+".sh",'w')
        for line in tpl:
            ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
            ol = ol.replace("SCRIPT","compare_celltypes_within_scMetaBrain_subsampled.py")
            ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
            ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
            ol = ol.replace("CMD",f"{ctA} {ctB}")
            ol = ol.replace("OUTDIR",outdir)
            fho.write(ol)
        fho.close()

logdir=curdir+"/logs/scmetabrain-downsampled-subset/"
jobdir=curdir+"/jobs/scmetabrain-downsampled-subset/"
outdir=curdir+"/output/scmetabrain-downsampled-subset/"

os.makedirs(jobdir, exist_ok=True)
os.makedirs(logdir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)
genesubset = "/path/to/meanexp-downsampled-merged-intersect-selection0.1.txt"
for ctA in celltypes:
    for ctB in celltypes:
        fho = open(jobdir+ctA+"-"+ctB+".sh",'w')
        for line in tpl:
            ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
            ol = ol.replace("SCRIPT","compare_celltypes_within_scMetaBrain_subsampled.py")
            ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
            ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
            ol = ol.replace("CMD",f"{ctA} {ctB}")
            ol = ol.replace("OUTDIR",outdir+" "+genesubset)
            fho.write(ol)
        fho.close()

# ########## MetaBrain
logdir=curdir+"/logs/metabrain/"
jobdir=curdir+"/jobs/metabrain/"
outdir=curdir+"/output/metabrain/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

celltypesA = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
celltypesB = ["basalganglia-EUR","cerebellum-EUR","cortex-EUR","cortex-EUR-rerun","hippocampus-EUR","spinalcord-EUR"]
for ctA in celltypesA:
    for ctB in celltypesB:
        if ctA != ctB:
            fho = open(jobdir+ctA+"-"+ctB+".sh",'w')
            for line in tpl:
                ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
                ol = ol.replace("SCRIPT","compare_MetaBrain_with_scMetaBrain.py")
                ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
                ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
                ol = ol.replace("CMD",f"{ctA} {ctB}")
                ol = ol.replace("OUTDIR",outdir)
                fho.write(ol)
            fho.close()
      

# ######### multiome
logdir=curdir+"/logs/multiome/"
jobdir=curdir+"/jobs/multiome/"
outdir=curdir+"/output/multiome/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

celltypesA = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
celltypesB = ["B","CD8T","CD4T","DC","NK","monocyte"]
for ctA in celltypesA:
    for ctB in celltypesB:
        if ctA != ctB:
            fho = open(jobdir+ctA+"-"+ctB+".sh",'w')
            for line in tpl:
                ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
                ol = ol.replace("SCRIPT","compare_multiomics_with_scMetaBrain.py")
                ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
                ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
                ol = ol.replace("CMD",f"{ctA} {ctB}")
                ol = ol.replace("OUTDIR",outdir)
                fho.write(ol)
            fho.close()

############## CIMA
# this analysis is unidirectional (not replicating CIMA eQTLs in scMetaBrain)
celltypesA = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
celltypesB = ["pre-Switched_Bm_JAM3","Plasma_IGHG1","gdT2_GZMH","CD8_CTL_GZMK","Transitional_B_MME","CD4_Tn_LIMS1","CD4_Tfh-like_CXCR5","gdT2_IL12RB2","NKT_IFNG","CD4_Th_TNFRSF11A","CD4_Treg_FCRL3","CD4_Tn_CCR7","Unswitched_Bm_CD1C","ncMono_FCGR3A","CD4_CTL_GZMH","CD4_Treg_FOXP3","gdT1_TRDV1","Terminal_NK_dim_CD160neg","CD4_Th17-like_RORC","DC1_CLEC9A","Unswitched_Bm_IL6","CD4_Th1-like_GZMK","cMono_CD14","gdT2_GZMK","Switched_Bm_IGHDneg","DC2_CD1C","CD8_Tcm_IFI44L","CD4_Th22-like_CCR10","CD4_Th_CCR4","MK_GP9","ncMono_C1QA","CD4_Th_CR1","Bn_IFIT3","DC_CSF2RA","CD8_Tn_SOX4","Switched_Bm_IGHE","Switched_activated_Bm_CD86","CD8_CTL_GZMB","Plasma_IGHA1","CD8_CTL_IFI44L","Bn_IL6","MAIT_SLC4A10","CD4_Tcm_IFI44L","Mature_NK_dim_FCGR3A","CD8_Tem_CCR7neg","HSPC_CD34","pDC_IRF4","NKT_NCR1","Plasmablast_MKI67","CD4_Tn_CXCR5","cMono_IL1B","CD8_Tn_CCR7","Transitional_NK_GZMK","CD4_Tn_SOX4","Proliferative_T_MKI67","cMono_CXCL10","intMono_GFRA2","Bn_TCL1A","cMono_IFI44L","ILC2_IL2RA","ncMono_IFIT1","CD4_Tr1_IL10","Atypical_Bm_ITGAX","Immature_T_CABP4","Proliferative_NK_MKI67","NK_bright_XCL1","CD4_Tem_CCR7neg","CD4_Th_LMNA","CD4_Tcm_CXCR5"]


logdir=curdir+"/logs/cima/"
jobdir=curdir+"/jobs/cima/"
outdir=curdir+"/output/cima/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

for ctA in celltypesA:
    for ctB in celltypesB:
        jobout = jobdir+ctA+"-"+ctB+".sh"
        print(jobout)
        fho = open(jobout,'w')
        for line in tpl:
            ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
            ol = ol.replace("SCRIPT","compare_CIMA_with_scMetaBrain.py")
            ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
            ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
            ol = ol.replace("CMD",f"{ctA} {ctB}")
            ol = ol.replace("OUTDIR",outdir)
            fho.write(ol)
        fho.close()
    
    # reverse lookup - disabled for now
    # fho = open(jobdir+"CIMA-"+ctA+".sh",'w')
    # for line in tpl:
    #     ol = line.replace("JOBNAME", "CIMA-vs-"+ctA)
    #     ol = ol.replace("SCRIPT","compare_CIMA_with_scMetaBrain.py")
    #     ol = ol.replace("OUTLOGFILE",logdir+"CIMA-vs-"+ctA+".log")
    #     ol = ol.replace("ERRLOGFILE",logdir+"CIMA-vs-"+ctA+".err")
    #     ol = ol.replace("CMD",f"CIMA {ctA}")
    #     fho.write(ol)
    # fho.close()




######## GTEx
celltypesA = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
celltypesB = ["Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Bladder","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_Cultured_fibroblasts","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood"]
logdir=curdir+"/logs/gtex/"
jobdir=curdir+"/jobs/gtex/"
outdir=curdir+"/output/gtex/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

for ctA in celltypesA:
    for ctB in celltypesB:
        jobout = jobdir+ctA+"-"+ctB+".sh"
        print(jobout)
        fho = open(jobout,'w')
        for line in tpl:
            ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
            ol = ol.replace("SCRIPT","compare_GTEx_with_scMetaBrain.py")
            ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
            ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
            ol = ol.replace("CMD",f"{ctA} {ctB}")
            ol = ol.replace("OUTDIR",outdir)
            fho.write(ol)
        fho.close()

######## SingleBrain

fh = open("./tpl-sb.tpl",'r')
tpl = fh.readlines()
fh.close()

# random FX data

celltypesA = ["SingleBrain_OPC","SingleBrain_Ext","SingleBrain_IN","SingleBrain_Ast","SingleBrain_OD","SingleBrain_MG","SingleBrain_End"]
matches = {
    "SingleBrain_OPC": "OPC",
    "SingleBrain_Ext": "EX",
    "SingleBrain_IN": "IN",
    "SingleBrain_Ast": "AST",
    "SingleBrain_OD": "OLI",
    "SingleBrain_MG": "MIC",
    "SingleBrain_End": "END"
}




for ctA in celltypesA:
    genelimitfile = f"/path/to/annotation/single_cell_eQTLsumstats/SingleBrain2026/data/{ctA}_eqtl_top_assoc.tsv.gz-genesWithSignificanceInRandomFxOnly.txt.gz"
    ctB = matches[ctA]
    logdir=curdir+"/logs/singlebrain-random/"
    jobdir=curdir+"/jobs/singlebrain-random/"
    outdir=curdir+"/output/singlebrain-random/"
    os.makedirs(jobdir,exist_ok=True)
    os.makedirs(logdir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)
    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","NO")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    fho = open(jobout,'w')
    print(jobout)        
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","NO")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()

    # random fx meta, limited to genes significant only in random fx meta
    logdir=curdir+"/logs/singlebrain-random-limited/"
    jobdir=curdir+"/jobs/singlebrain-random-limited/"
    outdir=curdir+"/output/singlebrain-random-limited/"
    os.makedirs(jobdir,exist_ok=True)
    os.makedirs(logdir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)

    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","NO")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT",genelimitfile)
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    fho = open(jobout,'w')
    print(jobout)        
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","NO")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT",genelimitfile)
        fho.write(ol)
    fho.close()


    # random FX data, flipped Z
    logdir=curdir+"/logs/singlebrain-random-flippedz/"
    jobdir=curdir+"/jobs/singlebrain-random-flippedz/"
    outdir=curdir+"/output/singlebrain-random-flippedz/"
    os.makedirs(jobdir,exist_ok=True)
    os.makedirs(logdir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)

    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","NO")
        ol = ol.replace("USEFLIPPEDZ","YES")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    fho = open(jobout,'w')
    print(jobout)        
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","NO")
        ol = ol.replace("USEFLIPPEDZ","YES")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()

    # random fx meta, limited to genes significant only in random fx meta
    logdir=curdir+"/logs/singlebrain-random-limited-flippedz/"
    jobdir=curdir+"/jobs/singlebrain-random-limited-flippedz/"
    outdir=curdir+"/output/singlebrain-random-limited-flippedz/"
    os.makedirs(jobdir,exist_ok=True)
    os.makedirs(logdir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)

    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","NO")
        ol = ol.replace("USEFLIPPEDZ","YES")
        ol = ol.replace("GENELIMIT",genelimitfile)
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    fho = open(jobout,'w')
    print(jobout)        
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","NO")
        ol = ol.replace("USEFLIPPEDZ","YES")
        ol = ol.replace("GENELIMIT",genelimitfile)
        fho.write(ol)
    fho.close()


    # fixed fx meta-analysis
    logdir=curdir+"/logs/singlebrain-fixed/"
    jobdir=curdir+"/jobs/singlebrain-fixed/"
    outdir=curdir+"/output/singlebrain-fixed/"

    os.makedirs(jobdir,exist_ok=True)
    os.makedirs(logdir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)

    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","YES")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    fho = open(jobout,'w')
    print(jobout)        
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","YES")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()


#############################
# Fujita replication
fh = open("./tpl-fujitabryois.tpl",'r')
tpl = fh.readlines()
fh.close()

celltypesA = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
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

logdir=curdir+"/logs/fujita-vs-woFujita/"
jobdir=curdir+"/jobs/fujita-vs-woFujita/"
outdir=curdir+"/output/fujita-vs-woFujita/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

for ctA in celltypesA:
    ctB = matches[ctA]

    jobname = "scMbWoFujita"+ctA+"-Fujita"+ctB
    jobout = jobdir+jobname+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", jobname)
        ol = ol.replace("SCRIPT","compare_Fujita_with_scMetaBrainWoFujita.py")
        ol = ol.replace("OUTLOGFILE",logdir+jobname+".log")
        ol = ol.replace("ERRLOGFILE",logdir+jobname+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","TRUE")
        fho.write(ol)
    fho.close()

    jobname = "Fujita"+ctA+"-scMbWoFujita"+ctB
    jobout = jobdir+jobname+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", jobname)
        ol = ol.replace("SCRIPT","compare_Fujita_with_scMetaBrainWoFujita.py")
        ol = ol.replace("OUTLOGFILE",logdir+jobname+".log")
        ol = ol.replace("ERRLOGFILE",logdir+jobname+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","FALSE")
        fho.write(ol)
    fho.close()

logdir=curdir+"/logs/fujita-vs-scMbFujita/"
jobdir=curdir+"/jobs/fujita-vs-scMbFujita/"
outdir=curdir+"/output/fujita-vs-scMbFujita/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

for ctA in celltypesA:
    ctB = matches[ctA]

    jobname = "scMbFujita"+ctA+"-Fujita"+ctB
    jobout = jobdir+jobname+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", jobname)
        ol = ol.replace("SCRIPT","compare_Fujita_with_scMetaBrainFujita.py")
        ol = ol.replace("OUTLOGFILE",logdir+jobname+".log")
        ol = ol.replace("ERRLOGFILE",logdir+jobname+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","TRUE")
        fho.write(ol)
    fho.close()

    jobname = "Fujita"+ctA+"-scMbFujita"+ctB
    jobout = jobdir+jobname+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", jobname)
        ol = ol.replace("SCRIPT","compare_Fujita_with_scMetaBrainFujita.py")
        ol = ol.replace("OUTLOGFILE",logdir+jobname+".log")
        ol = ol.replace("ERRLOGFILE",logdir+jobname+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","FALSE")
        fho.write(ol)
    fho.close()

# Bryois replication
fh = open("./tpl-fujitabryois.tpl",'r')
tpl = fh.readlines()
fh.close()

celltypesA = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
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

logdir=curdir+"/logs/bryois-vs-scMbWoBryois/"
jobdir=curdir+"/jobs/bryois-vs-scMbWoBryois/"
outdir=curdir+"/output/bryois-vs-scMbWoBryois/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

for ctA in celltypesA:
    ctB = matches[ctA]

    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_Bryois_with_scMetaBrainWoBryois.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","TRUE")
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_Bryois_with_scMetaBrainWoBryois.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","FALSE")
        fho.write(ol)
    fho.close()

logdir=curdir+"/logs/bryois-vs-scMbBryois/"
jobdir=curdir+"/jobs/bryois-vs-scMbBryois/"
outdir=curdir+"/output/bryois-vs-scMbBryois/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

for ctA in celltypesA:
    ctB = matches[ctA]

    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_Bryois_with_scMetaBrainBryois.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","TRUE")
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_Bryois_with_scMetaBrainBryois.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","FALSE")
        fho.write(ol)
    fho.close()


# diversecohorts
fh = open("./tpl-fujitabryois.tpl",'r')
logdir=curdir+"/logs/diversecohorts/"
jobdir=curdir+"/jobs/diversecohorts/"
outdir=curdir+"/output/diversecohorts/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

celltypesA = ["EX", "IN", "AST", "PER", "END", "OLI", "OPC", "MIC"]
for ctA in celltypesA:
    jobname = "DiverseCohorts-"+ctA+"-woDiverseCohorts"+ctA
    jobout = jobdir+jobname+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", jobname)
        ol = ol.replace("SCRIPT","compare_DiverseCohorts_with_scMetaBrainWoDiverseCohorts.py")
        ol = ol.replace("OUTLOGFILE",logdir+jobname+".log")
        ol = ol.replace("ERRLOGFILE",logdir+jobname+".err")
        ol = ol.replace("CMD",f"{ctA} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","TRUE")
        fho.write(ol)
    fho.close()

    jobname = "woDiverseCohorts-"+ctA+"-DiverseCohorts"+ctA
    jobout = jobdir+jobname+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", jobname)
        ol = ol.replace("SCRIPT","compare_DiverseCohorts_with_scMetaBrainWoDiverseCohorts.py")
        ol = ol.replace("OUTLOGFILE",logdir+jobname+".log")
        ol = ol.replace("ERRLOGFILE",logdir+jobname+".err")
        ol = ol.replace("CMD",f"{ctA} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("CTAISSCMB","FALSE")
        fho.write(ol)
    fho.close()

# diversecohorts vs SingleBrain
logdir=curdir+"/logs/diversecohorts/"
jobdir=curdir+"/jobs/diversecohorts/"
outdir=curdir+"/output/diversecohorts/"
os.makedirs(jobdir,exist_ok=True)
os.makedirs(logdir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)


######## SingleBrain
celltypesA = ["SingleBrain_OPC","SingleBrain_Ext","SingleBrain_IN","SingleBrain_Ast","SingleBrain_OD","SingleBrain_MG","SingleBrain_End"]
matches = {
    "SingleBrain_OPC": "OPC",
    "SingleBrain_Ext": "EX",
    "SingleBrain_IN": "IN",
    "SingleBrain_Ast": "AST",
    "SingleBrain_OD": "OLI",
    "SingleBrain_MG": "MIC",
    "SingleBrain_End": "END"
}

fh = open("./tpl-sb.tpl",'r')
tpl = fh.readlines()
fh.close()

# compare SingleBrain with DiverseCohorts
for ctA in celltypesA:
    genelimitfile = f"/path/to/annotation/single_cell_eQTLsumstats/SingleBrain2026/data/{ctA}_eqtl_top_assoc.tsv.gz-genesWithSignificanceInRandomFxOnly.txt.gz"
    logdir=curdir+"/logs/singlebrain-fixed-diversecohorts/"
    jobdir=curdir+"/jobs/singlebrain-fixed-diversecohorts/"
    outdir=curdir+"/output/singlebrain-fixed-diversecohorts/"
    os.makedirs(jobdir,exist_ok=True)
    os.makedirs(logdir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)
    
    ctB = matches[ctA]
    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_DiverseCohorts.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","YES")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_DiverseCohorts.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","YES")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()


# single brain vs scMetaBrain-SingleBrain
celltypesA = ["SingleBrain_OPC","SingleBrain_Ext","SingleBrain_IN","SingleBrain_Ast","SingleBrain_OD","SingleBrain_MG","SingleBrain_End"]
matches = {
    "SingleBrain_OPC": "OPC",
    "SingleBrain_Ext": "EX",
    "SingleBrain_IN": "IN",
    "SingleBrain_Ast": "AST",
    "SingleBrain_OD": "OLI",
    "SingleBrain_MG": "MIC",
    "SingleBrain_End": "END"
}

fh = open("./tpl-sb.tpl",'r')
tpl = fh.readlines()
fh.close()

# compare SingleBrain with DiverseCohorts
for ctA in celltypesA:
    genelimitfile = f"/path/to/single_cell_eQTLsumstats/SingleBrain2026/data/{ctA}_eqtl_top_assoc.tsv.gz-genesWithSignificanceInRandomFxOnly.txt.gz"
    logdir=curdir+"/logs/singlebrain-fixed-scmetabrainsinglebrain/"
    jobdir=curdir+"/jobs/singlebrain-fixed-scmetabrainsinglebrain/"
    outdir=curdir+"/output/singlebrain-fixed-scmetabrainsinglebrain/"
    os.makedirs(jobdir,exist_ok=True)
    os.makedirs(logdir,exist_ok=True)
    os.makedirs(outdir,exist_ok=True)
    
    ctB = matches[ctA]
    jobout = jobdir+ctA+"-"+ctB+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrainSingleBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
        ol = ol.replace("CMD",f"{ctA} {ctB}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","YES")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()

    jobout = jobdir+ctB+"-"+ctA+".sh"
    print(jobout)
    fho = open(jobout,'w')
    for line in tpl:
        ol = line.replace("JOBNAME", ctB+"-vs-"+ctA)
        ol = ol.replace("SCRIPT","compare_SingleBrain_with_scMetaBrainSingleBrain.py")
        ol = ol.replace("OUTLOGFILE",logdir+ctB+"-vs-"+ctA+".log")
        ol = ol.replace("ERRLOGFILE",logdir+ctB+"-vs-"+ctA+".err")
        ol = ol.replace("CMD",f"{ctB} {ctA}")
        ol = ol.replace("OUTDIR",outdir)
        ol = ol.replace("USEFIXEDFDR","YES")
        ol = ol.replace("USEFLIPPEDZ","NO")
        ol = ol.replace("GENELIMIT","\n")
        fho.write(ol)
    fho.close()

# # random FX data
# for ctA in celltypesA:
#     genelimitfile = f"/path/to/single_cell_eQTLsumstats/SingleBrain2026/data/{ctA}_eqtl_top_assoc.tsv.gz-genesWithSignificanceInRandomFxOnly.txt.gz"
#     logdir=curdir+"/logs/singlebrain-random-diversecohorts/"
#     jobdir=curdir+"/jobs/singlebrain-random-diversecohorts/"
#     outdir=curdir+"/output/singlebrain-random-diversecohorts/"
#     os.makedirs(jobdir,exist_ok=True)
#     os.makedirs(logdir,exist_ok=True)
#     os.makedirs(outdir,exist_ok=True)
    
#     ctB = matches[ctA]
#     jobout = jobdir+ctA+"-"+ctB+".sh"
#     print(jobout)
#     fho = open(jobout,'w')
#     for line in tpl:
#         ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
#         ol = ol.replace("SCRIPT","compare_SingleBrain_with_DiverseCohorts.py")
#         ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
#         ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
#         ol = ol.replace("CMD",f"{ctA} {ctB}")
#         ol = ol.replace("OUTDIR",outdir)
#         ol = ol.replace("USEFIXEDFDR","NO")
#         ol = ol.replace("USEFLIPPEDZ","NO")
#         ol = ol.replace("GENELIMIT","\n")
#         fho.write(ol)
#     fho.close()

#     # random fx meta, limited to genes significant only in random fx meta
#     logdir=curdir+"/logs/singlebrain-random-limited-diversecohorts/"
#     jobdir=curdir+"/jobs/singlebrain-random-limited-diversecohorts/"
#     outdir=curdir+"/output/singlebrain-random-limited-diversecohorts/"
#     os.makedirs(jobdir,exist_ok=True)
#     os.makedirs(logdir,exist_ok=True)
#     os.makedirs(outdir,exist_ok=True)

#     jobout = jobdir+ctA+"-"+ctB+".sh"
#     print(jobout)
#     fho = open(jobout,'w')

#     for line in tpl:
#         ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
#         ol = ol.replace("SCRIPT","compare_SingleBrain_with_DiverseCohorts.py")
#         ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
#         ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
#         ol = ol.replace("CMD",f"{ctA} {ctB}")
#         ol = ol.replace("OUTDIR",outdir)
#         ol = ol.replace("USEFIXEDFDR","NO")
#         ol = ol.replace("USEFLIPPEDZ","NO")
#         ol = ol.replace("GENELIMIT",genelimitfile)
#         fho.write(ol)
#     fho.close()

#     # random FX data, flipped Z
#     logdir=curdir+"/logs/singlebrain-random-flippedz-diversecohorts/"
#     jobdir=curdir+"/jobs/singlebrain-random-flippedz-diversecohorts/"
#     outdir=curdir+"/output/singlebrain-random-flippedz-diversecohorts/"
#     os.makedirs(jobdir,exist_ok=True)
#     os.makedirs(logdir,exist_ok=True)
#     os.makedirs(outdir,exist_ok=True)


#     jobout = jobdir+ctA+"-"+ctB+".sh"
#     print(jobout)
#     fho = open(jobout,'w')
#     for line in tpl:
#         ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
#         ol = ol.replace("SCRIPT","compare_SingleBrain_with_DiverseCohorts.py")
#         ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
#         ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
#         ol = ol.replace("CMD",f"{ctA} {ctB}")
#         ol = ol.replace("OUTDIR",outdir)
#         ol = ol.replace("USEFIXEDFDR","NO")
#         ol = ol.replace("USEFLIPPEDZ","YES")
#         ol = ol.replace("GENELIMIT","\n")
#         fho.write(ol)
#     fho.close()

#     # random fx meta, limited to genes significant only in random fx meta
#     logdir=curdir+"/logs/singlebrain-random-limited-flippedz-diversecohorts/"
#     jobdir=curdir+"/jobs/singlebrain-random-limited-flippedz-diversecohorts/"
#     outdir=curdir+"/output/singlebrain-random-limited-flippedz-diversecohorts/"
#     os.makedirs(jobdir,exist_ok=True)
#     os.makedirs(logdir,exist_ok=True)
#     os.makedirs(outdir,exist_ok=True)


#     jobout = jobdir+ctA+"-"+ctB+".sh"
#     print(jobout)
#     fho = open(jobout,'w')
#     for line in tpl:
#         ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
#         ol = ol.replace("SCRIPT","compare_SingleBrain_with_DiverseCohorts.py")
#         ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
#         ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
#         ol = ol.replace("CMD",f"{ctA} {ctB}")
#         ol = ol.replace("OUTDIR",outdir)
#         ol = ol.replace("USEFIXEDFDR","NO")
#         ol = ol.replace("USEFLIPPEDZ","YES")
#         ol = ol.replace("GENELIMIT",genelimitfile)
#         fho.write(ol)
#     fho.close()

#     # fixed fx meta-analysis
#     logdir=curdir+"/logs/singlebrain-fixed-diversecohorts/"
#     jobdir=curdir+"/jobs/singlebrain-fixed-diversecohorts/"
#     outdir=curdir+"/output/singlebrain-fixed-diversecohorts/"

#     os.makedirs(jobdir,exist_ok=True)
#     os.makedirs(logdir,exist_ok=True)
#     os.makedirs(outdir,exist_ok=True)

#     jobout = jobdir+ctA+"-"+ctB+".sh"
#     print(jobout)
#     fho = open(jobout,'w')
#     for line in tpl:
#         ol = line.replace("JOBNAME", ctA+"-vs-"+ctB)
#         ol = ol.replace("SCRIPT","compare_SingleBrain_with_DiverseCohorts.py")
#         ol = ol.replace("OUTLOGFILE",logdir+ctA+"-vs-"+ctB+".log")
#         ol = ol.replace("ERRLOGFILE",logdir+ctA+"-vs-"+ctB+".err")
#         ol = ol.replace("CMD",f"{ctA} {ctB}")
#         ol = ol.replace("OUTDIR",outdir)
#         ol = ol.replace("USEFIXEDFDR","YES")
#         ol = ol.replace("USEFLIPPEDZ","NO")
#         ol = ol.replace("GENELIMIT","\n")
#         fho.write(ol)
#     fho.close()