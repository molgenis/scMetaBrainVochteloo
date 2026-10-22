gene=<GENE>
chromosome=<CHROMOSOME>
gwas=<GWAS>
cell_type=<CELL_TYPE>
out_dir=<OUT_DIR>

gwas_path=<GWAS_PATH>
gwas_p_col=<GWAS_P_COL>
gwas_beta_col=<GWAS_BETA_COL>
gwas_se_col=<GWAS_SE_COL>
gwas_pos_col=<GWAS_POS_COL>
gwas_chrom_col=<GWAS_CHROM_COL>
gwas_ea_col=<GWAS_EA_COL>
gwas_oa_col=<GWAS_OA_COL>
gwas_snp_col=<GWAS_SNP_COL>
gwas_eaf_col=<GWAS_EAF_COL>
gwas_type=<GWAS_TYPE>
gwas_n=<GWAS_N>

eqtl_path=<EQTL_PATH>
sig_eqtl_path=<SIG_EQTL_PATH>
eqtl_n=<EQTL_N>

# Get GWAS data
python3 ./scripts/get_gwas_data.py ${gene} \
    ${chromosome} \
    ${gwas} \
    ${gwas_path} \
    ${gwas_p_col} \
    ${gwas_beta_col} \
    ${gwas_se_col} \
    ${gwas_pos_col} \
    ${gwas_chrom_col} \
    ${gwas_ea_col} \
    ${gwas_oa_col} \
    ${gwas_snp_col} \
    ${cell_type} \
    ${annotation_path}


# Merge with eQTLs
python3 ./scripts/merge.py ${gene} \
    ${chromosome} \
    ${cell_type} \
    ${eqtl_path} \
    ${gwas} \
    ${gwas}_${cell_type}_${gene}_gwas.txt \
    ${gwas_eaf_col} \
    ${annotation_path}


# Coloc
Rscript ./scripts/coloc_abf.R ${gene} \
    ${gwas}_${cell_type}_${gene}_merged.txt \
    ${gwas_n} \
    ${gwas_type} \
    ${cell_type} \
    ${gwas} \
    ${eqtl_n}


# Format output
python3 ./scripts/format_sv_output.py ${gene} \
    ${cell_type} \
    ${gwas}_${cell_type}_${gene}_gwas.txt \
    ${gwas}_${cell_type}_${gene}_reduced.txt \
    ${gwas}_${cell_type}_${gene}_sv_coloc_out.txt \
    ${gwas}_${cell_type}_${gene}_merged.txt \
    ${sig_eqtl_path} \
    ${gwas} \
    ${out_dir} \
    ${annotation_path}

rm ${gwas}_${cell_type}_${gene}_sv_coloc_out.txt
rm ${gwas}_${cell_type}_${gene}_merged.txt
rm ${gwas}_${cell_type}_${gene}_gwas.txt
rm ${gwas}_${cell_type}_${gene}_reduced.txt

