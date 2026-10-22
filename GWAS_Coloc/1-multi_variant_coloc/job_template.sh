gene=<GENE>
gwas=<GWAS>
cell_type=<CELL_TYPE>
eqtl_path=<EQTL_PATH>
sig_eqtl_path=<SIG_EQTL_PATH>
gwas_path=<GWAS_PATH>
gwas_summary_path=<GWAS_CS_SUMMARY_PATH>
out_dir=<OUT_DIR>
annotation_path=<ANNOTATION_PATH>

# Merge GWAS and eQTLs
python3 ./scripts/merge.py \
    ${gene} \
    ${eqtl_path} \
    ${gwas_path} \
    ${cell_type} \
    ${gwas}

# Coloc LBF()
Rscript ./scripts/coloc_lbf.R \
    ${gene} \
    ${gwas}_${cell_type}_${gene}_merged.txt \
    ${cell_type} \
    ${gwas}

# Format output
python3 ./scripts/format_mv_output.py \
    ${gene} \
    ${cell_type} \
    ${gwas_path} \
    ${eqtl_path} \
    ${gwas}_${cell_type}_${gene}_mv_coloc_out.txt \
    ${gwas}_${cell_type}_${gene}_merged.txt \
    ${sig_eqtl_path} \
    ${gwas} \
    ${gwas_summary_path} \
    ${annotation_path}


mv ${gwas}_${cell_type}_${gene}_mv_coloc_out_formatted.txt ${out_dir}
rm ${gwas}_${cell_type}_${gene}_mv_coloc_out.txt
rm ${gwas}_${cell_type}_${gene}_merged.txt