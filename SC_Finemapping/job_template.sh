gene=<GENE>
chromosome=<CHROMOSOME>
cell_type=<CELL_TYPE>
sample_size=<SAMPLE_SIZE>
eqtl_path=<EQTL_PATH>
sig_eqtl_path=<SIG_EQTL_PATH>
vcf_path=<VCF_PATH>
gte_path=<GTE_PATH>
out_dir=<OUT_DIR>


mkdir ${cell_type}_${gene}
cd ${cell_type}_${gene}

# Get eQTLs
python3 ./scripts/get_data.py ${gene} \
    ${chromosome} \
    ${eqtl_path} \
    ${sig_eqtl_path} \
    ${sample_size}

# Calculate LD matrix
python3 ./scripts/calculate_ld_matrix.py ${gene} \
    ${chromosome} \
    ${gene}_eqtls.txt \
    ${cell_type} \
    ${vcf_path} \
    ${gte_path}

# Run SuSiE
Rscript ./scripts/susie.R ${gene} \
    ${gene}_eqtls.txt \
    ${gene}_ld.txt \
    ${sample_size} \
    5

# Filter
python3 ./scripts/filter.py ${gene}_susie.txt \
        ${gene} \
        ${cell_type} \
        ${vcf_path} \
        ${sig_eqtl_path}

# Get number of credible sets
num_cs_found=$(zcat ${gene}_susie_filtered.txt.gz | awk -F'\t' 'NR>1 && $14 != "" && tolower($14)!="nan" {print $14}' | sort -u | wc -l)

# If no credible sets are found, repeat finemapping with 1 credible set
if [ "${num_cs_found}" -eq 0 ]; then
    rm ${gene}_susie_filtered.txt.gz
    rm ${gene}_susie.txt

    Rscript ./scripts/susie.R ${gene} \
        ${gene}_eqtls.txt \
        ${gene}_ld.txt \
        ${sample_size} \
        1
    
    python3 ./scripts/filter.py ${gene}_susie.txt \
        ${gene} \
        ${cell_type} \
        ${vcf_path} \
        ${sig_eqtl_path}
fi

# If 5 credible sets are found, repeat finemapping with 10 credible sets
if [ "${num_cs_found}" -eq 5 ]; then
    rm ${gene}_susie_filtered.txt.gz
    rm ${gene}_susie.txt

    Rscript ./scripts/susie.R ${gene} \
        ${gene}_eqtls.txt \
        ${gene}_ld.txt \
        ${sample_size} \
        10

    python3 ./scripts/filter.py ${gene}_susie.txt \
        ${gene} \
        ${cell_type} \
        ${vcf_path} \
        ${sig_eqtl_path}
fi 


mv ${gene}_susie_filtered.txt.gz ${out_dir}/${cell_type}
rm *

cd ..
rm -r ${cell_type}_${gene}