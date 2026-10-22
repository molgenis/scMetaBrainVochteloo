gene=<GENE>
ds1=<DS1>
path_ds1=<PATH_DS1>
ds1_sig_path=<SIG_PATH_DS1>
ds1_type=<DS1_TYPE>
ds2=<DS2>
path_ds2=<PATH_DS2>
ds2_sig_path=<SIG_PATH_DS2>
ds2_type=<DS2_TYPE>
annotation_path=<ANNOTATION_PATH>

python3 ./scripts/merge.py \
    ${ds1} \
    ${path_ds1} \
    ${ds2} \
    ${path_ds2} \
    ${gene} \
    ${annotation_path}

Rscript ./scripts/coloc_lbf.R \
    ${gene}_${ds1}_${ds2}.txt \
    ${gene} \
    ${ds1} \
    ${ds2}

python3 ./scripts/format_mv_output.py \\
    ${gene} \\
    ${ds1} \\
    ${ds2} \\
    ${path_ds1} \\
    ${path_ds2} \\
    ${gene}_${ds1}_${ds2}_mv_coloc_out.txt \\
    ${gene}_${ds1}_${ds2}.txt \\
    ${ds1_sig_path} \\
    ${ds2_sig_path} \\
    ${ds1_type} \\
    ${ds2_type}


rm ${gene}_${ds1}_${ds2}_sv_coloc_out.txt
rm ${gene}_${ds1}_${ds2}_mv_coloc_out.txt
rm ${gene}_${ds1}_${ds2}.txt