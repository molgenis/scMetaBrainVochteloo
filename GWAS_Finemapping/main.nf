nextflow.enable.dsl=2


process getSignificantGWASSnps {
    time '6h'
    memory '4 GB'
    cpus 1

    output:
    path "intervals.txt"

    script:
    """
    extract_sig_gwas_snps.py ${params.gwas_path} \
        ${params.gwas_chromosome_column} \
        ${params.gwas_position_column} \
        ${params.gwas_pvalue_column} \
        ${params.sig_cutoff} \
        ${params.window_size}
    """
}

process selectGenes {
    time '6h'
    memory '8 GB'
    cpus 1

    input:
    path intervalFile

    output:
    path "gene_list.txt"

    script:
    """
    get_gene_list.py ${params.annotation_path} \
        ${params.tested_genes_list_path}
    """
}

process runSusie {
    publishDir "${params.out_dir}/", mode: 'copy'

    time '24h'
    memory { 16.GB * 2 ** (task.attempt - 1) }
    cpus 1

    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    maxRetries 5

    input:
    tuple val(gene), val(chromosome)

    output:
    path "${gene}_susie_filtered.txt.gz"

    shell:
    '''
    # Get overlapping variants
    get_overlapping_variants.py !{gene} \
    !{chromosome} \
    !{params.tested_eqtl_variant_path}/chr!{chromosome}.txt.gz \
    !{params.ref_dir}/chr!{chromosome} \
    !{params.eqtl_vcf_path}

    # Get GWAS variants
    get_gwas_variants.py !{gene} \
        !{params.gwas_pvalue_column} \
        !{params.gwas_snp_column} \
        !{params.gwas_beta_column} \
        !{params.gwas_standard_error_column} \
        !{params.gwas_chromosome_column} \
        !{params.gwas_effect_allele_column} \
        !{params.gwas_other_allele_column} \
        !{params.gwas_position_column} \
        !{gene}_eqtl_variants_alleles_alpha_sorted.txt \
        !{params.gwas_path} \
        !{params.ref_dir}/chr!{chromosome}.snplist \
        !{params.gwas_allele_freq_column} \
        !{chromosome}

    # Calculate LD matrix
    calculate_ld_matrix.py !{gene} \
        !{gene}_gwas.txt \
        !{params.ref_dir}/chr!{chromosome}.vcf.gz

    # GWAS finemapping
    susie.R !{gene} \
        !{gene}_gwas.txt \
        !{gene}_ld.npy \
        !{params.gwas_n} \
        5

    # Filter results
    filter.py !{gene} \
        !{gene}_susie.txt \
        !{params.ref_dir}/chr!{chromosome}.vcf.gz \
        !{chromosome} \
        !{params.sig_cutoff}


    # Get the SusieRss_CS column number
    cs_col=$(zcat !{gene}_susie_filtered.txt.gz | head -n1 | awk -F'\t' '{for(i=1;i<=NF;i++) if(\$i=="SusieRss_CS"){print i; exit}}')
    
    # Get the number of unique credible sets
    #num_cs_found=$(zcat ${gene}_susie_filtered.txt.gz | awk -F'\t' 'NR>1 && $cs_col != "" && tolower($cs_col)!="nan" {print $cs_col}' | sort -u | wc -l)
    
    num_cs_found=$(zcat !{gene}_susie_filtered.txt.gz | \
    awk -F'\t' -v col="\$cs_col" \
    'NR>1 && $(col) != "" && tolower($(col))!="nan" {print $(col)}' | \
    sort -u | wc -l)
    echo $num_cs_found

    if [ "${num_cs_found}" -eq 0 ]; then
        rm !{gene}_susie_filtered.txt.gz
        rm !{gene}_susie.txt

        # GWAS finemapping
        susie.R !{gene} \
            !{gene}_gwas.txt \
            !{gene}_ld.npy \
            !{params.gwas_n} \
            1

        # Filter results
        filter.py !{gene} \
            !{gene}_susie.txt \
            !{params.ref_dir}/chr!{chromosome}.vcf.gz \
            !{chromosome} \
            !{params.sig_cutoff}
    fi

    if [ "${num_cs_found}" -eq 5 ]; then
        rm !{gene}_susie_filtered.txt.gz
        rm !{gene}_susie.txt

        # GWAS finemapping
        susie.R !{gene} \
            !{gene}_gwas.txt \
            !{gene}_ld.npy \
            !{params.gwas_n} \
            10

        # Filter results
        filter.py !{gene} \
            !{gene}_susie.txt \
            !{params.ref_dir}/chr!{chromosome}.vcf.gz \
            !{chromosome} \
            !{params.sig_cutoff}
    fi

    rm !{gene}_gwas.txt
    rm !{gene}_ld.npy
    rm !{gene}_susie.txt
    '''
}

workflow {
    // Get list file with significant GWAS SNPs
    getSignificantGWASSnps()

    // Select genes to be tested
    gene_list_channel = selectGenes(getSignificantGWASSnps.out)

    gene_channel = gene_list_channel.flatMap { file ->
        file.readLines()
            .findAll { it.trim() }
            .collect { line ->
                def (gene, chromosome) = line.split('\t')
                tuple(gene, chromosome)
            }
    }.filter { gene, chromosome ->
        !file("${params.out_dir}/${gene}_susie_filtered.txt.gz").exists()
    }

    runSusie(gene_channel)
}
