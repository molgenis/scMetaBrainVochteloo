library(ashr)
library(mashr)
library(corrplot)
library(glue)

m2 <- readRDS("../2-run_mashr/m2.rds")

discovery_datasets <- c(
    'AST',
    'END',
    'EX',
    'IN',
    'MIC',
    'OLI',
    'OPC',
    'PER'
)

top_effects <- read.csv('../1-data/strong/top_effects_all_cell_types.txt', sep='\t')

for (discovery_dataset in discovery_datasets) {
    m2_copy <- m2

    # Get a list top effects for current cell type
    top_effects_filtered <- top_effects[top_effects$cell_type == discovery_dataset,]
    discovery_snps <- top_effects_filtered[['snp_gene']]
    print(discovery_dataset)

    # Filter mashr object for top effects for current cell type
    m2_copy$result$PosteriorMean = m2$result$PosteriorMean[rownames(m2_copy$result$PosteriorMean) %in% discovery_snps, ]
    m2_copy$result$PosteriorSD = m2$result$PosteriorSD[rownames(m2_copy$result$PosteriorSD) %in% discovery_snps, ]
    m2_copy$result$NegativeProb = m2$result$NegativeProb[rownames(m2_copy$result$NegativeProb) %in% discovery_snps, ]
    m2_copy$result$lfsr = m2$result$lfsr[rownames(m2_copy$result$lfsr) %in% discovery_snps, ]

    print(str(m2_copy$result$PosteriorMean))

    # Get pairwise sharing values for cell type
    x = get_pairwise_sharing(m2_copy, factor=0, lfsr_thresh=0.05)

    # Write to file
    write.table(x,
                file = glue("output/{discovery_dataset}.txt"),
                sep = "\t",
                quote = FALSE,
                col.names = NA)

}


