###################################################
# differential targeting by sRNAs edgeR functions #
###################################################
filter_smRNA_by_exp <- function(D, min_samples, cpm_level, sample_prop){
  # save a copy of the original data
  d <- D
  # transform the counts to CPM
  D[,c(2:ncol(D))] <- cpm(D[,c(2:ncol(D))])
  # filter the data by the number of groups and the proportion of samples that have a CPM value above the cpm_level
  # - pivot the data to long format
  # - make a species column based on the first letter in the sample name
  # - group by peakID and species
  # - make three new columns:
  #   - cpm_filter_count: the number of samples with a CPM value above the cpm_level
  #   - cpm_filter_prop: the proportion of samples with a CPM value above the cpm_level
  #   - filterPASS: TRUE if the cpm_filter_count is greater than or equal to min_samples and the cpm_filter_prop is greater than or equal to sample_prop
  # - ungroup the data
  # - select the peakID, species, and filterPASS columns
  # - remove duplicates
  # - filter the data to only include peakIDs that have a filterPASS value larger than min_samples
  # - pivot the data to wide format
  # - return the original dataframe but only containing the peaks passing the filters
  k <- D  %>%  
    pivot_longer(!peakID, names_to = "samples", values_to = "CPM") %>% 
    mutate(species = substr(samples,1,1)) %>% 
    group_by(peakID, species) %>% 
    mutate(cpm_filter_count = sum(CPM > cpm_level),
           cpm_filter_prop = sum(CPM > cpm_level)/length(species),
           filterPASS = (cpm_filter_count >= 2 & cpm_filter_prop >= sample_prop)) %>%
    ungroup(species) %>% 
    dplyr::select(peakID, species, filterPASS) %>% 
    distinct() %>%
    filter(sum(filterPASS == FALSE) < min_samples) %>%
    pivot_wider(names_from = species, values_from = filterPASS)
  
  d <- d[d$peakID%in%k$peakID,]
  return(d)
}

####################
# TOP GO functions #
####################

enriched_GO <- function(go_db_fp, genes_list_fp, go_test_category, algorithm, statistic){
  # load GO database
  geneID2GO <- readMappings(go_db_fp)
  
  # extract gene IDs
  geneIDs <- names(geneID2GO)
  
  # load lit of genes to test
  genes2test <- read.table(genes_list_fp, sep="\t", header=F)
  # make a table with all genes marking genes to test with a one (other genes with a zero)
  fac_table <- as.factor(geneIDs) %in% genes2test$V1 %>% as.integer() %>% factor()
  names(fac_table) <- geneIDs
  
  # set up the GO enrichment analysis
  GOdata <- new("topGOdata", ontology = go_test_category, allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  # preform the GO enrichment
  GOtest <- runTest(GOdata, algorithm = algorithm, statistic = statistic)
  
  # extract result and adjust p-values
  #  - add column with the category (BP, MF or CC)
  allRes <- GenTable(GOdata, weight01_pval=GOtest, orderBy = "weight01", ranksOf = "weight01", numChar = 1000) %>%
    mutate(GO_category = go_test_category)
  
  # make a new column that contains the genes in enriched GO term
  # extract GO ids in the result data
  allRes_GOs <- allRes$ID
  
  # extract genes for each GO term
  allGO <- genesInTerm(GOdata)
  
  # - extract all go and genes associated with them
  # - make into a data frame
  # - pivot GO IDs into a column instead of one column per GO ID
  # - replace . in GO id with :
  # - split gene column by , and place each gene on it's own line
  # - remove all genes not in the genes to test list
  # - group by GO ID
  # - collapse genes into a , separated list
  # - remove duplicated rows
  goID_genes <- lapply(allGO[allRes$GO.ID][allRes[,1]], paste0, collapse = ", ") %>%
    as.data.frame() %>%
    pivot_longer(cols = starts_with("G"), names_to = "GO.ID", values_to = "genes") %>%
    mutate(GO.ID = gsub("\\.", ":", GO.ID)) %>%
    separate_longer_delim(genes, delim = ", ") %>%
    filter(genes %in% genes2test$V1) %>%
    group_by(GO.ID) %>%
    mutate(genes = paste0(genes, collapse = ";")) %>%
    distinct()

  # add test genes to result dataframe
  #  - rename columns
  #  - make adjusted p values into a umberical column
  #  - rearrange columns
  goOUT_genes <- merge(allRes, goID_genes) %>%
    dplyr::rename("GO_id" = "GO.ID",
                  "GO_term" = "Term",
                  "total_annotated_genes" = "Annotated",
                  "expected_genes" = "Expected",
                  "significant_genes" = "Significant",
                  "adj_pval" = "weight01_pval",
                  "id_significant_genes" = "genes")  %>%
    mutate(adj_pval = as.numeric(gsub("< ", "", adj_pval)))  %>%
    dplyr::select(GO_category, GO_id , GO_term, adj_pval, total_annotated_genes, expected_genes, significant_genes, id_significant_genes)

  return(goOUT_genes)
}

calculate_zscore_Ngenes_Npeaks <- function(fp_edgeR_DTpeaks, outGO_genesIDs, filter_pattern){
  DTpeaks <- read.table(fp_edgeR_DTpeaks, sep="\t", header=T) %>% 
    filter(FDR <= 0.05 & grepl(filter_pattern, region)) %>% na.omit() %>% dplyr::rename(genes = "geneID")

  d <- distinct(merge(outGO_genesIDs, DTpeaks, all=F)) %>%
    mutate(logFC = as.numeric(logFC)) %>%
    group_by(ID) %>%
    mutate(zscore = sum(logFC)/sqrt(length(ID)),
           n_peaks = length(unique(peakID)),
           n_genes = length(unique(genes))) %>%
    ungroup(ID) %>%
    dplyr::select("peakID", "genes", "cdsID",	"region", "category",	"ID",	"Term",	"adj_pval",	"zscore",	"n_peaks",	"n_genes") %>%
    distinct()
  return(d)
}