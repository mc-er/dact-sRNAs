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

enriched_GO <- function(id2go, fp_testGenes, testGroup){
  geneIDs <- names(id2go)
  genes2test <- read.table(fp_testGenes, sep="\t", header=F)
  fac_table <- as.factor(geneIDs) %in% genes2test$V1 %>% as.integer() %>% factor()
  names(fac_table) <- geneIDs
  
  GOdata <- new("topGOdata", ontology = testGroup, allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  GOtest <- runTest(GOdata, statistic = "fisher")
  
  allRes <- GenTable(GOdata, weight01_pval=GOtest, orderBy = "weight01", ranksOf = "weight01", topNodes = 100, numChar = 1000) %>%
    mutate(category = testGroup) %>%
    dplyr::rename(ID = "GO.ID")
  
  allGO <- genesInTerm(GOdata)
  
  return(allRes)
}

list_genes_for_GO_ID <- function(id2go, go_out){
  GO2geneID <- inverseList(id2go)
  d <- lapply(GO2geneID[go_out$ID][go_out[,1]], paste0, collapse = ", ") %>%
    as.data.frame() %>%
    pivot_longer(cols = starts_with("G"), names_to = "ID", values_to = "genes") %>%
    mutate(ID = gsub("\\.", ":", ID))
  return(d)
} 

add_geneID_to_GOout <- function(go_out, GO_gene_df){
  d <- merge(go_out, GO_gene_df) %>%
    dplyr::rename(adj_pval = "weight01_pval") %>%
    dplyr::select(category, ID , Term, adj_pval, genes) %>%
    mutate(adj_pval = as.numeric(gsub("< ", "", adj_pval)),
           genes = strsplit(genes, ", ")) %>%
    unnest(genes) %>%
    distinct()
  return(d)
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