# packages
library(tidyverse)
library(ape)
library(parallel)

setwd(Sys.getenv("NG_MANUSCRIPT_ROOT", "~/Desktop/ng_mount/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/"))

# distance matrix
#distmat <- readRDS("./probe_creation/2025-01-14_regentrans/data/dist_mat_formatted.RDS") 
distmat <- readRDS("./data/phylogeny/source/dist_mat.RDS")

# tree
tree <- read.tree("./data/phylogeny/source/cognac_nj.tre")

# relevant metadata
ng_amrsearch <- read_csv("../../data/pathogenwatch_metadata/Neisseria_gonorrhoeae_amrsearch.csv") %>% subset(`Genome Name` %in% tree$tip.label)
ng_snps <- readRDS("../../data/ng_pathogenwatch_files/ng_snps.RDS") %>% subset(`Genome.Name` %in% tree$tip.label)
ng_meta <- left_join(ng_amrsearch[,c(-1,-3,-4)], ng_snps[,c(-1,-3,-4)], by = join_by(`Genome Name` == `Genome.Name`))
colnames(ng_meta)[1] <- "Isolate"
rownames(ng_meta) <- ng_meta$Isolate

# ng snp gene grouping

snp_data <- ng_snps[,5:66]
prefixes <- unique(str_extract(colnames(snp_data), "^[^_]+"))

collapsed_snps <- as.data.frame(sapply(prefixes, function(prefix) {
  matched_cols <- grep(paste0("^", prefix, "_"), colnames(snp_data), value = TRUE)
  if (length(matched_cols) > 0) {
    rowSums(snp_data[, matched_cols, drop = FALSE], na.rm = TRUE)
  } else {
    snp_data[, prefix, drop = FALSE]
  }
}))

rownames(collapsed_snps) <- rownames(ng_snps)
collapsed_snps_binary <- collapsed_snps %>% mutate_all(~ as.factor(ifelse(. > 0, "1", "0")))
colnames(collapsed_snps_binary)[colnames(collapsed_snps_binary) == "X16S"] <- "16S_rDNA"
colnames(collapsed_snps_binary)[colnames(collapsed_snps_binary) == "X23S"] <- "23S_rDNA"

ng_meta_gene_grouped <- left_join(ng_meta[, 1:9], collapsed_snps_binary %>% as.data.frame %>% mutate(Isolate = rownames(.))) %>% as.data.frame()
rownames(ng_meta_gene_grouped) <- ng_meta_gene_grouped$Isolate

# ---- functions ----

#get_tn_clusters_SNPthresh - Takes as input: 1) dist_mat - A pairwise genetic distance matrix
#                      2) thresh  - A genetic distance threshold
#
#Returns distance clusters where all pairs of isolates within a cluster are below the threshold
get_tn_clusters_SNPthresh <- function(dist_mat, thresh)
{
  dist_clust = hclust(as.dist(dist_mat))
  clusters = cutree(dist_clust,h = thresh)
  return(clusters)
}#end get_tn_clusters_SNPthresh

num_clusters_hclust <- function(cluster_list){
  clusters_id <- table(cluster_list) %>% subset(. >1) %>% names 
  cluster_str = subset(cluster_list,cluster_list %in% clusters_id) 
  no_cluster <- clusters_id %>% length
  no_cluster_ct <- cluster_str %>% length( )
  not_in_cluster_ct <- cluster_list %>% unlist %>% subset(!.%in% clusters_id) %>% length
  total_cluster_singleton_ct <- no_cluster + not_in_cluster_ct
  cluster_size_med <- table(cluster_str) %>% median
  cluster_size_mean <- table(cluster_str) %>% mean
  cluster_size_max <- table(cluster_str) %>% max
  cluster_size_range <- table(cluster_str) %>% range %>% paste0(collapse = "-")
  clustering_data <- cbind.data.frame(no_cluster,total_cluster_singleton_ct=total_cluster_singleton_ct,no_cluster_ct,not_in_cluster_ct=not_in_cluster_ct,cluster_size_med=cluster_size_med,cluster_size_mean=cluster_size_mean,cluster_size_max=cluster_size_max,cluster_size_range)
  return(clustering_data)
}

collapse_cl <- function(metadata, cluster_str, vars) {
  metadata %>% subset(rownames(.) %in% names(cluster_str)) -> meta_sub
  as.data.frame(cluster_str) %>% `colnames<-`("cluster") -> clusters_str_df
  clusters_str_df$Isolate <- rownames(clusters_str_df)
  left_join(meta_sub, clusters_str_df) -> meta_sub_cl
  res <- meta_sub_cl %>% group_by(across(c(vars, "cluster"))) %>% summarise(count = n(), isolates = paste0(Isolate, collapse = ","),.groups = 'drop')
  return(res)
}

expand_ids <- function(thresh, clr_collapse) {
  clr_collapse %>% subset(threshold == thresh) %>% select(isolates, cl_prof_id) %>% separate_rows(isolates, sep = ",") -> isolate_id_df
  isolate_id_str <- isolate_id_df$cl_prof_id %>% `names<-`(isolate_id_df$isolates) 
  return(isolate_id_str)
}

## -- Implementation --
thresholds = seq(from = 25, to = 1500, by = 25)
colnames(ng_meta_gene_grouped[,-1]) -> vars # variables of interest

# parallelized loop: 
clustering_results <- mclapply(thresholds, function(i) {
  clusters <- get_tn_clusters_SNPthresh(distmat, thresh = i)
  cluster_stats <- num_clusters_hclust(clusters)
  cluster_stats$threshold <- i
  cluster_collapse <- collapse_cl(ng_meta_gene_grouped, clusters, vars) %>% mutate(threshold = i)
  cluster_collapse <- cluster_collapse %>% group_by(cluster) %>% arrange(desc(count)) %>% mutate(cl_prof_id = paste0(cluster, ".", row_number())) %>% ungroup()
  prof_ids <- expand_ids(i, cluster_collapse)
  cluster_stats$distinct_profiles = nrow(cluster_collapse)
  res = list(cluster_str = clusters, cluster_stats = cluster_stats, cluster_collapse = cluster_collapse, profile_str = prof_ids)
  return(res)
}, mc.cores=8) %>% `names<-`(thresholds) 

lapply(clustering_results, FUN = function(x) {
  x[["cluster_stats"]] 
}) %>% do.call(rbind, .) -> cluster_stats_df
saveRDS(cluster_stats_df, "./data/phylogeny/source/distinct_profiles/cluster_stats.RDS")

lapply(clustering_results, FUN = function(x) {
  x[["cluster_collapse"]]
}) %>% do.call(rbind, .) -> cluster_collapse_df
saveRDS(cluster_collapse_df, "./data/phylogeny/source/distinct_profiles/cluster_collapse.RDS")

lapply(clustering_results, FUN = function(x) {
  x[["cluster_str"]] 
}) %>% do.call(cbind.data.frame, .) -> cluster_str_df
saveRDS(cluster_str_df, "./data/phylogeny/source/distinct_profiles/cluster_str.RDS")

lapply(clustering_results, FUN = function(x) {
  x[["profile_str"]] 
}) %>% do.call(cbind.data.frame, .) -> profile_id_df
saveRDS(profile_id_df, "./data/phylogeny/source/distinct_profiles/profile_id.RDS")

# save the whole list clr result
saveRDS(clustering_results, "./data/phylogeny/source/distinct_profiles/cluster_results_list.RDS")
saveRDS(list(script = "scripts/representative_isolates/03__distinct_profiles/02__cluster_gen.R", original_script = "../../2024_NG_enrichment_probe_creation/2025-02-22_distinct_profiles/scripts/cluster_gen.R", generated_at = Sys.time()), "./data/phylogeny/source/distinct_profiles/provenance.RDS")
