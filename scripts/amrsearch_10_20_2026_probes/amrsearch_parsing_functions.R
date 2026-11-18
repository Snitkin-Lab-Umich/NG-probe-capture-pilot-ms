# AMR-SEARCH Parsing Scripts
# Date: 09/29/26
# Author: Kyle Gontjes
get_amrsearch_data <- function(path){
  isolate <- basename(path) %>% gsub(".jsn","",.)
  amrsearch_data <- jsonlite::read_json(path)
  resistance_profile <- get_resistance_profile(amrsearch_data$resistanceProfile)
  resistance_profile$isolate <- isolate
  genotype_profile <- get_mutation_data(amrsearch_data$matches)
  genotype_profile$isolate <- isolate
  
  list(
       resistance_profile = resistance_profile,
       genotype_profile = genotype_profile)
}

get_resistance_profile <- function(resistance_profile){
  resistance_profile_data <- lapply(resistance_profile,FUN=function(x){
    antibiotic <- x$agent$name
    antibiotic_abv <-  x$agent$key
    antibiotic_type <-  x$agent$type
    state <- x$state 
    antibiotic_determinants_acquired <-lapply(x$determinants$acquired,FUN=function(y){y$gene}) %>% unlist %>% paste0(.,collapse=';')
    antibiotic_determinants_variant <- lapply(x$determinants$variants,FUN=function(y){paste0(y$gene,"_",y$variant)}) %>% unlist %>% paste0(.,collapse=';')
    antibiotic_determinants <- c(antibiotic_determinants_acquired,antibiotic_determinants_variant) %>% subset(.!='') %>% paste0(collapse=';')
    cbind.data.frame(antibiotic,antibiotic_abv,antibiotic_type,state,antibiotic_determinants,antibiotic_determinants_acquired,antibiotic_determinants_variant)
  }) %>% do.call(rbind,.)
  return(resistance_profile_data)
}

get_mutation_data <- function(matches){
  matches_data <- lapply(matches,FUN=function(x){
    # Gene level data
    gene_data <- data.frame(x[names(x)!='resistanceVariants']) 
    gene_data$entry <- "gene"
    gene_data$genotype_name <- gene_data$refId
    # Variant data
    resistance_variables <- x$resistanceVariants
    if(length(resistance_variables)>0){
      variant_data <- lapply(resistance_variables,FUN=function(y){
        dataset <- data.frame(y)
        
        dataset$entry <- 'variant'
        dataset$refId <- x$refId
        dataset$queryId <- x$queryId
        dataset$strand <- x$strand
        dataset$reversed <- x$reversed
        
        dataset
      }) %>% do.call(rbind.data.frame,.)
      
      variant_data$refId <-  x$refId
      variant_data$genotype_name <- paste0(variant_data$refId,"_",variant_data$name)
      final_data <-    bind_rows(gene_data,variant_data)
      final_data
    } else {
      gene_data
    }
  }) %>% do.call(bind_rows,.)
  return(matches_data)
}
