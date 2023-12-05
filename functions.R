## load library 
library(shiny)
library(igraph)
library(networkD3)
library(htmlwidgets)
library(enrichplot)
library(DT)
library("ggnewscale")
library(GOSemSim)
library(ggraph)
library(readr)
library(GO.db)
library(clusterProfiler)
library(visNetwork)

## here load all functions
indx_lab = function(cell_type){
  indx_lab = NULL
  for (ii in 1:length(cell_type)){
    x = cell_type[ii]
    for (jj in 1:length(cell_type)){
      y = cell_type[jj]
      lab = c(ii, jj, paste( x, " | ", y))
      indx_lab <- rbind(indx_lab, lab)
    }}
  colnames(indx_lab) <- c("idx", "idy", "lab")
  return(indx_lab)
}

load_data <- function(data_dir){
  edge_table <- read_csv(paste(data_dir, "edge_table.csv", sep="/"))
  node_table <- read_csv(paste(data_dir,"node_table.csv", sep="/"))
  node_table = as.data.frame(node_table)
  edge_table = as.data.frame(edge_table)
  rownames(node_table) <- node_table$name.copy
  ligand_receptor = lapply(edge_table$`shared name`, function(x){
    return(strsplit(x, split = '_'))
  })
  ligand_receptor = unlist(ligand_receptor, recursive = FALSE)
  ligand_receptor = do.call(rbind, ligand_receptor)
  edge_table$ligand <- ligand_receptor[,1]
  edge_table$receptor <- ligand_receptor[,2]
  edge_table$name <- edge_table$`shared name`
  
  # read the informations
  human_lr_pair <- read.delim(paste(data_dir, "human_lr_pair.txt", sep="/"))
  human_gene_info <- read.delim(paste(data_dir, "human_gene_info.txt", sep="/"))
  human_gene2ensembl <- read.delim(paste(data_dir, "human_gene2ensembl.txt", sep="/"))
  
  ### EXTRACT LR INFORMATION 
  lr_genes = c(unique(human_lr_pair$ligand_gene_symbol), unique(human_lr_pair$receptor_gene_symbol))
  ligand_receptor_description = c(rep('ligand', length(unique(human_lr_pair$ligand_gene_symbol))),
                                  rep('receptor', length(unique(human_lr_pair$receptor_gene_symbol))))
  lr_description = cbind(lr_genes, ligand_receptor_description)
  
  lr_desc_function = lapply(lr_genes, function(x){
    pos = match(x, human_gene_info$Symbol)
    return(human_gene_info[pos[1], c("description", "type_of_gene")])
  })
  lr_desc_function = do.call(rbind, lr_desc_function)
  lr_description = cbind(lr_description, lr_desc_function)
  rownames(lr_description) = lr_description$lr_genes
  
  interaction_df = edge_table
  in_db = lapply(1:dim(interaction_df)[1], function(x){
    g1_fn = match(interaction_df$ligand[x], lr_description$lr_genes)
    g2_fn = match(interaction_df$receptor[x], lr_description$lr_genes)
    if (is.na(g1_fn) | is.na(g2_fn)){
      interaction = 'L_OR_R_NOT_IN_DB'
    } else {
      if (lr_description[g1_fn,]$ligand_receptor_description==lr_description[g2_fn,]$ligand_receptor_description){
        interaction = 'L_L_OR_R_R'
      } else {
        if (lr_description[g1_fn,]$ligand_receptor_description == "receptor"){
          sel_rec = match(lr_description[g1_fn, ]$lr_genes, human_lr_pair$receptor_gene_symbol)
          if (lr_description[g2_fn,]$lr_genes == human_lr_pair$ligand_gene_symbol[sel_rec]){
            interaction = "FOUND"
            evidence = human_lr_pair$evidence[sel_rec]
          } else {
            interaction = "NOT_FOUND" 
          }
        } else {
          sel_rec = match(lr_description[g2_fn, ]$lr_genes, human_lr_pair$receptor_gene_symbol)
          if (lr_description[g1_fn,]$lr_genes == human_lr_pair$ligand_gene_symbol[sel_rec]){
            interaction = "FOUND"
          } else {
            interaction = "NOT_FOUND" 
          }
        }
      }
    }
    return(interaction)
  })
  interaction_df = cbind(interaction_df, in_db= unlist(in_db))
  return(interaction_df)
}

vec_cell_type <- function(interaction_df){
  return(unique(c(interaction_df$from, interaction_df$to)))
}

selected_data <- function(interaction_df, tab_indx_lab, row_indx, cell_type){
  ii = as.numeric(tab_indx_lab[row_indx,1])
  jj = as.numeric(tab_indx_lab[row_indx,2])
  x = cell_type[ii]
  y = cell_type[jj]
  df_x = interaction_df[which(interaction_df$from==x), ]  
  df_x_y = df_x[which(df_x$to==y), ]
  return(df_x_y)
}


enrichment_analysis <- function(df_x_y, pvaluecutoff = 0.01,
                                qvaluecutoff = 0.01){
  organism = "org.Hs.eg.db"
  df_x_y = cbind(df_x_y$ligand, df_x_y$receptor, df_x_y) 
  g = graph_from_data_frame(d=df_x_y, directed = TRUE)
  ego = enrichGO(V(g)$name, OrgDb = organism, 
                 ont="ALL", keyType = "SYMBOL",
                 pAdjustMethod = "BH", #universe = geneList,
                 pvalueCutoff = pvaluecutoff, 
                 qvalueCutoff = qvaluecutoff)
  return(ego)
}

universe_geneList <- function(interaction_df){
  return(unique(c(interaction_df$source, interaction_df$target)))
}

extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}

list2df <- function(inputList) {
  ldf <- lapply(1:length(inputList), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  do.call('rbind', ldf)
}

update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  return(n)
}

g_cnet <- function(ego, df_x_y, n_cat = 5, col_cel_type){
  geom_edge <- geom_edge_link
  geneSets <- extract_geneSets(ego, n_cat)
  g <- list2graph(geneSets)
  edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  V(g)$color <- "#B3B3B3"
  V(g)$color[1:n] <- "#E5C494"
  V(g)$color[V(g)$name %in% df_x_y$ligand] = col_cel_type[unique(df_x_y$from)]
  V(g)$color[V(g)$name %in% df_x_y$receptor] = col_cel_type[unique(df_x_y$to)]
  V(g)$category <- "BF"
  V(g)$category[V(g)$name %in% df_x_y$ligand] <- "from_gene"
  V(g)$category[V(g)$name %in% df_x_y$receptor] <- "to_gene"
  return(g)
}

color_sel_celltype <- function(cell_type){
  col_cel_type <- c("lightpink", "azure",
                    "lemonchiffon", "lightblue",
                    "lightsalmon", "mistyrose", "navajowhite",
                    "plum1", "seagreen1")
  names(col_cel_type) <- cell_type
  return(col_cel_type)
}

construct_complex_interactive_network <- function(g){
  node_table = as.data.frame(get.vertex.attribute(g))
  node_table$id = node_table$name
  node_table$label = node_table$name
  node_table$value = node_table$size
  node_table$group = node_table$category
  edge_table = as.data.frame(get.edgelist(g))
  colnames(edge_table) = c("from", "to")
  visNetwork(node_table, edge_table)  %>% visIgraphLayout()
}


Symbol_entreId <- function(geneSymbol){
  geneList_conversion = bitr(geneSymbol, fromType="SYMBOL", 
                             toType="ENTREZID", OrgDb="org.Hs.eg.db")
}
