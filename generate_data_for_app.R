## Author: Marianyela Petrizzelli
## generate the data to be visualized:
## it may take some time 
source(paste(getwd(), 'functions.R', sep = "/"))

data_dir = ""
cell_cell_interaction_table_file = ""
set_pvalues = c(0.01, 0.02, 0.05)
set_qvalues = c(0.01, 0.02, 0.05)

build_interaction_networks(paste0(data_dir, '/', cell_cell_interaction_table_file, sep=""))
interaction_tab = load_data(data_dir = data_dir)
## select the cell types of interest 
cell_type = vec_cell_type(interaction_tab) # here we consider all that are present in the cell_cell_interaction_table_file
tab_indx_lab = indx_lab(cell_type)
lab_interaction = tab_indx_lab[, "lab"]

kk = list()
for(pvalue in set_pvalues){
  for(qvalue in set_qvalues){
    for (cc_int in lab_interaction){ 
      row_indx = grep(cc_int, lab_interaction, fixed = T)
      data_selection = selected_data(interaction_tab, tab_indx_lab, row_indx, cell_type)
      ego_enrich = enrichment_analysis(data_selection, 
                              pvaluecutoff = pvalue,
                              qvaluecutoff = qvalue)
      kk[[paste0(cc_int,"_", pvalue, "_", qvalue)]] = ego_enrich
  }}
}

d = godata(
  OrgDb = "org.Hs.eg.db",
  keytype = "SYMBOL",
  ont = "BP",
  computeIC = TRUE
)

save.image(paste0(data_dir, "/.RData", sep = ""))