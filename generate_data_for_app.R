## save plots so that the algorithm is faster 
source(paste(getwd(), 'functions.R', sep = "/"))

interaction_tab = load_data(data_dir = "/cloud-home/I0558001/Signaling")
cell_type = vec_cell_type(interaction_tab)[c(4, 7, 23, 26, 24, 6, 22, 30, 9)]
tab_indx_lab = indx_lab(cell_type)
lab_interaction = tab_indx_lab[, "lab"]

kk = list()
for(pvalue in c(0.05)){
  for(qvalue in c(0.01, 0.02, 0.05)){
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
