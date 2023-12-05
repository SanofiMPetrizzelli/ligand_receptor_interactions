# A dynamical network representation of cell-cell communication strategies

This folder provides a set of codes to build and run a ShinyApp in R which allows for visualization of gene set enrichment (GSE) plots. 

The codes have been built to visualize GSE of ligand-receptor interactions among myeloid cells but can be easely updated to build the visualization for your own dataset. 

To this end, update the **generate_data_for_app.R** by setting your data directory (```data_dir``` variable). Set the ligand-receptor interaction file name (```cell_cell_interaction_table_file```, .csv format, comma separated)  with columns ```from``` and ```to``` indicating the interacting cell types, and ```ligand``` and ```receptor``` with the ligand from the first cell type and receptor from the second cell type under consideration.   

I recommend to run the script **generate_data_for_app.R** in order to generate the data first (it may require some time). It will generate an .RData object that is subsequently charged to the R environment for visualization by running the **app.R** script. 

