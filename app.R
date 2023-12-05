#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Author: Marianyela Petrizzelli
source(paste(getwd(), 'functions.R', sep = "/"))
load(paste(getwd(),'.RData', sep = "/"))


# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("OpenTarget - Network rappresentations"),
    sidebarLayout(
      sidebarPanel(
        # drop down menu 
        selectInput("data", "Select cell type interaction", 
                      choices = as.vector(lab_interaction)),
        # slider input for pvalue 
        selectInput("pval","p-value cutoff:", choices=set_pvalues, selected = 0.01),# width= "50%")
         
        selectInput("qval","q-value cutoff:", choices=set_qvalues, selected = 0.01),# width ="50%")
        
        sliderInput("ncat","Number of BF", min=1, max= 15, value = 5)
    ),
    # Show a plot of the generated distribution
    mainPanel(
        tabsetPanel(
          tabPanel(title = "GO terms",
                   visNetworkOutput("cnetGO", width="100%", height = "700px")),
          tabPanel(title = "Enrichment plot",
                   plotOutput("dotplot", width = "100%", height = "700px")),
          tabPanel(title = "Enrichment map",
                 plotOutput("emaplot", width = "100%", height = "700px"))
          )
)))
    
    


# Define server logic required to draw a histogram
server <- function(input, output) {
    output$intnet <- renderPlot
  
    output$cnetGO <- renderVisNetwork({
        row_indx = grep(input$data, lab_interaction, fixed = T)
        data_selection = selected_data(interaction_tab, tab_indx_lab, row_indx, cell_type)
        col_cel_type = color_sel_celltype(cell_type)
        g = g_cnet(kk[[paste0(input$data, "_", input$pval, "_", input$qval)]], 
                   data_selection, n_cat = input$ncat, col_cel_type)
        construct_complex_interactive_network(g)
    })
    output$dotplot <- renderPlot({
      dotplot(kk[[paste0(input$data, "_", input$pval, "_", input$qval)]], showCategory=input$ncat) + ggtitle("Top-contributing genes")})
    
    output$emaplot <- renderPlot({
      ego_pos2 <- pairwise_termsim(kk[[paste0(input$data, "_", input$pval, "_", input$qval)]], method="Wang", semData = d)
      emapplot(ego_pos2, showCategory = input$ncat, cex_label_category =0.9)})
    
}


# Run the application 
shinyApp(ui = ui, server = server)

