####################################################################################

# CNV shiny app for Middle Eastern Genomes Project                                 #

# R Version 4.1.0                                                                  #

# Authors: Luke Peacock m1606864@sgul.ac.uk & Alan Pittman apittman@sgul.ac.uk     #

# RStudio version 1.4.1717                                                         #

####################################################################################

#library(rsconnect)
#rsconnect::deployApp('S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/R/shiny/app_final')

#file.copy("S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db", "S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/R/shiny/app_final/CNV_data.db")
#file.copy("S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/ROH_calls/bcftools_roh_output", "S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/R/shiny/app_final/bcftools_roh_output")

#library(BiocManager)
#repos <- BiocManager::repositories()
#getOption(toString(repos))

library(DBI)
library(RSQLite)
library(DT)
library(RCircos)
library(biomaRt)
library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)

#library(BiocManager)
#repos <- BiocManager::repositories()
#options()



CNV_db <- dbConnect(RSQLite::SQLite(), 'CNV_data.db')
ID_group <- dbGetQuery(CNV_db, 'SELECT distinct ID from fullcombined')
ID_roh <- dbGetQuery(CNV_db, 'SELECT distinct sample_ID from rohcombined')
PI_ID <- dbGetQuery(CNV_db, 'SELECT DISTINCT PI FROM fullcombined INNER JOIN sample_info ON fullcombined.ID = sample_info.sampleID')
proj_ID <- dbGetQuery(CNV_db, 'SELECT DISTINCT projectID FROM fullcombined INNER JOIN sample_info ON fullcombined.ID = sample_info.sampleID')

#dbDisconnect(CNV_db)

unlink("www/RCircos_*.pdf")
unlink("*.csv")

ui <- dashboardPage(skin = "blue",
                    dashboardHeader(title = "MEGP Browser"),
                    
                    dashboardSidebar(
                      sidebarMenu(menuItem("How to use", tabName = "help", icon = icon("question-circle")),
                                  menuItem("CNV Workflow", tabName = "workflow", icon = icon("dna")),
                                  menuItem("How to interpret tables", tabName = "instruct_tables", icon = icon("question-circle")),
                                  menuItem("Download CSVs of data", tabName = "csvs", icon = icon("download")),
                                  menuItem("Search by Sample", tabName = "search", icon = icon("search")),
                                  menuItem("Search CNVs in ROH", tabName = "roh", icon = icon("search"))
                                  )
                    ),
                                     
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "help",
                                fluidRow(
                                  box(width = 12,
                                      div(class="centre", h1("Middle Eastern Genomes Project", align = "middle"),
                                          h2("Variant Browser", align = "middle"),
                                          br(),
                                          (p(HTML('<p><center><img src="sgul_logo.jpg" style="width:290px;height:150px""/></center></p>'))),
                                          h3(tags$a(href="https://www.genomemed.org/megp", "More Information About the MEGP"), align = "middle")),
                                      br(),
                                      div(align = "middle",
                                        h3("The navigation panel to the left allows you to view the different tabs available:", align = "centre"),
                                        br(),
                                        h4("Tab 1) 'How to Use' Page - this page!"),
                                        h4("Tab 2) 'Search by Sample' Page - this allows you to select or type a sample ID to view its CNV calls. This outputs an interactive CNV table and a Circos plot which can be downloaded.",
                                        h4("Tab 3) 'Search by Gene' Name - this allows you to search the whole CNV database for the dataset. This outputs an interactive CNV table - showing the sample IDs - and a Circos plot which can be downloaded."))
                                  )
                                )
                        )),
                        
                        tabItem(tabName = "csvs",
                                fluidRow(
                                  box(width = 12,
                                      h2("Download full data set:", align = "middle"),
                                      downloadButton("full_set_query", "Download full dataset in csv format"),
                                      downloadButton("roh_set_query", "Download CNVs in ROH dataset in csv format"),
                                      downloadButton("roh_regions", "Download raw ROH data for all samples in bed format"),
                                      downloadButton("OMIM_all", "Download OMIM annotated CNVs - full dataset")
                                  ),
                                  box(width = 12,
                                    h2("Select your project ID or Principle Investigator to download a subset of CNV data in csv format:", align = "middle")),
                                  box(width = 12,
                                      h4("Search by Principle Investigator (PI):", align = "middle"),
                                      selectInput(
                                    'PI',
                                    'Select your PI:',
                                    PI_ID$PI,
                                    selected = NULL),
                                    downloadButton("csv_query_PI", "Download 'PI' csv from full dataset"),
                                    h5("Expand the table below to see the sample query output:", align = "centre")),
                                    box(width = 12,
                                        shinycssloaders::withSpinner(DT::dataTableOutput("PI_tab")), collapsible = TRUE, collapsed = TRUE),
                                    box(width = 12,
                                    h4("Search by Project ID:", align = "middle"),
                                    selectInput(
                                      'proj',
                                      'Select your project ID:',
                                      proj_ID$projectID,
                                      selected = NULL),
                                    h5("Query from full dataset:", align = "middle"),
                                    downloadButton("csv_query_proj", "Download 'Project ID' csv from full dataset"),
                                    h5("Expand the table below to see the sample query output:", align = "centre")),
                                    box(width = 12,
                                        shinycssloaders::withSpinner(DT::dataTableOutput("proj_tab")), collapsible = TRUE, collapsed = TRUE),
                                  
                                  
                                  box(width = 12,
                                      h2("Select your project ID or Principle Investigator to download a subset/all of ROH data in csv format:", align = "middle"),
                                      h5("These ROHs have been calculated based on VCF files using bcftools.", align = "middle")),

                                  box(width = 12,
                                      h4("Search by Principle Investigator (PI):", align = "middle"),
                                      selectInput(
                                        'PI_roh',
                                        'Select your PI:',
                                        PI_ID$PI,
                                        selected = NULL),
                                      downloadButton("csv_query_PI_roh", "Download 'PI' csv from full dataset"),
                                      h5("Expand the table below to see the sample query output:", align = "centre")),
                                  
                                  box(width = 12,
                                      shinycssloaders::withSpinner(DT::dataTableOutput("PI_tab_roh")), collapsible = TRUE, collapsed = TRUE),
                                  
                                  box(width = 12,
                                      h4("Search by Project ID:", align = "middle"),
                                      selectInput(
                                        'proj_roh',
                                        'Select your project ID:',
                                        proj_ID$projectID,
                                        selected = NULL),
                                      h5("Query from full dataset:", align = "middle"),
                                      downloadButton("csv_query_proj_roh", "Download 'Project ID' csv from full dataset"),
                                      h5("Expand the table below to see the sample query output:", align = "centre")),
                                  
                                  box(width = 12,
                                      shinycssloaders::withSpinner(DT::dataTableOutput("proj_tab_roh")), collapsible = TRUE, collapsed = TRUE),
                                  
                                  
                                  )
                                  ),
                        
                        tabItem(tabName = "instruct_tables",
                                fluidRow(
                                  box(width = 12,
                                      div(class="centre", h2("How to interpret the CNV tables generated", align = "middle"),
                                          br(),
                                          (p(HTML('<p><center><img src="table_instructions.png" style="width:auto;height:auto"/></center></p>'))),
                                          br()
                                          ),
                                      div(align = "middle",
                                          h3("Description of each column:", align = "centre"),
                                          br(),
                                          HTML('<p><b>Type)</b> Shows the type of Copy Number Variant - i.e. duplication/deletion.</p>'),
                                          HTML('<p><b>nexons)</b> Shows the number of exons that this CNV spans.</p>',
                                          HTML('<p><b>Start)</b> The start point of the CNV - chromosome coordinates.</p>'),
                                          HTML('<p><b>End)</b> The end point of the CNV - chromosome coordinates.</p>'),
                                          HTML('<p><b>chromosome)</b> The chromosome containing the CNV.</p>'),
                                          HTML('<p><b>BF)</b> Bayes Factor: Can be interpreted as follows (note that both positive and negative values can be significant!):</p>'),
                                          HTML('<p><b>REF:</b> Oct2017 "Sociological Methods & Research" 48(3):004912411772971</>'),
                                          HTML('<p><center><img src="BF_help.png" style="width:300px;height:auto""/></center></p>'),
                                          HTML('<p><b>reads.ratio)</b> The ratio of the expected read counts vs the observed read counts. For example, 1.0 suggests wildtype; 0.0 suggests deletion of this region; 2.0 suggests duplication of this region.</p>'),
                                          HTML('<p><b>Conrad.hg19)</b> Known <i>common</i> CNVs in the detected CNV region - from the Conrad et al 2010 Nature paper. </p>'),
                                          HTML('<p><b>exons.hg19)</b> hg19 reference genome exon names in the detected CNV region.</p>'),
                                          )
                                      )
                                  )
                                )),                        
                                                
                        tabItem(tabName = "roh",
                                fluidRow(shinyjs::useShinyjs(), shinyalert::useShinyalert(),
                                         box(width = 6,
                                             div(class="centre", h4("CNVs in Regions of Homozygosity", align = "centre"),
                                                 br(),
                                                 h5("Enter the sample ID to view the SQL table output", align = "centre")
                                             )
                                         ),
                                         box(selectInput(
                                           'sample_roh',
                                           'Select your sample ID:',
                                           ID_roh$sample_ID,
                                           selected = NULL)),
                                         box(width = 12,
                                             collapsible = TRUE,
                                             h3("Interactive Table of CNVs in ROH", align = "centre"),
                                             downloadButton("download_roh_tab", "Download csv of table"),
                                             shinycssloaders::withSpinner(DT::dataTableOutput("roh_out"))
                                             # downloadButton(export_data, "Download")
                                         )
                                         #box()
                                         )),
                        
                        tabItem(tabName = "search",
                                fluidRow(shinyjs::useShinyjs(),
                                  box(width = 6,
                                      h5("NOTE: other sidebar tabs will not work whilst the this tab is generating a circos plot", align = "centre"),
                                     div(class="centre", h3("Search for all CNVs in a sample", align = "centre"),
                                  br(),
                                   h5("Enter the sample ID to view the SQL table output", align = "centre")
                                             )),
                                  box(selectInput(
                                    'sample',
                                    'Select your sample ID:',
                                    ID_group$ID,
                                    selected = NULL)),
                                  box(width = 12,
                                      collapsible = TRUE,
                                      h3("Interactive Table of CNVs", align = "centre"),
                                      downloadButton("download_tab", "Download csv of table"),
                                      shinycssloaders::withSpinner(DT::dataTableOutput("table_out")),
                                      ),
                                  
                                  box(h3("Generate a Circos Plot:", align = "centre"),
                                    actionButton("generate", label = "Click to Generate Circos Plot")
                                    ),
                                  box(plotOutput("circos"), width = 1, height = 0, collapsed = TRUE
                                    ),
                                )
                    
                       
                          
                        )
                      )
                    )
                    )



server <- function(input, output, session) {
  
  shinyalert::shinyalert(title = "Downloaded Files", text = paste("circos plots and associated gene files downloaded within this session will be automatically deleted from the R directory next time this app is launched. Make sure to save any files you want to keep in a different directory!", sep = ""), type = "warning", timer = 70000)
  
  
  output$roh_regions <- downloadHandler(filename = "all_ROHs.bed", content = "bcftools_roh_output")
  
  output$full_set_query <- downloadHandler(
    filename = function() {paste("full_CNV_dataset_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(data_A <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM fullcombined INNER JOIN sample_info ON fullcombined.ID = sample_info.sampleID INNER JOIN frequencies ON fullcombined.ID_CNV_id_type = frequencies.CNV_ID'
    ), file, row.names = FALSE)}
  )
  
  output$roh_set_query <- downloadHandler(
    filename = function() {paste("full_CNV_dataset_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(data_B <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM rohcombined INNER JOIN sample_info ON rohcombined.ROH_sample_ID = sample_info.sampleID INNER JOIN frequencies ON rohcombined.ROH_sample_ID_CNV_ID_CNV_type = frequencies.CNV_ID'
    ), file, row.names = FALSE)}
  )
  
  output$OMIM_all <- downloadHandler(
    filename = function() {paste("OMIM_full_CNV_dataset_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(data_C <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM OMIM_full INNER JOIN sample_info ON OMIM_full.sample_ID = sample_info.sampleID INNER JOIN frequencies ON OMIM_full.CNV_ID = frequencies.CNV_ID'
    ), file, row.names = FALSE)}
  )
  
  output$csv_query_PI <- downloadHandler(
    filename = function() {paste(input$PI, "_all_CNVs_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(data_t <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM fullcombined INNER JOIN sample_info ON fullcombined.ID = sample_info.sampleID WHERE PI = ?',
      params = input$PI, 
    ), file, row.names = FALSE)}
  )
  
  output$csv_query_proj <- downloadHandler(
    filename = function() {paste(input$proj, "_all_CNVs_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(data_t <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM fullcombined INNER JOIN sample_info ON fullcombined.ID = sample_info.sampleID WHERE projectID = ?',
      params = input$proj, 
    ), file, row.names = FALSE)}
  )

  output$PI_tab <- DT::renderDataTable(
    #colnames = c('CNV Type', 'Number of exons', 'Start', 'End', 'Chromosome', 'BF', 'reads.ratio', 'Conrad.hg19', 'exons.hg19'),
    data_tab <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM fullcombined INNER JOIN sample_info ON fullcombined.ID = sample_info.sampleID WHERE PI = ?',
      params = input$PI, 
    )
  )
  
  output$proj_tab <- DT::renderDataTable(
    #colnames = c('CNV Type', 'Number of exons', 'Start', 'End', 'Chromosome', 'BF', 'reads.ratio', 'Conrad.hg19', 'exons.hg19'),
    data_tab <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM fullcombined INNER JOIN sample_info ON fullcombined.ID = sample_info.sampleID WHERE projectID = ?',
      params = input$proj, 
    )
  )
  
  
  
  output$csv_query_PI_roh <- downloadHandler(
    filename = function() {paste(input$PI_roh, "_in_ROH_CNVs_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(data_t <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM rohcombined INNER JOIN sample_info ON rohcombined.sample_ID = sample_info.sampleID WHERE PI = ?',
      params = input$PI_roh, 
    ), file, row.names = FALSE)}
  )
  
  output$csv_query_proj_roh <- downloadHandler(
    filename = function() {paste(input$proj_roh, "_in_ROH_CNVs_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(data_t <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM rohcombined INNER JOIN sample_info ON rohcombined.sample_ID = sample_info.sampleID WHERE projectID = ?',
      params = input$proj_roh, 
    ), file, row.names = FALSE)}
  )
  
  output$PI_tab_roh <- DT::renderDataTable(
    #colnames = c('CNV Type', 'Number of exons', 'Start', 'End', 'Chromosome', 'BF', 'reads.ratio', 'Conrad.hg19', 'exons.hg19'),
    data_tab <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM rohcombined INNER JOIN sample_info ON rohcombined.sample_ID = sample_info.sampleID WHERE PI = ?',
      params = input$PI_roh, 
    )
  )
  
  output$proj_tab_roh <- DT::renderDataTable(
    #colnames = c('CNV Type', 'Number of exons', 'Start', 'End', 'Chromosome', 'BF', 'reads.ratio', 'Conrad.hg19', 'exons.hg19'),
    data_tab <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM rohcombined INNER JOIN sample_info ON rohcombined.sample_ID = sample_info.sampleID WHERE projectID = ?',
      params = input$proj_roh, 
    )
  )
  
  
  
  output$table_out <- DT::renderDataTable(
    #colnames = c('CNV Type', 'Number of exons', 'Start', 'End', 'Chromosome', 'BF', 'reads.ratio', 'Conrad.hg19', 'exons.hg19'),
    data_tab <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT type, nexons, start, end, chromosome, BF, \"reads.ratio\", \"Conrad.hg19\", \"exons.hg19\" FROM fullcombined WHERE ID = ?',
      params = input$sample
    )
  )
  
  output$download_tab <- downloadHandler(
    filename = function() {paste(input$sample, "_CNVs_", Sys.Date(), ".csv", sep = "")},
    content = function(file) {write.csv(data_t <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM fullcombined WHERE ID = ?',
      params = input$sample
    ), file, row.names = FALSE)}
  )
  
  
  output$roh_out <- DT::renderDataTable(
    #colnames = c('CNV Type', 'Number of exons', 'Start', 'End', 'Chromosome', 'BF', 'reads.ratio', 'Conrad.hg19', 'exons.hg19'),
    data_roh <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT CNV_type, nexons, CNV_start, CNV_end, CNV_chromosome, BF, \"reads.ratio\", \"Conrad.hg19\", \"exons.hg19\" FROM rohcombined WHERE sample_ID = ?',
      params = input$sample_roh
    )
  )
  
  output$download_roh_tab <- downloadHandler(
    filename = function() {paste(input$sample, "_CNVs_in_ROH_", Sys.Date(), ".csv")},
    content = function(file) {write.csv(data_roh <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM rohcombined WHERE sample_ID = ?',
      params = input$sample_roh
    ), file, row.names = FALSE)}
  )
  
  output$proj_tab <- DT::renderDataTable(
    #colnames = c('CNV Type', 'Number of exons', 'Start', 'End', 'Chromosome', 'BF', 'reads.ratio', 'Conrad.hg19', 'exons.hg19'),
    data_tab <- dbGetQuery(
      conn = CNV_db,
      statement = 'SELECT * FROM fullcombined INNER JOIN sample_info ON fullcombined.ID = sample_info.sampleID WHERE projectID = ?',
      params = input$proj, 
    )
  )

  observeEvent(input$generate, {
    
  output$circos <- renderPlot({
    
    shinyalert::shinyalert(title = "Request received", text = paste("Your circos plot request has been received on the server side. Please do not press button again until the next pop-up window appears. This and other tabs will temporarily not work as the server generates the circos plot!", sep = ""), type = "success", timer = 10000)
    
    ID_sample <- input$sample
    
    #load data
    CNV_db <- dbConnect(RSQLite::SQLite(), 'CNV_data.db')
    data_sam <- dbGetQuery(conn = CNV_db,
                           statement = 'SELECT chromosome, start, end, type FROM [fullcombined] WHERE ID = ?',
                           params = ID_sample)
    dbDisconnect(CNV_db)
    
    #make del and dup data frames
    data_sam["height"] <- 1
    data_del <- data_sam[!grepl("duplication", data_sam$type),]
    data_dup <- data_sam[!grepl("deletion", data_sam$type),]
    
    
    #biomart: genes
    #dataset: hsapiens_gene_ensembl
    #version: GRCh37.p13
    #attributes <- listAttributes(mart)
    #filters <- listFilters(mart)
    
    #connect to biomart
    
    mart <- useEnsembl(biomart = "genes",
                                         dataset = "hsapiens_gene_ensembl",
                                         GRCh = "37")
    
    #make data frames & fetch genes and positions for DELETIONS
    for (CNV_del in 1:nrow(data_del)) {
      temp <- paste(data_del$chromosome, ":", data_del$start, ":", data_del$end, sep = "")
    }
    del_pos <- data.frame(positions = c(temp))
    
    for (del in 1:nrow(del_pos)) {
      del_genelist <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name", "phenotype_description"),
                            filters = c("chromosomal_region"),
                            values = c(del_pos$positions),
                            mart = mart)
    }
    
    #make data frames & fetch genes and positions for DUPLICATIONS
    for (CNV_dup in 1:nrow(data_dup)) {
      temp <- paste(data_dup$chromosome, ":", data_dup$start, ":", data_dup$end, sep = "")
    }
    dup_pos <- data.frame(positions = c(temp))
    
    for (del in 1:nrow(del_pos)) {
      dup_genelist <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name", "phenotype_description"),
                            filters = c("chromosomal_region"),
                            values = c(dup_pos$positions),
                            mart = mart)
    }
    
    
    #set column 1 to 'chr'x
    data_del <- transform(data_del, chromosome = sprintf('chr%s', chromosome))
    data_dup <- transform(data_dup, chromosome = sprintf('chr%s', chromosome))
    del_genelist <- transform(del_genelist, chromosome_name = sprintf('chr%s', chromosome_name))
    dup_genelist <- transform(dup_genelist, chromosome_name = sprintf('chr%s', chromosome_name))
    
    #set column orders and names
    dup_genelist = dup_genelist[ , c("chromosome_name", "start_position", "end_position", "external_gene_name", "phenotype_description")]    
    del_genelist = del_genelist[ , c("chromosome_name", "start_position", "end_position", "external_gene_name", "phenotype_description")]    
    
    #data_sam=data_sam[-c(2:72),]
    #data_sam = data_sam [-1,]
    #names(data_sam) <- NULL
    
    #load cytoband data
    data("UCSC.HG19.Human.CytoBandIdeogram")
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
    chr.exclude <- c("chrX", "chrY")
    
    RCircos.Set.Core.Components(cyto.info, chr.exclude,
                                tracks.inside=9, tracks.outside=3);
    
    rcircos.params <- RCircos.Get.Plot.Parameters()
    rcircos.cyto <- RCircos.Get.Plot.Ideogram()
    rcircos.position <- RCircos.Get.Plot.Positions()
    RCircos.List.Plot.Parameters()
    
    
    out.file <- paste("www/RCircos_", "sample", ".pdf", sep = "")
    #out.file <- "C:/Users/Luke Peacock/Documents/summer_project/R/shiny/app_final/www/RCircos_test.pdf"
    pdf(file=out.file, height=8, width=8, compress=TRUE)
    #png(file=out.file, height=8, width=8, compress=TRUE)
    
    #recordPlot()
    
    RCircos.Set.Plot.Area()
    
    
    RCircos.Chromosome.Ideogram.Plot()
    
    #gene labels for DELETIONS
    #data(del_genelist)
    name.col <- 4
    side <- "in"
    track.num <- 2
    RCircos.Gene.Connector.Plot(del_genelist, track.num, side)
    track.num <- 3
    RCircos.Gene.Name.Plot(del_genelist, name.col,track.num, side)
    
    RCircos.Chromosome.Ideogram.Plot()
    
    #gene labels for DUPLICATIONS
    #data(del_genelist)
    name.col <- 4
    side <- "out"
    track.num <- 2
    RCircos.Gene.Connector.Plot(dup_genelist, track.num, side)
    track.num <- 3
    RCircos.Gene.Name.Plot(dup_genelist, name.col,track.num, side)
    
    RCircos.Chromosome.Ideogram.Plot()
    
    # OLD
    # data(RCircos.Gene.Label.Data)
    # name.col <- 4
    # side <- "in"
    # track.num <- 2
    # RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
    # track.num <- 3
    # RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col,track.num, side)
    # 
    # RCircos.Chromosome.Ideogram.Plot()
    
    
    
    #data(data_sam)
    #PLOT deletions (red)
    RCircos.Polygon.Plot(polygon.data=data_del, data.col=5, track.num=1, side="in", is.sorted = TRUE, polygon.col = "#E10600", border.col = '#E10600')
    RCircos.Polygon.Plot(polygon.data=data_dup, data.col=5, track.num=1, side="out", is.sorted = TRUE, polygon.col = "#66FF00", border.col = '#66FF00')
    
    # close graphic device
    dev.off()
    
    file_name <- paste("www/RCircos_",input$sample,".pdf", sep = "")
    file.copy("www/RCircos_sample.pdf", file_name)
    
    dup_filename <- paste(input$sample, "_genes_with_duplications.csv", sep = "")
    write.csv(dup_genelist, dup_filename, row.names = FALSE)
    del_filename <- paste(input$sample, "_genes_with_deletions.csv", sep = "")
    write.csv(del_genelist, del_filename, row.names = FALSE)
    
    #show_pdf <- system2('open', args = paste("www/RCircos_", input$sample, ".pdf"), wait = FALSE)
    shinyalert::shinyalert(title = "Circos Plot Information", text = paste("Your circos plot for sample ", input$sample, " is ready for download and will open on exiting this pop-up window.", sep = ""), type = "info", callbackR = function() {
      system2('open', args = paste("www/RCircos_", input$sample, ".pdf", sep = ""), wait = FALSE)
      shell.exec(paste("", input$sample, "_genes_with_duplications.csv", sep = ""))
      shell.exec(paste("", input$sample, "_genes_with_deletions.csv", sep = ""))
      })
    #system2('open', args = paste("www/RCircos_", input$sample, ".pdf", sep = ""), wait = FALSE)

  })
  })
}

shinyApp(ui, server)

