source(file = "library_functions.r")
source(file = "SignalingProfiler.R")

##### UI #####

ui <- fluidPage(
  title = "MultiOmicsXplorer",
  useShinyjs(),
  tags$head(
    tags$title("MultiOmicsXplorer"),
    tags$style(HTML("
      body {background-color: #F5F5F5;} 
      .logo-container {text-align: center; margin-top: 20px; margin-bottom: 5px;} 
      .navbar-default {background-color: #FFFFFF; border: none; border-radius: 0;}
      .navbar-default .navbar-nav {
        float: none;
        display: flex;
        justify-content: center;
      }
      .navbar-default .navbar-nav > li > a {
        color: #333333;
        font-weight: bold;
      } 
      .navbar-default .navbar-nav > li > a:hover {
        background-color: #EEEEEE; 
        color: #33478e;
      } 
      .navbar-default .navbar-nav > li > a .fa-home {
        margin-right: 5px;
      } 
      
      .plot-container { max-width: 800px; margin: 0 auto; } 
          .gene-selection-box {
            padding: 10px 20px;  
            margin: 10px 0; 
            border-radius: 8px;  
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);  
            text-align: center; 
            font-size: 14px;  
            cursor: pointer;  
            border: 2px solid #D56E2A;  
            transition: background-color 0.3s ease, transform 0.3s ease, border-color 0.3s ease; 
          }
          
      .gene-selection-box:hover {
        background-color: #D56E2A;
        border-color: #D56E2A;  
      }
      
      .gene-selection-box:active {
        background-color: #FFB74D;  
        transform: scale(0.95);  
        border-color: #FFB74D; 
      }
      
     
      .gene-selection-box.clicked {
        animation: resetColor 1s forwards;  
        }
      
      @keyframes resetColor {
        0% {
          background-color: #FFB74D; 
          }
        100% {
          background-color: #D56E2A; 
        }
      }
      
       .tooltip-wrapper {
      position: relative;
      display: inline-block;
    }

    .tooltip-wrapper .tooltip-text {
      visibility: hidden;
      width: 240px;
      background-color: #333;
      color: #fff;
      text-align: left;
      border-radius: 5px;
      padding: 8px;
      position: absolute;
      z-index: 1;
      bottom: 125%; 
      left: 50%;
      transform: translateX(-50%);
      opacity: 0;
      transition: opacity 0.2s ease-in-out;
      font-size: 0.9em;
    }

    .tooltip-wrapper:hover .tooltip-text {
      visibility: visible;
      opacity: 1;
    }

    .tooltip-text::after {
      content: '';
      position: absolute;
      top: 100%;
      left: 50%;
      margin-left: -5px;
      border-width: 5px;
      border-style: solid;
      border-color: #333 transparent transparent transparent;
    }
  
    "))
  ),
  
  div(class = "logo-container",
      img(src = "logo.png", height = "200px", class = "logo")  
  ),
  
  
  #### Navbar ####
  navbarPage(
    title = "",
    tabPanel(value = "home", icon("home"),  
             fluidRow(
               column(width = 6, 
                      h2("Welcome to MultiOmicsXplorer!", style = "font-size: 2.5em; font-weight: bold; color: #1A2448;"), 
                      h3("This tool analyzes phosphoproteomic, proteomic, and transcriptomic data 
                               from various types of tumors at different stages, allowing comparisons 
                             between tumors and stages.", style = "font-size: 1.5em; color: #1A2448;"),
                      p(),
                      h3("In the ", span("OncoXplorer-CPTAC Data Analysis", style = "font-weight: bold;"), 
                         " section, you can compare tumor data from the ", 
                         a("CPTAC portal", href = "https://proteomics.cancer.gov/programs/cptac", target = "_blank", style = "font-weight: bold; color: #5a9bd5;text-decoration: underline;"), 
                         " across different stages and tumor types.", 
                         style = "font-size: 1.5em; color: #1A2448;"),
                      p(),
                      h3("The ", span("Extract Protein Activity from Your Data", style = "font-weight: bold;"), 
                         " section allows you to use SignalingProfiler to extract protein activity uploading your data.", 
                         style = "font-size: 1.5em; color: #1A2448;"),
                      h3("Click here for the ", 
                         a("tutorial", href = "tutorial.html", target = "_blank", style = "font-weight: bold; color: #5a9bd5;text-decoration: underline;"), 
                         " to learn how to use the app.", 
                         style = "font-size: 1.5em; color: #1A2448;"),
                      align = "left"
                    
               ),
               #column(width = ),  
               column(
                 width = 6,
                 img(
                   src = "MultiOmicsXplorer_home.png",
                   style = "width: 100%; height: auto; display: block; margin: 0 auto; max-width: none;"
                 )
             )
    )),   
    #### OncoXplorer - CPTAC data analysis ####
    tabPanel("OncoXplorer - CPTAC data analysis",
             tabsetPanel(
               tabPanel("About the data",
                        div(
                          style = "padding: 20px;",
                          fluidRow(
                            column(
                              width = 7,  
                              h2(strong("CPTAC data")),  
                              h3("Data Processing and Harmonization"),
                              p("We use data that have been downloaded using the tutorial described in the following paper:"),
                              tags$a(href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8022323/", 
                                     "Simplified and Unified Access to Cancer Proteogenomic Data", 
                                     target = "_blank"),
                              p("We collected clinical data - used to obtain information about the tumors stages, and raw 
                                transcriptomics proteomics and phosphoproteomics data."),
                              h3("Pre-processing Steps"),
                              p("The raw data underwent several filtering and transformation steps to ensure consistency and 
                                usability:"),
                              tags$ul(
                                tags$li(strong("Transcriptomics: "), "Filtered out rows with more than 80% zero values and removed 
                                        non-coding elements."),
                                tags$li(strong("Proteomics: "), "Updated protein identifiers to match Uniprot ID."),
                                tags$li(strong("Phosphoproteomics: "), "Updated Uniprot ID and sequence window based on the 
                                        phosphorylated residue.")
                              ),
                              h3("Harmonization Process"),
                              p("To standardize the data across different omics layers, we applied the following transformations:"),
                              tags$ul(
                                tags$li("Missing values (NA) in proteomics and phosphoproteomics data were imputed."),
                                tags$li("Z-score normalization was applied to all three omics datasets to facilitate comparisons.")
                              ),
                            ),  
                            
                            column(
                              width = 5,  
                              div(
                                style = "text-align: center; margin-top: 20px;", 
                                img(src = "sample_omics.png", style = "max-width: 500px; width: 105%; height: auto;")
                              ),
                              p(),
                              p(),
                              p(
                                "If you need to harmonize your data, you can easily access the ", 
                                a("PatientProfiler", href = "https://github.com/SaccoPerfettoLab/PatientProfiler.git", target = "_blank"), 
                                " package", 
                                a("(Lombardi et al., 2025)", href = "https://www.biorxiv.org/content/10.1101/2025.01.31.635886v1", target = "_blank"), 
                                
                                style = "font-size: 1em;font-weight: bold; color: #1A2448; text-align: center;"
                              )
                            )
                          )
                        )
               ),
               #### Analyte abundance panel ####
               tabPanel("Explore Analyte Abundance",
                        fluidRow(
                          column(width = 4, 
                                 radioButtons(inputId = "data_type",
                                              label = "Select one 'omic' layer:",
                                              choices = c("Phosphoproteomics", "Proteomics", "Transcriptomics"),
                                              selected = NULL),
                                 selectizeInput("gene", "Select one analyte:", choices = "- select", multiple = FALSE),
                                 radioButtons(
                                   inputId = "analyte_filter",
                                   label = tagList(
                                     "Select analyte filter:",
                                     icon("info-circle", id = "info_analyte_filter")
                                   ),
                                   choices = c("All Analytes" = "all", "Significant Analytes" = "significant"),
                                   selected = "all"
                                 ),
                                 
                                 bsPopover(id = "info_analyte_filter",
                                           title = "Analyte Filter",
                                           content = "Choose whether to display all analytes or only those with significant variation (Z-score < -1.96 or > +1.96).",
                                           placement = "right",
                                           trigger = "hover"),
                                 selectInput("comparison_type", "Select comparison type:", 
                                             choices = c("Tumor vs Tumor", "Tumor vs Normal"), 
                                             selected = "Tumor vs Tumor"),  
                                 
                                 conditionalPanel(
                                   condition = "input.comparison_type == 'Tumor vs Tumor'",
                                   selectizeInput(
                                     "tumors",
                                     "Select one or more tumors:",
                                     choices = c("Brca", "Ccrcc", "Coad", "Gbm", "Hnscc", "Lscc", "Luad", "Ov", "Pdac","Ucec", "all"),
                                     multiple = TRUE
                                   ),
                                   selectizeInput(
                                     "stages",
                                     "Select one, more or all stages:",
                                     choices = c(1:4, "all"),
                                     multiple = TRUE,
                                     options = list(id = "stage_selector")
                                   )
                                 ),
                                 conditionalPanel(
                                   condition = "input.comparison_type == 'Tumor vs Normal'",
                                   selectizeInput(
                                     "tumors_norm",
                                     "Select one tumor:",
                                     choices = c("Ccrcc", "Coad", "Gbm", "Hnscc", "Lscc", "Luad", "Ov", "Pdac","Ucec", "all"),
                                     multiple = FALSE
                                   ),
                                   selectizeInput(
                                     "stages",
                                     "Select one, more or all stages:",
                                     choices = c(1:4, "all"),
                                     multiple = TRUE,
                                     options = list(id = "stage_selector")
                                   )),
                                 
                                 checkboxInput("show_points", "Show Distribution", value = FALSE),
                                 checkboxInput("show_significance", "Show Significance (Wilcoxon test)", value = FALSE),
                                 
                                 div(
                                   actionButton("generate_plot", "Generate boxplot",class = "gene-selection-box" ),
                                   tags$div(
                                     style = "margin-top: 5px; font-size: 13px; color: #777;",
                                     icon("info-circle"), 
                                     "Remember to click this button every time to update the plot and the table."
                                   )
                                 )
                          ),
                          
                          column(width = 8,  
                                 box(
                                   width = 12, 
                                   style = "height: 500px;",
                                   withSpinner(
                                     plotOutput("combined_plot"),
                                     type = 4,
                                     color = "#3399FF",
                                     size = 0.5
                                   ),
                                   class = "plot-container",
                                   downloadButton(
                                     "download_pdf", 
                                     "Download plot (PDF)", 
                                     class = "btn btn-default btn-sm", 
                                     style = "margin-top: 10px;",
                                     onclick = "function() { Shiny.onInputChange('download_pdf', Date.now()); }"
                                   ),
                                   uiOutput("signor_info"),
                                   verbatimTextOutput("debug_link")
                                 )
                          ),
                          
                          column(width = 12,
                                 dataTableOutput("data_table"),
                                 selectInput("file_type", "Select file type for download:", 
                                             choices = c("xlsx", "csv", "tsv")),
                                 downloadButton("download_data", "Download data")
                          )
                        )
               ),
               
    
    #### Protein activity panel ####
    tabPanel("Explore Protein Activity",
             fluidRow(
               column(width = 4,
                      selectizeInput("protein", "Select one protein:", choices = NULL, multiple = FALSE),
                   
                      selectInput("comparison_type_protein", "Select comparison type:", 
                                  choices = c("Tumor vs Tumor", "Tumor vs Normal"), 
                                  selected = "Tumor vs Tumor"),  
                      
                      conditionalPanel(
                        condition = "input.comparison_type_protein == 'Tumor vs Tumor'",
                        selectizeInput(
                          "tumors_prot",
                          "Select one or more tumors:",
                          choices = c("Brca", "Ccrcc", "Coad", "Gbm", "Hnscc", "Lscc", "Luad", "Ov", "Pdac", "Ucec","all"),
                          multiple = TRUE
                        ),
                        selectizeInput(
                          "stages_prot",
                          "Select one, more or all stages:",
                          choices = c(1:4, "all"),
                          multiple = TRUE,
                          options = list(id = "stage_selector_prot")
                        )
                      ),
                      
                      conditionalPanel(
                        condition = "input.comparison_type_protein == 'Tumor vs Normal'",
                        selectizeInput(
                          "tumors_prot",
                          "Select one tumor:",
                          choices = c("Ccrcc", "Coad", "Gbm", "Hnscc", "Lscc", "Luad", "Ov", "Pdac", "Ucec", "all"),
                          multiple = FALSE
                        ),
                        selectizeInput(
                          "stages_prot",
                          "Select one, more or all stages:",
                          choices = c(1:4, "all"),
                          multiple = TRUE,
                          options = list(id = "stage_selector_prot")
                        )
                      ),
                      
                      checkboxInput("show_points_protein", "Show Distribution", value = FALSE),
                      checkboxInput("show_significance_protein", "Show Significance (Wilcoxon test)", value = FALSE),
                      
                      div(
                        actionButton("generate_protact_plot", "Generate boxplot", class = "gene-selection-box"),
                        tags$div(
                          style = "margin-top: 5px; font-size: 13px; color: #777;",
                          icon("info-circle"), 
                          "Remember to click this button every time to update the plot and the table."
                        )
                      )
                      
               ),
               
               column(width = 8,
                      box(
                        width = 12,
                        style = "height: 500px;",  
                        
                        withSpinner(
                          plotOutput("combined_plot_prot"),
                          type = 4,
                          color = "#3399FF",
                          size = 0.5
                        ),
                        class = "plot-container",
                        downloadButton(
                          "download_pdf_prot", 
                          "Download plot (PDF)", 
                          class = "btn btn-default btn-sm", 
                          style = "margin-top: 10px;",
                          onclick = "function() { Shiny.onInputChange('download_pdf_prot', Date.now()); }"
                        )
                      )
               ),
               
               column(width = 12,
                      dataTableOutput("data_table_prot"),
                      selectInput("file_type_prot", "Select file type for download:", 
                                  choices = c("xlsx", "csv", "tsv")),
                      downloadButton("download_data_prot", "Download data")
               )
             )))
    ),
    
    #### Extract protein activity from your data ####
    tabPanel("Extract Protein Activity from your data ",
             tabsetPanel(
               #### UI about SP ####
               tabPanel("About SignalingProfiler",
                        
                        h3(strong("SignalingProfiler Tool")),
                        
                        div(
                          style = "display: flex; justify-content: space-between; align-items: flex-start;",
                          
                          div(
                            style = "flex: 1;",
                            
                            p("This application uses the SignalingProfiler pipeline to perform protein activity inference as described in the following paper:"),
                            
                            tags$a(
                              href = "https://www.nature.com/articles/s41540-024-00417-6", 
                              "SignalingProfiler: A Tool for the Quantitative Analysis of Proteomics Data (Venafra V., Sacco F., Perfetto L., 2024)", 
                              target = "_blank"
                            ),
                            
                            hr(),
                            
                            tags$a(
                              href = "https://html-preview.github.io/?url=https://github.com/SaccoPerfettoLab/SignalingProfiler/blob/dev/vignettes/SignalingProfiler_tutorial_Bioc.html", 
                              "Click here to access the complete tutorial", 
                              target = "_blank"
                            ),
                            
                            hr(),
                            
                            p("The following three sample dataframes are provided to illustrate the proper format for providing omics 
                              dataframes to input to SignalingProfiler."),
                            
                            p("To ensure analysis, the columns must be renamed as they are in the examples."),
                            
                            p("All three omics dataframes can be loaded to initiate the analysis, or if one or two are missing, 
                              the analysis will be performed in the most optimal way based on the loaded dataframe(s)."),
                            
                            h4("Transcriptomics Data Example"),
                            tableOutput("example_transcriptomics"),
                            
                            hr(),
                            
                            h4("Proteomics Data Example"),
                            tableOutput("example_proteomics"),
                            
                            hr(),
                            
                            h4("Phosphoproteomics Data Example"),
                            tableOutput("example_phosphoproteomics")
                          ),
                          
                          div(
                            style = "flex: 0 0 auto; margin-left: 20px;",
                            img(src = "sp_logo.png", style = "max-width: 200px;")
                          )
                        )
               )
               ,
               
               #### UI run analysis SP ####
               tabPanel("Run Analysis",
                        fluidRow(
                          column(4,
                                 h4("Upload data for one or multiple samples"),
                                 
                                 p(style = "text-align: justify;",
                                   "Note: File names should follow the pattern ",
                                   strong("P_<ID>"), " for proteomics, ",
                                   strong("Ph_<ID>"), " for phosphoproteomics, and ",
                                   strong("T_<ID>"), " for transcriptomics. For example, if the sample ID is ",
                                   strong("C3L00010"), ", you should name the proteomics file: ",
                                   em("P_C3L00010.tsv or .csv or .xlsx")
                                 ),
                                 
                                 fileInput("upload_files", "Upload Data Files", multiple = TRUE),
                                 
                                 div(
                                   style = "background-color: #fff3cd; border-left: 5px solid #f0ad4e; padding: 15px; margin-top: 10px; font-size: 1em; border-radius: 8px; box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);",
                                   HTML(
                                     "⚠️ <strong>Warning:</strong> Please ensure your data contains columns in the order shown in the previous section.
                                    <br> In addition, it is also important to have in your data correct UniProt IDs and valid amino acids (S, T, or Y). 
                                    Otherwise, the data will undergo a <strong>partial cleanup</strong> and some entries might be removed.<br>
                                    For proper data preprocessing, refer to the <a href='https://github.com/SaccoPerfettoLab/PatientProfiler' target='_blank' style='color: #f39c12;'>PatientProfiler</a> package."
                                   )
                                 ),
                                 
                                 hr(),
                                
                                 #### organism selection sp ####
                                 div(
                                   style = "margin-bottom: 15px;",
                                   
                                   tags$label("Select Organism:"),
                                   
                                   radioButtons("organism", NULL,
                                                choices = c("👤 Human" = "human", "🐭 Mouse" = "mouse"),
                                                selected = "human",
                                                inline = TRUE),
                                   
                                   tags$div(
                                     style = "font-size: 0.9em; color: #666; margin-top: 8px; display: flex; gap: 5px;",
                                     icon("info-circle", id = "info_analyte_filter"),
                                     
                                     "Select Human or Mouse to apply organism-specific regulatory models during activity inference."
                                     
                                   )
                                 ),
                                 
                                 
                                 radioButtons("param_choice", "Choose Parameter Settings:",
                                              choices = c("Use Default Parameters", "Select Parameters"),
                                              selected = "Use Default Parameters"),
                                 
                                 div(
                                   style = "margin-top: 10px; display: flex; gap: 10px; align-items: center;",
                                   
                                   actionButton("run_analysis", "Run Analysis", class = "btn-primary"),
                                   
                                   tags$div(
                                     class = "tooltip-wrapper",
                                     actionButton("clear_data", "Clear Data", class = "btn-danger"),
                                     tags$div(class = "tooltip-text", 
                                              "This will clear all current data and results. You'll need to re-upload the files before running a new analysis.")
                                   )
                                 ),
                                 
                                 tags$p(
                                   style = "font-size: 1.1em; margin-top: 10px; display: flex; align-items: center; gap: 8px;",
                                   icon("info-circle"),
                                   "Press 'Clear Data' before starting a new analysis."
                                 ),
                                 
                                 conditionalPanel(
                                   condition = "input.param_choice == 'Select Parameters'",
                                   hr(),
                                   h4("TF Activity Footprint Analysis"),
                                   numericInput("reg_minsize_tf", "Regulon Minimum Size", value = 10, min = 1),
                                   checkboxInput("exp_sign_tf", "Use Expression Sign", value = FALSE),
                                   checkboxInput("hypergeom_corr_tf", "Apply Hypergeometric Correction", value = TRUE),
                                   checkboxInput("GO_annotation_tf", "Apply GO Annotation", value = TRUE),
                                   checkboxInput("correct_proteomics_tf", "Correct Proteomics", value = FALSE),
                                   
                                   hr(),
                                   h4("Kinase-Phosphatase Activity Footprint Analysis"),
                                   numericInput("reg_minsize_kin", "Regulon Minimum Size", value = 5, min = 1),
                                   checkboxInput("exp_sign_kin", "Use Expression Sign", value = FALSE),
                                   checkboxInput("integrated_regulons_kin", "Use Integrated Regulons", value = TRUE),
                                   checkboxInput("hypergeom_corr_kin", "Apply Hypergeometric Correction", value = TRUE),
                                   checkboxInput("GO_annotation_kin", "Apply GO Annotation", value = TRUE),
                                   checkboxInput("correct_proteomics_kin", "Correct Proteomics", value = FALSE),
                                   
                                   hr(),
                                   h4("PhosphoScore Computation"),
                                   checkboxInput("activatory", "Select Activatory Phosphosites Only", value = TRUE),
                                   checkboxInput("GO_annotation_phospho", "Apply GO Annotation", value = TRUE),
                                   hr()
                                 )
                          
                          
                          ),
                          column(8,
                                 box(
                                   width = 12,
                                   
                                   h3("Protein Activity Inference Results", style = "font-weight: bold; text-align: center; margin-bottom: 20px;"),
                                   
                                   withSpinner(
                                     DTOutput("results_table"),
                                     type = 4, color = "#3399FF", size = 0.5
                                   ),
                                   
                                   div(
                                     style = "text-align: center; margin-top: 15px;",
                                     downloadButton("download_results", "Download Results", class = "btn-primary")
                                   ),
                                   
                                   p(
                                     style = "font-weight: bold; color: #3ac554; margin-top: 10px; text-align: center;",
                                     "Once you have obtained the results, you can use them in the 'Plot Results' section to visualize the data."
                                   ),
                                   
                                   hr(style = "border-top: 2px solid #D3D3D3; margin-top: 20px;"),
                                   
                                   h4("Summary Table", style = "text-align: center; margin-bottom: 10px;"),
                                   
                                   withSpinner(
                                     DTOutput("summary_table"),
                                     type = 4, color = "#3399FF", size = 0.5
                                   ),
                                   
                                   div(
                                     style = "text-align: center; margin-top: 15px;",
                                     downloadButton("download_summary", "Download Summary", class = "btn-info")
                                   ),
                                   
                                   br(), br(),
                                   
                                   hr(style = "border-top: 2px solid #D3D3D3;")
                                 )
                          )
                          
                          
                          
                        )
               ),
               #### UI SP PLOT ####
               tabPanel("Plot Results",
                        h3("View the Top 15 Up and Down-Regulated Proteins Based on Their Predicted Activity"),
                        
                        p(style = "text-align: justify;", 
                          "Here you can upload your results from the protein activity inference to visualize your data 
                           in a barplot, separated by molecular function."
                        ),
                        
                        fileInput("results_file", "Upload Results File", accept = c(".csv", ".tsv", ".xlsx")),
                        
                        radioButtons("sample_selection", "Select Sample Data:", choices = c("All Samples", "Specific Sample")),
                        p(style = "color:blue;", "(Selecting 'All Samples' will show the mean predicted activity across your cohort)"),
                        
                        uiOutput("sample_selector_ui"),
                        
                        actionButton("generate_plot", "Generate Plot"),
                        downloadButton("download_plot", "Download Plot"),
                        
                        tags$br(),
                        
                        withSpinner(
                          plotOutput("protein_plot", height = "1000px"),
                          type = 4,
                          color = "#3399FF",
                          size = 0.5
                        ),
                        
                        tags$br(),
                        
                        tags$div(style = "font-size: 14px;", 
                                 "Molecular Functions (mf):",
                                 tags$ul(
                                   tags$li("tf: Transcription Factor"),
                                   tags$li("kin: Kinase"),
                                   tags$li("phos: Phosphatase"),
                                   tags$li("other: Other")
                                 )
                        )
               )
               )),
    
    #### Help ####
    tabPanel("Help",
             fluidRow(
               column(width = 6, offset = 3,
                      h2("📘 Help and Documentation"),
                      hr(),
                      h3("❓ Common Issues"),
                      tags$ul(
                        tags$li(
                          strong("Q:"), " I can't see the plots after selecting the analytes.",
                          br(),
                          strong("A:"), " Make sure you have selected at least one tumor and one stage. Also, ensure that the data for the selected analyte is available."
                        ),
                        tags$li(
                          strong("Q:"), " The data download button is not working.",
                          br(),
                          strong("A:"), " Please check your internet connection and ensure that your browser allows downloads."
                        )
                      ),
                      br(),
                      h3("💡 Tutorial"),
                      p("For a detailed tutorial on how to use MultiOmicsXplorer, please visit: ",
                        a("MultiOmicsXplorer Tutorial", 
                          href = "https://perfettolab.bio.uniroma1.it/MultiOmicsXplorer/tutorial.html", 
                          target = "_blank")
                      ),
                      br(),
                      h3("📨 Contact Information"),
                      p("For further assistance or inquiries, feel free to contact us:"),
                      tags$ul(
                        tags$li("PerfettoLab of Bioinformatics, University of Rome \"La Sapienza\""),
                        tags$li("Head of the Lab: Dr. Livia Perfetto"),
                        tags$li("Primary curator of MultiOmicsXplorer: Dr. Eleonora Meo"),
                        tags$li("Email: ", a("eleonorameo.hp@gmail.com", href = "mailto:eleonorameo.hp@gmail.com"))
                      ),
                      br(),
                      h3("🔬 About the Lab"),
                      p("PerfettoLab specializes in bioinformatics research and analysis, with a focus on cancer genomics and proteomics. 
                     We are dedicated to developing tools like MultiOmicsXplorer to facilitate data analysis for cancer research."),
                      p("Visit our website: ",
                        a("PerfettoLab Homepage", href = "https://sites.google.com/view/perfettolab/home", target = "_blank")
                      ),
                      br()
               )
             )
    )
  ),
  
  #### Footer ####
  tags$footer(
    div(
      style = "width: 100%; background-color: transparent; padding: 40px 20px 10px 20px; font-size: 14px; color: #939393; 
      margin-top: 100px !important; margin-bottom: 10px; display: flex; flex-direction: column; align-items: flex-start;",  
      img(src = "banner_sapienza.png", height = "170px", style = "margin-bottom: 10px;"),  
      
      div(
        style = "text-align: left;", 
        "For support, contact: ",
        a("eleonorameo.hp@gmail.com", href = "mailto:eleonorameo.hp@gmail.com", style = "color: #939393; font-weight: bold; 
          text-decoration: none;")
      )
    ),
    
    div(
      style = "width: 100%; background-color: #F5F5F5; padding: 5px; text-align: center; font-size: 12px; color: black;",
      "© 2025 MultiOmicsXplorer | All rights reserved | Developed by PerfettoLab, \"La Sapienza\" University of Rome."
    )
  )
  
  

  
  
)

##### SERVER #####

server <- function(input, output, session) {
  print(R.version)
  
  show_stats <- reactiveVal(FALSE)
  
  #### Load data abundance ####
  data_reactive <- reactive({
    req(input$generate_plot)
    
    data_type <- tolower(input$data_type)
    gene <- input$gene
    comparison_type <- input$comparison_type
    
    # upload data
    data <- load_data(data_type, gene)
    
    tumor_types <- input$tumors
    if ("all" %in% tumor_types) {
      tumor_types <- unique(data$Tumor_Type)
    }
    
    if (comparison_type == "Tumor vs Tumor") {
      stages <- input$stages
    } else if (comparison_type == "Tumor vs Normal") {
      stages <- NULL  
    } else {
      return(NULL)
    }
    
    return(data)
  })
  
  
  #### Generation plot abundance ####
  p <- eventReactive(input$generate_plot,{
    
    shinyjs::addClass(id = "generate_plot", class = "clicked")
    
    shinyjs::delay(1000, {
      shinyjs::removeClass(id = "generate_plot", class = "clicked")
    })
  
    if (is.null(input$gene) || input$gene == "" || input$gene == "- select") {
      return(
        ggplot() +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "Please select an analyte", size = 6, color = "red")
      )
    }
    
    data_type <- tolower(input$data_type)
    gene <- input$gene
    comparison_type <- input$comparison_type
    show_significance <- input$show_significance  
    
    if (is.null(gene) || is.null(comparison_type)) {
      return(NULL)
    }
    
    data <- data_reactive()
    
    req(data)
 
    if (comparison_type == "Tumor vs Tumor") {
      tumor_types <- input$tumors

      if ("all" %in% tumor_types) {
        tumor_types <- unique(data$Tumor_Type)
      }
      
      stages <- input$stages
      plot_obj <- generate_plot_abundance(data, 
                                          gene, 
                                          tumor_types, 
                                          stages, 
                                          data_type, 
                                          comparison_type, 
                                          analyte_option = input$analyte_filter,
                                          input$show_points, 
                                          show_significance)
    } else if (comparison_type == "Tumor vs Normal") {
      tumor_types <- input$tumors_norm[1]
      
      stages <- input$stages
      plot_obj <- generate_plot_abundance(data, 
                                          gene, 
                                          tumor_types, 
                                          stages, 
                                          data_type, 
                                          comparison_type, 
                                          analyte_option = input$analyte_filter,
                                          input$show_points,
                                          show_significance)
    }
    
    if (is.null(plot_obj$data) || nrow(plot_obj$data) == 0) {
      showNotification("No data available for the selected tumor types and stages.", type = "warning", duration = 5)
      return(
        ggplot() +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "No data available for this selection", size = 6, color = "orange")
      )
    } else { 
      req(plot_obj)
      return(plot_obj)  
    }
   
   
  })
  
  #### Output plot analyte abundance ####
  output$combined_plot <- renderPlot({
    plot_obj <- p()  
    if (!is.null(plot_obj)) {
      return(plot_obj)  
    } else {
      message("Plot not available.")  
      return(NULL)
    }
  })
  

  
  #### Generation of Signor link ####
  observeEvent(input$generate_plot, {
    plot_obj <- p()  
    req(plot_obj)  
    
    output$signor_info <- renderUI(NULL)  
    
    if (input$data_type == "Phosphoproteomics" && !is.null(input$gene)) {
      
      analyte <- input$gene
      if (grepl("_", analyte)) {
        parts <- unlist(strsplit(analyte, "_"))
        protein_name <- parts[1]
        aa_dict <- list("S" = "Ser", "T" = "Thr", "Y" = "Tyr")
        residue_code <- substr(parts[2], 1, 1)
        position <- substr(parts[2], 2, nchar(parts[2]))
        
        if (residue_code %in% names(aa_dict)) {
          residue <- paste0(aa_dict[[residue_code]], position)
          
          uniprot_data <- read_tsv("fosfosit_uniprot_table.tsv")
          selected_row <- uniprot_data %>% dplyr::filter(Name == analyte)
          
          if (nrow(selected_row) > 0) {
            uniprot_id <- selected_row$UNIPROT[1]
            
            network_link <- query_ptm_residue(uniprot_id, protein_name, residue)
            
            output$signor_info <- renderUI({
              if (!is.null(network_link)) {
                tags$a(href = network_link, 
                       target = "_blank", 
                       style = "color: #FD7F20; font-weight: bold;", 
                       "View in Signor for detailed upstream and downstream interactions at this phosphorylation site.", 
                       tags$img(src = "signor_logo.png", height = "30px", style = "vertical-align: middle; margin-left: 5px;"))
              } else {
                tags$p("No network link found for this phosphorylation site.", style = "color: orange;")
              }
            })
          } else {
            cat("Protein not found in UniProt table.\n")  # Debug log
            output$signor_info <- renderUI({
              tags$p("Protein not found in UniProt table.", style = "color: red;")
            })
          }
        }
      }
    }
  })
  
  
  observeEvent(input$generate_plot, {
    show_stats(FALSE)
  })
  
  observeEvent(input$generate_stats_plot, {
    show_stats(TRUE)
  })
  
  #### Analyte selection for data type (omic layer) ####
  observeEvent(input$data_type, {
    data_type_uc <- tools::toTitleCase(input$data_type)
    data_type_prefix <- switch(data_type_uc,
                               "Proteomics" = "Prot",
                               "Phosphoproteomics" = "Phospho",
                               "Transcriptomics" = "Transc")
    genes_file <- paste0(data_type_prefix, "_genes.tsv")
    
    if (fs::file_exists(genes_file)) {
      gene_list <- readr::read_tsv(genes_file)$gene_name
      updateSelectizeInput(session, "gene", choices = gene_list, selected = "- select", server = TRUE)
    } else {
      updateSelectizeInput(session, "gene", choices = NULL, selected = NULL, server = TRUE)
    }
  })
  
  #### Download plot abundance ####
  output$download_pdf <- downloadHandler(
    filename = function() {
      clean_name <- gsub("[^A-Za-z0-9_]", "_", input$gene)
      paste0("analyte_abundance_plot_", clean_name, "_", Sys.Date(), ".pdf")    
    },
    content = function(file) {
      pdf_file <- tempfile(fileext = ".pdf")  
      ggsave(pdf_file, plot = p(), device = "pdf", width = 12, height = 6)  
      file.rename(pdf_file, file) 
    },
    contentType = "application/pdf"  
  )
  
  
  #### Output table abundance ####
  table_data <- eventReactive(input$generate_plot, {
    req(input$data_type, input$gene, input$comparison_type)
    
    data_type <- tolower(input$data_type)
    gene <- input$gene
    comparison_type <- input$comparison_type
    analyte_option <- input$analyte_filter
    
    data <- load_data(data_type, gene)
    req(data)  
    
    if (comparison_type == "Tumor vs Tumor") {
      tumor_types <- input$tumors
      if ("all" %in% tumor_types) {
        tumor_types <- unique(data$Tumor_Type)
      }
      stages <- input$stages
      
      if (is.null(tumor_types) || length(tumor_types) == 0) {
        return(NULL)
      }
      
      filtered_data <- if ("all" %in% stages) {
        data %>% filter(Name == gene, Tumor_Type %in% tumor_types, !grepl("\\.N$", Patient_ID))
      } else {
        data %>% filter(Name == gene, Tumor_Type %in% tumor_types, Stage %in% stages, !grepl("\\.N$", Patient_ID))
      }
      
      if (analyte_option == "significant") {
        filtered_data <- filtered_data %>% filter(z_score < -1.96 | z_score > 1.96)
      }
      
      if (nrow(filtered_data) == 0) return(NULL)
      
      return(filtered_data)
      
    } else if (comparison_type == "Tumor vs Normal") {
      tumor_type <- input$tumors_norm  
      stages <- input$stages
      if (is.null(stages)) stages <- unique(data$Stage)  
      
      if (is.null(tumor_type)) return(NULL)
      
      filtered_data <-
        generate_table_norm_ab(data, gene, tumor_type, stages, analyte_option)
      
      
      if (is.null(filtered_data) || nrow(filtered_data) == 0) return(NULL)
      
      return(filtered_data)
    }
    
    return(NULL)
  })
  
  output$data_table <- renderDataTable({
    req(table_data())  
    datatable(table_data())  
  })
  
  #### Download table abundance ####
  output$download_data <- downloadHandler(
    filename = function() {
      clean_name <- gsub("[^A-Za-z0-9_]", "_", input$gene)
      paste0("analyte_abundance_table_", clean_name, "_", Sys.Date(), ".", input$file_type)
    },
    content = function(file) {
      data_to_download <- table_data()  
      
      if (nrow(data_to_download) == 0) {
        return(NULL)
      }
      
      if (input$file_type == "xlsx") {
        openxlsx::write.xlsx(data_to_download, file)
      } else if (input$file_type == "csv") {
        write.csv(data_to_download, file, row.names = FALSE)
      } else if (input$file_type == "tsv") {
        write.table(data_to_download, file, sep = "\t", row.names = FALSE)
      }
    }
  )
  
  
  #### Load data protein activity ####
  observe({
    protein_choices <- load_protein_choices()
    updateSelectizeInput(session, "protein", choices = protein_choices, server = TRUE)
  })
  
  data_reactive_prot <- reactive({
    req(input$generate_protact_plot)
    protein <- input$protein
    comparison_type <- input$comparison_type_protein
    
    data <- load_data_protein(protein)
    
    if (comparison_type == "Tumor vs Tumor") {
      tumor_types <- input$tumors_prot
      if ("all" %in% tumor_types) {
        tumor_types <- unique(data$Tumor_Type)
      }
      stages <- input$stages_prot
    } else if (comparison_type == "Tumor vs Normal") {
      tumor_types <- input$tumors_prot
      stages <- NULL
    } else {
      return(NULL)
    }
    
    return(data)

  })
  
  #### Generation plot protein activity ####
  p_prot <- eventReactive(input$generate_protact_plot,{
    
    shinyjs::addClass(id = "generate_protact_plot", class = "clicked")
    
    shinyjs::delay(1000, {
      shinyjs::removeClass(id = "generate_protact_plot", class = "clicked")
    })
    
    if (is.null(input$protein) || input$protein == "") {
      return(
        ggplot() +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "Please select a protein", size = 6, color = "red")
      )
    }
    
    protein <- input$protein
    comparison_type_prot <- input$comparison_type_protein
    show_points_protein <- input$show_points_protein
    
    if (is.null(protein) || is.null(comparison_type_prot)) {
      return(NULL)
    }
    
    data <- data_reactive_prot()
    
    if (input$generate_protact_plot > 0) {
      if (comparison_type_prot == "Tumor vs Tumor") {
        req(input$protein, input$tumors_prot, input$stages_prot)
        tumors <- input$tumors_prot
        if ("all" %in% tumors) {
          tumors <- unique(data$Tumor_Type)
        }
        stages <- input$stages_prot
        plot_obj <- generate_plot_protein_activity(
          data, 
          protein, 
          tumors, 
          stages, 
          comparison_type_prot, 
          show_points_protein = show_points_protein, 
          show_significance = input$show_significance_protein
        )
      } else if (comparison_type_prot == "Tumor vs Normal") {
        req(input$protein, input$tumors_prot, input$stages_prot)
        tumors <- input$tumors_prot[1]
        stages <- input$stages_prot
        plot_obj <- generate_plot_protein_activity(
          data, 
          protein, 
          tumors, 
          stages, 
          comparison_type_prot, 
          show_points_protein = show_points_protein, 
          show_significance = input$show_significance_protein
        )
      }
    } 
    
    if (is.null(plot_obj$data) || nrow(plot_obj$data) == 0) {
      showNotification("No data available for the selected tumor types and stages.", type = "warning", duration = 5)
      return(
        ggplot() +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "No data available for this selection", size = 6, color = "orange")
      )
    } else { 
      req(plot_obj)
      return(plot_obj)  
    }
    
  })
  
  
  output$combined_plot_prot <- renderPlot({
    p_prot()
  })
  
  
  #### Download plot activity ####
  output$download_pdf_prot <- downloadHandler(
    filename = function() {
      clean_name <- gsub("[^A-Za-z0-9_]", "_", input$protein)
      paste0("protein_activity_plot_", clean_name, "_", Sys.Date(), ".pdf")    
      },
    content = function(file) {
      pdf_file <- tempfile(fileext = ".pdf")
      ggsave(pdf_file, plot = p_prot(), device = "pdf", width = 12, height = 6)
      file.rename(pdf_file, file)  
    },
    contentType = "application/pdf"  
  )
  
  
  #### Output table protein activity ####
  table_data_prot <- eventReactive(input$generate_protact_plot, {
    req(input$protein, input$comparison_type_protein)
    
    protein <- input$protein
    comparison_type <- input$comparison_type_protein
    analyte_option <- input$analyte_filter_prot
    
    data <- data_reactive_prot()
    req(data)  
    
    if (comparison_type == "Tumor vs Tumor") {
      tumor_types <- input$tumors_prot
      if ("all" %in% tumor_types) {
        tumor_types <- unique(data$Tumor_Type)
      }
      stages <- input$stages_prot
      
      if (is.null(tumor_types) || length(tumor_types) == 0) {
        return(NULL)
      }
      
      filtered_data <- if ("all" %in% stages) {
        data %>% filter(Name == protein, Tumor_Type %in% tumor_types, !grepl("\\.N$", Patient_ID))
      } else {
        data %>% filter(Name == protein, Tumor_Type %in% tumor_types, Stage %in% stages, !grepl("\\.N$", Patient_ID))
      }
      
      if (nrow(filtered_data) == 0) return(NULL)
      
      return(filtered_data)
      
    } else if (comparison_type == "Tumor vs Normal") {
      tumor_type <- input$tumors_prot[1]
      stages <- input$stages_prot
      if (is.null(stages)) stages <- unique(data$Stage)
      
      if (is.null(tumor_type)) return(NULL)
      
      filtered_data <- generate_table_norm_ac(data, protein, tumor_type, stages)
      
      if (is.null(filtered_data) || nrow(filtered_data) == 0) return(NULL)
      
      return(filtered_data)
    }
    
    return(NULL)
  })
  
  output$data_table_prot <- renderDataTable({
    req(table_data_prot())  
    datatable(table_data_prot())  
  })
  
  #### Download table protein activity ####
  output$download_data_prot <- downloadHandler(
    filename = function() {
      clean_name <- gsub("[^A-Za-z0-9_]", "_", input$protein)
      paste0("protein_activity_table_", clean_name, "_", Sys.Date(), ".", input$file_type_prot)
    },
    
    content = function(file) {
      data_to_download <- table_data_prot()  
      
      if (nrow(data_to_download) == 0) {
        return(NULL)
      }
      
      if (input$file_type_prot == "xlsx") {
        openxlsx::write.xlsx(data_to_download, file)
      } else if (input$file_type_prot == "csv") {
        write.csv(data_to_download, file, row.names = FALSE)
      } else if (input$file_type_prot == "tsv") {
        write.table(data_to_download, file, sep = "\t", row.names = FALSE)
      }
    }
  )
  
  
  
  #### SIGNALING PROFILER ####
  #### Example data for about tab ####
  output$example_transcriptomics <- renderTable({
    data.frame(
      gene_name = c("A4GALT", "AAAS", "AACS", "AAED1", "AAGAB", "AAMDC", "AAMP", "AAR2", "AARS", "AARS2"),
      difference = c(-0.496, -0.850, -0.691, 1.33, 0.594, 0.581, 0.124, -0.646, 0.0457, -1.56),
      logpval = c(NA, 2.39, NA, NA, NA, NA, NA, 1.94, NA, NA),
      significant = c(NA, "+", NA, NA, NA, NA, NA, NA, NA, NA)
    )
  })
  
  output$example_proteomics <- renderTable({
    data.frame(
      gene_name = c("AAAS", "AACS", "AAGAB", "AAK1", "AAMDC", "AAMP", "AAR2", "AARS1", "AARS2", "AARSD1"),
      UNIPROT = c("Q9NRG9;F8VZ44", "Q86V21;Q86V21-2;Q86V21-3", "Q6PD74", "Q2M2I8-2;Q2M2I8;E9PG46", "Q9H7C9;E9PNP3;E9PR47;E9PIQ4", 
                  "C9JG97;Q13685;C9JEH3;C9JTS3", "Q9Y312;A2A2Q9", "P49588", "Q5JTZ9", "Q9BTE6;Q9BTE6-2"),
      difference = c(-4.38, 0.239, 0.246, 0.578, 0.178, 0.288, 2.64, 0.254, -3.19, -0.988),
      logpval = c(1.99, NA, NA, NA, NA, NA, 2.04, NA, NA, NA),
      significant = c("+", NA, NA, NA, NA, NA, "+", NA, NA, NA)
    )
  })
  
  output$example_phosphoproteomics <- renderTable({
    data <- data.frame(
      UNIPROT = c("Q5T3U5", "Q5T3U5", "O14639", "O14639", "O14639", "P60709", "Q5BKX5", "Q5BKX5", "O75078", "M0QY12"),
      gene_name = c("ABCC10", "ABCC10", "ABLIM1", "ABLIM1", "ABLIM1", "ACTB", "ACTMAP", "ACTMAP", "ADAM11", "ADAMTS10"),
      aminoacid = c("S", "S", "S", "S", "S", "S", "S", "S", "S", "T"),
      position = c(852, 854, 431, 435, 452, 239, 6, 14, 687, 4),
      sequence_window = c("EGLEEEQSTSGRLLQ", "LEEEQSTSGRLLQEE", "GESPRTLSPTPSAEG", "RTLSPTPSAEGYQDV", "RMIHRSTSQGSINSP", 
                          "SSSSLEKSYELPDGQ", "__MTSPCSPPLKPPI", "PPLKPPISPPKTPVP", "SGERRICSHHGVCSN", "____MGPTSVLRAGL"),
      difference = c(1.85, 2.14, -0.308, -0.308, 2.42, 0.808, -0.605, -0.605, 0.915, 1.94),
      logpval = c(1.44, 1.40, NA, NA, NA, 2.75, NA, NA, NA, NA),
      significant = c("+", "+", NA, NA, NA, "+", NA, NA, NA, NA)
    )
    
    data$position <- format(data$position, big.mark = "") 
    data$difference <- format(data$difference, nsmall = 2)  
    
    return(data)
  })
  
  

  #### Upload data for SP analysis ####
  sample_data_list <- reactiveVal(list())
  
  # Create session-specific temp file path -> this is important to not overlap analyses!!!
  user_temp_dir <- tempdir()
  file_path <- file.path(user_temp_dir, paste0("user_data_", as.integer(Sys.time()), "_", sample(1000:9999, 1), ".rds"))

  #### Reactive values for results and summaries ####
  all_results <- reactiveVal(data.frame(
    UNIPROT = character(),
    gene_name = character(),
    mf = character(),      
    predicted_activity = numeric(),
    method = character(),
    Sample_ID = character()
    
  ))
  
  all_summaries <- reactiveVal(data.frame(
    Sample_ID = character(), 
    tf_count = numeric(), 
    kin_count = numeric(), 
    phos_count = numeric(), 
    other_count = numeric()
  ))
  
  observeEvent(input$upload_files, {
    req(input$upload_files)
    
    files <- input$upload_files
    file_paths <- files$datapath
    file_names <- files$name
    
    current_data <- list()
    
    for (i in seq_along(file_paths)) {
      filepath <- file_paths[i]
      filename <- file_names[i]
      
      sample_id <- sub("^[^_]+_(.+?)\\..+$", "\\1", filename)
      
      if (grepl("^T_", filename)) {
        omic_type <- "Transcriptomics"
      } else if (grepl("^P_", filename)) {
        omic_type <- "Proteomics"
      } else if (grepl("^Ph_", filename)) {
        omic_type <- "Phosphoproteomics"
      } else {
        shinyjs::alert("File type not recognized. Please upload files with correct naming conventions.")
        stop("File type not recognized")
      }
      
      data <- read_data(filepath)
      
      check_result <- check_columns(data, omic_type)
      if (check_result != TRUE) {
        shinyjs::alert(paste("Error in file:", filename, "\n", check_result))
        return(NULL)
      }
      
      if (is.null(current_data[[sample_id]])) {
        current_data[[sample_id]] <- list(
          Transcriptomics = NULL,
          Proteomics = NULL,
          Phosphoproteomics = NULL
        )
      }
      
      current_data[[sample_id]][[omic_type]] <- data
    }
    
    saveRDS(current_data, file_path)
    sample_data_list(current_data)
  })
  
  
  #### Run SP analysis ####
  
  observeEvent(input$run_analysis, {
    shinyjs::disable("run_analysis")
    notification_id <- showNotification("Running analysis, please wait...", type = "message", duration = NULL)
    
    params <- list(
      
      organism = input$organism,
      
      reg_minsize_tf = ifelse(input$param_choice == "Use Default Parameters", 10, input$reg_minsize_tf),
      exp_sign_tf = ifelse(input$param_choice == "Use Default Parameters", FALSE, input$exp_sign_tf),
      hypergeom_corr_tf = ifelse(input$param_choice == "Use Default Parameters", TRUE, input$hypergeom_corr_tf),
      GO_annotation_tf = ifelse(input$param_choice == "Use Default Parameters", TRUE, input$GO_annotation_tf),
      correct_proteomics_tf = ifelse(input$param_choice == "Use Default Parameters", FALSE, input$correct_proteomics_tf),
      
      reg_minsize_kin = ifelse(input$param_choice == "Use Default Parameters", 5, input$reg_minsize_kin),
      exp_sign_kin = ifelse(input$param_choice == "Use Default Parameters", FALSE, input$exp_sign_kin),
      integrated_regulons_kin = ifelse(input$param_choice == "Use Default Parameters", TRUE, input$integrated_regulons_kin),
      hypergeom_corr_kin = ifelse(input$param_choice == "Use Default Parameters", TRUE, input$hypergeom_corr_kin),
      GO_annotation_kin = ifelse(input$param_choice == "Use Default Parameters", TRUE, input$GO_annotation_kin),
      correct_proteomics_kin = ifelse(input$param_choice == "Use Default Parameters", FALSE, input$correct_proteomics_kin),
      
      activatory = ifelse(input$param_choice == "Use Default Parameters", TRUE, input$activatory),
      GO_annotation_phospho = ifelse(input$param_choice == "Use Default Parameters", TRUE, input$GO_annotation_phospho)
      
    )
    
    all_results_data <- data.frame()
    temp_summaries <- data.frame()
    
    sample_data <- sample_data_list()
    
    if (is.null(sample_data) || length(sample_data) == 0 || all(sapply(sample_data, is.null))) {
      showNotification("No sample data available.", type = "error")
      shinyjs::enable("run_analysis")
      removeNotification(notification_id)
      return(NULL)
    }
    
    for (sample_id in names(sample_data)) {
      trans_data <- sample_data[[sample_id]][["Transcriptomics"]]
      prot_data <- sample_data[[sample_id]][["Proteomics"]]
      phospho_data <- sample_data[[sample_id]][["Phosphoproteomics"]]
      
      sample_params <- params
      
      if (is.null(prot_data)) {
        if (params$correct_proteomics_tf || params$correct_proteomics_kin) {
          showNotification(
            paste("Sample", sample_id, ": proteomics data not provided. Correction of proteomics will be skipped for this sample."),
            type = "warning", duration = 5
          )
          sample_params$correct_proteomics_tf <- FALSE
          sample_params$correct_proteomics_kin <- FALSE
        }
      }
      
      analysis_output <- perform_analysis(trans_data, prot_data, phospho_data, sample_params, sample_id)
      
      if (!is.null(analysis_output$combined_results)) {
        all_results_data <- bind_rows(all_results_data, analysis_output$combined_results)
        temp_summaries <- bind_rows(temp_summaries, analysis_output$molecular_function_summary)
      }
    }
    
    all_results(all_results_data)
    all_summaries(temp_summaries)
    
    shinyjs::enable("run_analysis")
    removeNotification(notification_id)
  })
  
  #### Remove temp file ####
  session$onSessionEnded(function() {
    if (file.exists(file_path)) {
      file.remove(file_path)
    }
  })
  
  #### output SP and summary tables ####
  output$results_table <- renderDataTable({
  
    removeNotification("running_analysis")  
    all_results()
  })
  
  output$summary_table <- renderDataTable({
  
    all_summaries()
  })
  
  # Download handler for complete results
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("signaling_profiler_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(all_results(), file, row.names = FALSE)
    },
    contentType = "text/csv"  
  )
  
  # Download handler for summary table
  output$download_summary <- downloadHandler(
    filename = function() {
      paste0("summary_table_mf_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(all_summaries(), file, row.names = FALSE)
    },
    contentType = "text/csv"  
  )
  
  #### SP clear data ####
  observeEvent(input$clear_data, {
    sample_data_list(list())
    
    if (file.exists(file_path)) {
      file.remove(file_path)
    }
    
    all_results <- reactiveVal(data.frame(
      UNIPROT = character(),
      gene_name = character(),
      mf = character(),      
      predicted_activity = numeric(),
      method = character(),
      Sample_ID = character()
      
    ))
    
    all_summaries(data.frame(Sample_ID = character(), tf_count = numeric(), 
                             kin_count = numeric(), phos_count = numeric(), 
                             other_count = numeric()))
    
    showNotification("All data has been cleared.", type = "message")
    
    output$results_table <- renderDataTable({ all_results() })
    output$summary_table <- renderDataTable({ all_summaries() })
  })
  
  
  
  #### SP PLOT ####
  
  # Reactive variable for memorize results data
  results_data <- reactiveVal(NULL)
  
  # Upload data from results of predicted activity
  observeEvent(input$results_file, {
    req(input$results_file)
    
    # Read results data file
    file_path <- input$results_file$datapath
    results_df <- read.csv(file_path, stringsAsFactors = FALSE)  
    
    # Memorize results data in the reactive variable
    results_data(results_df)
  })
  
  # Dynamic UI for ID sample input
  output$sample_selector_ui <- renderUI({
    if (input$sample_selection == "Specific Sample") {
      sample_ids <- unique(results_data()$Sample_ID)
      selectInput("sample_id", "Select Sample ID", choices = sample_ids)
    } else {
      NULL
    }
  })
  
  
  #### Generate sp plot ####
  observeEvent(input$generate_plot, {
    req(results_data())  
    
    sample_id <- if (input$sample_selection == "Specific Sample") {
      req(input$sample_id)  
      input$sample_id
    } else {
      NULL
    }
    
    cat("Sample ID: ", sample_id, "\n")
    
    plot <- create_top_proteins_plot(results_data(), sample_id)
    
    output$protein_plot <- renderPlot({
      print(plot)  
    })
    
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("protein_activity_plot", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        plot_width <- session$clientData$output_protein_plot_width
        plot_height <- session$clientData$output_protein_plot_height
        
        ggsave(file, plot = plot, 
               width = plot_width / 72,  
               height = plot_height / 72, 
               units = "in", dpi = 300, bg = "white")
      }
    )
    
    
  })
  
  
}


# Combine UI and server in a Shiny obj
shinyApp(ui = ui, server = server)

