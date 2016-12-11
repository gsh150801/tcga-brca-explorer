library(shiny)

fluidPage(
  titlePanel("TCGA Breast Cancer Data Explorer"),
  fluidRow(
    column(12,
      p(
        "Basic exploration of TCGA breast cancer data retrieved from",
        a("cBioPortal",
          href = "http://www.cbioportal.org/",
          target = "_blank")),
      p("Please adhere to the", 
        a("TCGA publication guidelines", 
          href = "http://cancergenome.nih.gov/abouttcga/policies/publicationguidelines",
          target = "_blank"), 
        "when using TCGA data in your publications.", 
        br(),
        "Please cite",
        a("Gao et al. Sci. Signal. 2013", 
          href = "http://www.ncbi.nlm.nih.gov/pubmed/23550210",
          target = "_blank"), 
        "and", 
        a("Cerami et al. Cancer Discov. 2012", 
          href = "http://cancerdiscovery.aacrjournals.org/content/2/5/401.abstract",
          target = "_blank"), 
        "when publishing results based on cBioPortal."))
  ),  
  h3("Preparations"),
  fluidRow(
    column(4,
      wellPanel(
        textInput("query_str", "Gene set", value = "CDKN2A RB1 TP53"),
        actionButton("retrieve_button", "Retrieve TCGA data from cBioPortal")), 
      p(textOutput("retrieved_genes")) 
    ),
    column(4,
      wellPanel(
        uiOutput("var_y_ui"),
        uiOutput("var_x_ui")
      ))), 
  h3("Graphs"),
  p('To save a figure to file, left-click/ctrl-click on the image and "Save Image As..." (or similar, depending on web browser).'),
  checkboxInput("mark_mut", "Mark somatic point mutations in all graphs", value = FALSE),
  fluidRow(
    column(4,
      p("Somatic point mutations"),
      plotOutput("plot1"), 
      checkboxInput("show_mut", "Show somatic point mutations", value = FALSE),
      tableOutput("tab1")),
    column(4,
      p("Putative copy-number alterations (CNA)"),
      plotOutput("plot2")),
    column(4,
      p("mRNA expression"),
      plotOutput("plot3"), 
      p(textOutput("text3")),
      selectInput("smooth_method3", "Smoother",
        choices = c(
          "(none)" = "(none)", 
          "Linear regression" = "lm", 
          "Local polynomial regression (loess)" = "loess"), 
        selected = "loess"))
  ),
  hr(),
  p(
    "© 2016 John Lövrot",
    br(),
    "This data retrieval and visualisation tool is licensed under a",
    a("Creative Commons Attribution 4.0 International License",
      href = "http://creativecommons.org/licenses/by/4.0/",
      target = "_blank"),
    br(),
    "The source code is available at",
    a("github.com/lovrot/tcga-brca-explorer",
      href = "https://github.com/lovrot/tcga-brca-explorer",
      target = "_blank"),
    br(),
    "Version 0.0.0.9008")
)
