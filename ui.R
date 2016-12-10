library(shiny)
library(ggplot2)

ggplot2::theme_set(theme_classic() +
    theme(axis.line.x = element_blank()) + 
    theme(axis.line.y = element_blank()))

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
  fluidRow(column(12, h3("Data"))),
  fluidRow(
    column(4,
      wellPanel(
        textInput("var1", textOutput("var1_lbl"), value = "FGFR4"),
        textInput("var2", textOutput("var2_lbl"), value = "ERBB2"),
        actionButton("retreive_button", "Retrive TCGA data from cBioPortal"),
        checkboxInput("flip_axes", "Flip axes", value = FALSE)
      ))), 
  fluidRow(column(12, h3("Graphs"))),
  fluidRow(
    column(4,
      p("Mutations"),
      plotOutput("plot1")),
    column(4,
      p("Putative copy-number alterations (CNA)"),
      plotOutput("plot2")),
    column(4,
      p("mRNA expression"),
      plotOutput("plot3"), 
      p(textOutput("text3")),
      selectInput("smooth_method3", "Smoother",
        choices = c("(none)", "Linear regression", "Local polynomial regression (loess)"), 
        selected = "Local polynomial regression (loess)"))
  ),
  p('To save a figure to file, left-click/ctrl-click on the image and "Save Image As..." (or similar, depending on web browser).'),
  hr(),
  p(
    "© 2016 John Lövrot",
    br(),
    "This work is licensed under a",
    a("Creative Commons Attribution 4.0 International License",
      href = "http://creativecommons.org/licenses/by/4.0/",
      target = "_blank"),
    br(),
    "The source code is available at",
    a("github.com/lovrot/tcga-brca-explorer",
      href = "https://github.com/lovrot/tcga-brca-explorer",
      target = "_blank"),
    br(),
    "Version 0.0.0.9001")
)
