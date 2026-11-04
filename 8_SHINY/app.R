# =============================================================================
# Selection Bias in GWAS: Method Guide
# Companion to "Navigating Selection Bias in Genome-Wide Association Studies"
# Requires shiny_brain.csv (or shiny_brain.xlsx) in the same folder as app.R
# =============================================================================

library(shiny)
library(shinythemes)
library(dplyr)

# -----------------------------------------------------------------------------
# Load the method database ("brain")
# -----------------------------------------------------------------------------
read_brain = function(path_csv = "shiny_brain.csv", path_xlsx = "shiny_brain.xlsx") {
  
  if (file.exists(path_csv)) {
    lines = readLines(path_csv, warn = FALSE, encoding = "UTF-8")
    lines = lines[nzchar(trimws(lines))]
    lines[1] = sub("^\ufeff", "", lines[1])                 # strip byte-order mark
    
    # Repair rows saved by Excel as a single quoted cell
    if (all(grepl('^".*"$', lines))) {
      lines = sub('^"(.*)"$', "\\1", lines)
      lines = gsub('""', '"', lines)
    }
    
    sep = if (grepl(";", lines[1]) && !grepl(",", lines[1])) ";" else ","
    df = read.csv(text = lines, sep = sep, stringsAsFactors = FALSE,
                   na.strings = c("NA", ""), check.names = FALSE)
    
  } else if (file.exists(path_xlsx) && requireNamespace("readxl", quietly = TRUE)) {
    df = as.data.frame(readxl::read_excel(path_xlsx, sheet = 1, na = c("NA", "")),
                        stringsAsFactors = FALSE)
  } else {
    stop("Cannot find shiny_brain.csv (or shiny_brain.xlsx with the readxl package installed).")
  }
  
  df[] = lapply(df, function(x) if (is.character(x)) trimws(x) else x)
  
  needed = c("Method", "Purpose", "Data_level", "Selection_form", "Required_inputs",
              "Requirements", "Key_assumption", "Not_appropriate_when",
              "Paper", "Repo_script", "Software_link", "Software_status")
  missing = setdiff(needed, names(df))
  if (length(missing) > 0)
    stop("shiny_brain is missing columns: ", paste(missing, collapse = ", "))
  
  df
}

methods_db = read_brain()

# -----------------------------------------------------------------------------
# Lookup values (must match codes and labels used in shiny_brain)
# -----------------------------------------------------------------------------
input_labels = c(
  "Participation indicator (participants vs non-participants)" = "participation_indicator",
  "Population reference data (e.g. census or register)"        = "population_reference",
  "Covariate data on non-participants (sampling frame)"        = "nonparticipant_data",
  "A variable that predicts participation but not the trait"   = "exclusion_variable",
  "GWAS of participation (summary statistics)"                 = "participation_gwas",
  "Sex-stratified GWAS"         = "sex_stratified",
  "Disease incidence GWAS (for progression analyses)"          = "incidence_gwas",
  "Familial (e.g. parental) outcome data"                      = "familial_outcomes",
  "Population disease prevalence"                              = "disease_prevalence",
  "A negative control phenotype"                               = "negative_control"
)
label_of = setNames(names(input_labels), input_labels)

selection_forms = c("Non-representative participation",
                     "Sex-differential participation",
                     "Case-control secondary phenotype",
                     "Index-event (progression among cases)",
                     "Survival to recruitment")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
make_links = function(x, label = "Link") {
  if (is.na(x) || trimws(x) == "" || toupper(trimws(x)) == "NA") return(HTML("—"))
  links = unlist(strsplit(trimws(x), "\\s+"))
  links = links[grepl("^https?://", links)]
  if (length(links) == 0) return(HTML("—"))
  HTML(paste0('<a href="', links, '" target="_blank">', label, "</a>", collapse = " | "))
}

# Required_inputs: ";" = needs all of these, "|" = needs any one of these
missing_inputs = function(req, have) {
  if (is.na(req) || req == "") return(character(0))
  groups = strsplit(req, ";")[[1]]
  miss = character(0)
  for (g in groups) {
    opts = strsplit(trimws(g), "\\|")[[1]]
    if (!any(opts %in% have))
      miss = c(miss, paste(unname(label_of[opts]), collapse = " or "))
  }
  miss
}

show_error = function(e) {
  tags$p(paste("Error:", conditionMessage(e)), style = "color:red;font-weight:bold")
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui = fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Selection Bias in GWAS: Method Guide"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Step 1: Purpose"),
      radioButtons("purpose", "Do you want to detect or correct selection bias?",
                   choices = list("Diagnose" = "Diagnostic", "Correct" = "Corrective")),
      hr(),
      
      h4("Step 2: Data"),
      radioButtons("data_level", "What GWAS data do you have?",
                   choices = list("Individual-level genotype and phenotype" = "Individual",
                                  "Summary statistics only" = "Summary")),
      hr(),
      
      h4("Step 3: Form of selection"),
      selectInput("selection_form", "Which form of selection do you suspect?",
                  choices = as.list(c("Not sure", selection_forms))),
      hr(),
      
      h4("Available inputs (optional)"),
      checkboxGroupInput("inputs", "Tick what you have:", choices = as.list(input_labels)),
      
      helpText("This tool lists methods applicable to your setting. It does not rank them;",
               "see the paper and repository for assumptions and worked examples.")
    ),
    
    mainPanel(
      h3("Applicable methods"),
      uiOutput("notes"),
      uiOutput("method_table"),
      hr(),
      tags$p(HTML(paste0(
        "Companion to <em>Navigating Selection Bias in Genome-Wide Association Studies</em>. ",
        'Code and worked examples: <a href="[GitHub URL]" target="_blank">GitHub repository</a>.'
      ))),
      tags$p(HTML("Created by <strong>Ciarrah Jane Shannon Barry</strong>."),
             style = "font-size: 0.9em; color: #555;")
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server = function(input, output) {
  
  filtered = reactive({
    df = methods_db %>% filter(Purpose == input$purpose)
    
    # Individual-level users can also generate summary statistics, so see all methods
    if (input$data_level == "Summary")
      df = df %>% filter(Data_level %in% c("Summary", "Either"))
    
    if (input$selection_form != "Not sure")
      df = df %>% filter(Selection_form %in% c(input$selection_form, "Any"))
    
    df
  })
  
  output$notes = renderUI({
    tryCatch({
      df = filtered()
      notes = list()
      
      if (input$purpose == "Corrective" && any(grepl("Index-event", df$Selection_form)))
        notes = c(notes, list(tags$p(
          "DHO, Slope-Hunter and the MR-Horse adaptation identify the bias under different ",
          "untestable assumptions. Apply all three and compare.")))
      
      if (input$purpose == "Corrective" && input$data_level == "Summary" &&
          input$selection_form %in% c("Not sure", "Non-representative participation"))
        notes = c(notes, list(tags$p(
          "Summary-level options for non-representative participation are limited and ",
          "require a participation GWAS or sex-stratified results.")))
      
      if (input$purpose == "Diagnostic")
        notes = c(notes, list(tags$p(
          "No single diagnostic is conclusive; use several in combination.")))
      
      if (length(notes) > 0) div(class = "alert alert-info", notes) else NULL
    }, error = show_error)
  })
  
  output$method_table = renderUI({
    tryCatch({
      df = filtered()
      
      if (nrow(df) == 0)
        return(tags$p("No methods match these choices.", style = "color:red;font-weight:bold"))
      
      rows = lapply(seq_len(nrow(df)), function(i) {
        miss = missing_inputs(df$Required_inputs[i], input$inputs)
        all_req = missing_inputs(df$Required_inputs[i], character(0))
        
        need = if (length(input$inputs) == 0) {
          if (length(all_req) == 0) HTML("—")
          else tags$span(paste(all_req, collapse = "; "), style = "color:#555")
        } else if (length(miss) == 0) {
          tags$span("All listed inputs available", style = "color:green")
        } else {
          tags$span(paste(miss, collapse = "; "), style = "color:#b35900")
        }
        
        tags$tr(
          tags$td(tags$strong(df$Method[i])),
          tags$td(df$Selection_form[i]),
          tags$td(df$Requirements[i]),
          tags$td(df$Key_assumption[i]),
          tags$td(df$Not_appropriate_when[i]),
          tags$td(need),
          tags$td(make_links(df$Paper[i], "Paper")),
          tags$td(make_links(df$Repo_script[i], "Script")),
          tags$td(make_links(df$Software_link[i], "Software"))
        )
      })
      
      headers = c("Method", "Selection addressed", "Requires", "Key assumption",
                   "Not appropriate when", "You still need", "Paper", "Repository", "Software")
      
      tags$table(
        class = "table table-striped table-hover",
        tags$thead(tags$tr(lapply(headers, tags$th))),
        tags$tbody(rows)
      )
    }, error = show_error)
  })
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
shinyApp(ui = ui, server = server)

## To redeploy: unlink("rsconnect", recursive = TRUE); rsconnect::deployApp()
