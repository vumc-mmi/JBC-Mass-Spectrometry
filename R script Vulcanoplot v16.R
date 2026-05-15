# Set 30 MB upload limit
options(shiny.maxRequestSize = 30 * 1024^2) 

#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load Libraries
packages_to_install <- c("shiny", "ggplot2", "dplyr", "readr", "limma", "ggrepel", "DT", "colourpicker")
new_packages <- packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(limma)
library(ggrepel)
library(DT)
library(colourpicker) # For color inputs

# Define UI
ui <- fluidPage(
  
  # Application title
  titlePanel("Volcano Plot Generator"),
  
  # Sidebar with user inputs
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$h3("Instructions"),
      
      # File Upload
      tags$h4("1. Upload Data"),
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Group Selection
      tags$h4("2. Define Groups"),
      p("Select the 'LFQ intensity' columns for each group."),
      uiOutput("group1_selector"),
      uiOutput("group2_selector"),
      
      # Processing Options
      tags$h4("3. Processing Options"),
      checkboxInput("median_normalize", "Median Normalize LFQ Intensities", value = TRUE),
      
      # Tryptic Site Filtering
      tags$hr(),
      tags$h5("Tryptic Site Filtering"),
      p("Select which peptide types to include in the analysis. 
        (Requires 'Amino acid before' and 'Last amino acid' columns)."),
      checkboxGroupInput("tryptic_filter", "Include peptides that are:",
                         choices = c("Both N- and C-term tryptic (K/R)" = "both",
                                     "Only N-term tryptic (K/R)" = "N_only",
                                     "Only C-term tryptic (K/R)" = "C_only",
                                     "Neither N- nor C-term tryptic" = "none"),
                         selected = c("both", "N_only", "C_only", "none"), # Default to all
                         inline = TRUE),
      tags$hr(),
      
      # Valid Value Filtering
      tags$h5("Valid Value Filtering"),
      
      # Radio buttons for filter logic
      radioButtons("filter_logic", "Filter based on valid values:",
                   choices = c("In at least one group" = "one_group",
                               "In both groups" = "both_groups"),
                   selected = "one_group"),
      
      # Numeric input for minimum values
      numericInput("min_valid_values", "Minimum valid values:", 
                   value = 2, min = 1, step = 1),
      tags$hr(),
      
      # Missing Data Handling
      tags$h5("Missing Data Handling"),
      # Radio buttons to choose missing data method
      radioButtons("missing_data_method", "After filtering, handle remaining missing values by:",
                   choices = c("Keeping NAs (for Limma)" = "filter",
                               "Imputation (Down-shift)" = "impute"),
                   selected = "impute"),
      
      # Conditional panel for Imputation
      conditionalPanel(
        condition = "input$missing_data_method == 'impute'",
        p("Imputation is performed on the *filtered* data set."),
        sliderInput("impute_downshift", "Imputation Down-shift (log2 units)",
                    min = 1.0, max = 3.0, value = 1.8, step = 0.1),
        sliderInput("impute_width", "Imputation Width (sd)",
                    min = 0.1, max = 1.0, value = 0.3, step = 0.1)
      ),
      
      
      # Plotting Parameters
      tags$h4("4. Plotting Parameters"),
      
      # Axis Limits
      tags$h5("Manual Axis Limits"),
      p("Leave blank for automatic limits."),
      fluidRow(
        column(6, numericInput("xmin", "X-min", value = NA)),
        column(6, numericInput("xmax", "X-max", value = NA))
      ),
      fluidRow(
        column(6, numericInput("ymin", "Y-min", value = NA)),
        column(6, numericInput("ymax", "Y-max", value = NA))
      ),
      tags$hr(),
      
      # Fold Change Threshold
      numericInput("fc_threshold", "Log2 Fold Change Threshold", 
                   value = 3, min = 0, step = 0.1),
      
      # P-Value Threshold
      numericInput("pval_threshold", "Adjusted P-Value Threshold (-log10)", 
                   value = 1.3, min = 0, step = 0.1), # 1.3 corresponds to p=0.05
      
      # Color Pickers for Up/Down
      colourpicker::colourInput("up_color", "Up-regulated Color", value = "#d62728"), # Red
      colourpicker::colourInput("down_color", "Down-regulated Color", value = "#1f77b4"), # Blue
      
      sliderInput("dot_size_multiplier", "Dot Size Multiplier",
                  min = 0.5, max = 3, value = 2, step = 0.1),
      
      # Control for dot edge width
      sliderInput("dot_edge_width", "Dot Edge Width", 
                  min = 0.1, max = 3, value = 1.5, step = 0.1),
      
      # Default edge color (always visible)
      colourpicker::colourInput("default_edge_color", "Default Edge Color", value = "grey70"),
      
      
      # Labeling Options
      tags$h4("5. Labeling Options"),
      # Column to use for labels
      uiOutput("label_col_selector"),
      
      # Number of top hits to label
      numericInput("top_n_labels", "Label Top N Hits", 
                   value = 0, min = 0, step = 1),
      
      # Overlap Control
      sliderInput("max_overlaps", "Label Overlap Limit",
                  min = 0, max = 100, value = 20, step = 1),
      p("Lower values remove overlapping labels to reduce clutter."),
      
      # Duplicate Simplification
      checkboxInput("simplify_duplicates", "Replace all labels with numbers", value = FALSE),
      p("If checked, labels are replaced by numeric IDs and a key is added."),
      tags$hr(),
      
      # PE/PPE Highlighting
      tags$h5("PE/PPE Highlighting"),
      p("Only highlights *significant* proteins from this class."),
      checkboxInput("highlight_pe_ppe", "Highlight PE/PPE proteins?", value = TRUE),
      
      # Conditional dropdown for PE/PPE column
      conditionalPanel(
        condition = "input$highlight_pe_ppe == true",
        uiOutput("pe_ppe_col_selector"),
        colourpicker::colourInput("pe_ppe_edge_color", "PE/PPE Edge Color", value = "limegreen")
      ),
      
      # Manual Highlighting
      tags$hr(),
      tags$h5("Manual Highlighting"),
      p("Only highlights/labels *significant* genes from this list."),
      uiOutput("manual_highlight_col_selector"), # Select column to search
      uiOutput("manual_highlight_gene_selector"), # Multi-select box for genes
      colourpicker::colourInput("manual_highlight_color", "Manual Highlight Color", value = "#FFD700"), # Gold
      
      # Action button to run
      tags$hr(),
      actionButton("run_analysis", "Generate Plot", class = "btn-primary", style = "width: 100%;"),
      
      # Download buttons
      br(), br(),
      uiOutput("download_plot_button"), # For PDF
      br(),
      uiOutput("download_button") # For CSV
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("Volcano Plot", 
                 p("Plot will appear after uploading data, selecting groups, and clicking 'Generate Plot'."),
                 plotOutput("volcano_plot", height = "700px")
        ),
        tabPanel("Results Table", 
                 p("Results table will appear here."),
                 DT::dataTableOutput("results_table")
        ),
        tabPanel("Uploaded Data", 
                 p("A preview of your uploaded file."),
                 DT::dataTableOutput("contents")
        )
      )
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  
  # Reactive: Read Uploaded Data
  rawData <- reactive({
    req(input$file1)
    tryCatch({
      df <- readr::read_csv(input$file1$datapath)
      df
    }, error = function(e) {
      stop(safeError(e))
    })
  })
  
  # Render UI: Data Preview Table
  output$contents <- DT::renderDataTable({
    req(rawData())
    DT::datatable(rawData(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Render UI: Dynamic Column Selectors
  
  # Get all column names
  all_cols <- reactive({
    req(rawData())
    colnames(rawData())
  })
  
  # Get only LFQ intensity columns
  lfq_cols <- reactive({
    req(all_cols())
    grep("LFQ intensity", all_cols(), value = TRUE)
  })
  
  # Get annotation columns (non-LFQ)
  anno_cols <- reactive({
    req(all_cols())
    grep("LFQ intensity", all_cols(), value = TRUE, invert = TRUE)
  })
  
  # Selector for Group 1
  output$group1_selector <- renderUI({
    req(lfq_cols())
    checkboxGroupInput("group1_cols", "Select Group 1 Columns",
                       choices = lfq_cols(),
                       selected = lfq_cols()[1:3]) # Default to first 3
  })
  
  # Selector for Group 2
  output$group2_selector <- renderUI({
    req(lfq_cols())
    # Choices are all LFQ columns *not* in group 1
    available_cols <- setdiff(lfq_cols(), input$group1_cols)
    checkboxGroupInput("group2_cols", "Select Group 2 Columns",
                       choices = available_cols,
                       selected = available_cols[1:3]) # Default to next 3
  })
  
  # Selector for Label Column
  output$label_col_selector <- renderUI({
    req(anno_cols())
    # Try to find "Gene name" or "Gene locus" as default
    default_label <- "Gene name"
    if (!("Gene name" %in% anno_cols())) {
      default_label <- anno_cols()[1] # Fallback to first anno col
    }
    
    selectInput("label_col", "Column for Labels",
                choices = anno_cols(),
                selected = default_label)
  })
  
  # Selector for PE/PPE Column
  output$pe_ppe_col_selector <- renderUI({
    req(anno_cols())
    # Try to find "Gene name" as default
    default_pe <- "Gene name"
    if (!("Gene name" %in% anno_cols())) {
      default_pe <- anno_cols()[1] # Fallback
    }
    
    selectInput("pe_ppe_col", "Column to find PE/PPE in",
                choices = anno_cols(),
                selected = default_pe)
  })
  
  # Selector for Manual Highlight Column
  output$manual_highlight_col_selector <- renderUI({
    req(anno_cols())
    default_highlight_col <- "Gene name"
    if (!("Gene name" %in% anno_cols())) {
      default_highlight_col <- anno_cols()[1] # Fallback
    }
    selectInput("manual_col", "Column to search for genes",
                choices = anno_cols(),
                selected = default_highlight_col)
  })
  
  # Multi-select for Manual Highlight Genes
  output$manual_highlight_gene_selector <- renderUI({
    req(rawData(), input$manual_col)
    
    # Get unique, non-NA gene names from the selected column
    gene_choices <- rawData() %>%
      pull(!!sym(input$manual_col)) %>%
      unique() %>%
      na.omit() %>%
      sort()
    
    selectizeInput("manual_genes", "Select Genes to Highlight",
                   choices = gene_choices,
                   multiple = TRUE,
                   options = list(
                     placeholder = 'Type to search for a gene...',
                     plugins = list('remove_button')
                   ))
  })
  
  # Reactive: Run Differential Expression
  de_results <- eventReactive(input$run_analysis, {
    
    req(rawData(), input$group1_cols, input$group2_cols)
    
    data <- rawData()
    
    # Data Preparation
    
    # Get selected LFQ columns
    g1_cols <- input$group1_cols
    g2_cols <- input$group2_cols
    all_lfq_cols <- c(g1_cols, g2_cols)
    
    # Get annotation columns (all columns *except* LFQs)
    anno_cols_names <- setdiff(colnames(data), grep("LFQ intensity", colnames(data), value = TRUE))
    
    # Create a minimal dataframe
    df <- data[, c(anno_cols_names, all_lfq_cols)]
    
    # Tryptic Site Filtering
    tryptic_keep_rows <- rep(TRUE, nrow(df)) 
    if ("Amino acid before" %in% colnames(df) && "Last amino acid" %in% colnames(df)) {
      
      # Create logical vectors for each condition
      # Handle NAs by treating them as non-tryptic (FALSE)
      n_term_tryptic <- ifelse(is.na(df$`Amino acid before`), FALSE, df$`Amino acid before` %in% c("K", "R"))
      c_term_tryptic <- ifelse(is.na(df$`Last amino acid`), FALSE, df$`Last amino acid` %in% c("K", "R"))
      
      # Logic TRUE/FALSE
      is_both <- (n_term_tryptic) & (c_term_tryptic)
      is_N_only <- (n_term_tryptic) & (!c_term_tryptic)
      is_C_only <- (!n_term_tryptic) & (c_term_tryptic)
      is_none <- (!n_term_tryptic) & (!c_term_tryptic)
      
      # Initialize a logical vector to keep no rows
      tryptic_keep_rows <- rep(FALSE, nrow(df))
      
      # Add rows to keep based on user selection
      if ("both" %in% input$tryptic_filter) {
        tryptic_keep_rows <- tryptic_keep_rows | is_both
      }
      if ("N_only" %in% input$tryptic_filter) {
        tryptic_keep_rows <- tryptic_keep_rows | is_N_only
      }
      if ("C_only" %in% input$tryptic_filter) {
        tryptic_keep_rows <- tryptic_keep_rows | is_C_only
      }
      if ("none" %in% input$tryptic_filter) {
        tryptic_keep_rows <- tryptic_keep_rows | is_none
      }
      
      if (sum(tryptic_keep_rows) == 0) {
        stop("No proteins remaining after Tryptic Site filtering. Please select at least one tryptic type.")
      }
      
    } else if (!is.null(input$tryptic_filter) && length(input$tryptic_filter) < 4) {
      # Update warning message
      showNotification("Tryptic site columns ('Amino acid before', 'Last amino acid') not found. Skipping this filter.", 
                       type = "warning", duration = 10)
    }
   
    
    # Apply tryptic filter. df_tryptic contains RAW (non-logged) data
    df_tryptic <- df[tryptic_keep_rows, ]
    
    # Data Prep for Limma (NA, Log2, Normalize)
    df_limma <- df_tryptic
    
    # Convert 0s to NAs (limma works with NA, not 0)
    df_limma[all_lfq_cols][df_limma[all_lfq_cols] == 0] <- NA
    
    # Log2 Transform
    df_limma[all_lfq_cols] <- log2(df_limma[all_lfq_cols])
    
    # Median Normalization (Optional)
    if (input$median_normalize) {
      df_limma[all_lfq_cols] <- lapply(df_limma[all_lfq_cols], function(col) {
        col - median(col, na.rm = TRUE)
      })
    }
    
    # Valid Values Filtering
    min_val <- input$min_valid_values
    
    # Count valid (non-NA) values in each group from the log-transformed data
    valid_g1 <- rowSums(!is.na(df_limma[, g1_cols]))
    valid_g2 <- rowSums(!is.na(df_limma[, g2_cols]))
    
    # Apply filter logic based on user choice
    filter_keep_index <- NULL
    if (input$filter_logic == "one_group") {
      # Keep rows where *at least one* group has the minimum number of valid values
      filter_keep_index <- (valid_g1 >= min_val) | (valid_g2 >= min_val)
    } else {
      # Keep rows where *both* groups have the minimum number of valid values (default)
      filter_keep_index <- (valid_g1 >= min_val) & (valid_g2 >= min_val)
    }
    
    # df_filtered is the log-transformed, normalized, value-filtered data
    df_filtered <- df_limma[filter_keep_index, ]
    
    if (nrow(df_filtered) == 0) {
      stop("No proteins remaining after valid value filtering. Try changing the filter settings.")
    }
    
    # Handle Missing Data (Impute or Keep NAs)
    
    df_for_limma <- NULL
    
    if (input$missing_data_method == "filter") {
      # Method 1: Limma Filtering (Keep NAs) 
      df_for_limma <- df_filtered
      
    } else {
      # Method 2: Imputation (on the filtered data)
      expr_data <- as.data.frame(df_filtered[, all_lfq_cols])
      
      # Calculate global mean and sd from ALL log-transformed, tryptic-filtered data (before value filter)
      global_mean <- mean(as.matrix(df_limma[, all_lfq_cols]), na.rm = TRUE)
      global_sd <- sd(as.matrix(df_limma[, all_lfq_cols]), na.rm = TRUE)
      
      # Calculate parameters for the down-shifted distribution
      impute_mean <- global_mean - (input$impute_downshift * global_sd)
      impute_sd <- global_sd * input$impute_width
      
      # Get total number of missing values *in the filtered data*
      n_missing <- sum(is.na(expr_data))
      
      if (n_missing > 0) {
        # Generate random values from the new distribution
        set.seed(1) # for reproducibility
        imputed_values <- rnorm(n_missing, mean = impute_mean, sd = impute_sd)
        
        # Replace NAs with imputed values
        expr_data[is.na(expr_data)] <- imputed_values
      }
      
      # Combine back with annotation data
      df_for_limma <- cbind(df_filtered[, anno_cols_names, drop = FALSE], expr_data)
    }
    
    # Differential Expression with limma
    
    # Create the design matrix
    n_g1 <- length(g1_cols)
    n_g2 <- length(g2_cols)
    design <- model.matrix(~ 0 + factor(c(rep(1, n_g1), rep(2, n_g2))))
    colnames(design) <- c("Group1", "Group2")
    
    # Extract the expression data (now processed)
    expr_data_processed <- df_for_limma[, all_lfq_cols]
    
    # Fit the linear model
    fit <- lmFit(expr_data_processed, design)
    
    # Create the contrast matrix (Group 2 vs Group 1)
    contrast_matrix <- makeContrasts(Group2 - Group1, levels = design)
    
    # Apply contrasts
    fit_contrast <- contrasts.fit(fit, contrast_matrix)
    
    # Apply eBayes empirical Bayes moderation
    fit_bayes <- eBayes(fit_contrast)
    
    # Get the results table
    results <- topTable(fit_bayes, coef = 1, number = Inf, sort.by = "P")
    
    # Calculate Dot Size (Mean LFQ) [FIXED LOGIC]
    
    # Add row_id (rownames) to the raw tryptic-filtered data *before* filtering
    df_tryptic$row_id <- rownames(df_tryptic)
    
    # Get the row_ids of the data that passed the value filter
    # df_filtered has the correct rows, and its rownames are the row_ids
    passed_value_filter_ids <- rownames(df_filtered)
    
    # Filter the raw tryptic data to get *only* the rows that passed
    raw_for_size_df <- df_tryptic[df_tryptic$row_id %in% passed_value_filter_ids, ]
    
    # Calculate MeanLFQ on correctly filtered raw data
    raw_for_size_values <- raw_for_size_df[, all_lfq_cols]
    raw_for_size_values[raw_for_size_values == 0 | is.na(raw_for_size_values)] <- NA
    
    # Calculate means and add back to the df with row_ids
    raw_for_size_df$MeanLFQ <- rowMeans(raw_for_size_values, na.rm = TRUE)
    
    # Create a small df to merge by row_id
    size_data_to_merge <- raw_for_size_df[, c("row_id", "MeanLFQ")]
    
    
    # Format Final Table
    
    # Get row_ids from limma results and annotation data
    results$row_id <- rownames(results) # 'results' from topTable
    df_for_limma$row_id <- rownames(df_for_limma) # 'df_for_limma'
    
    # Merge annotations
    anno_data_to_merge <- df_for_limma[, c("row_id", anno_cols_names)] %>% distinct()
    
    final_results <- merge(anno_data_to_merge, 
                           results, 
                           by = "row_id")
    
    # Merge MeanLFQ data
    final_results <- merge(final_results,
                           size_data_to_merge,
                           by = "row_id")
    
    # Clean up and add helpful columns
    final_results <- final_results %>%
      select(-row_id) %>% # Remove row_id after all merges
      rename(
        logFC = logFC,
        P.Value = P.Value,
        adj.P.Val = adj.P.Val
      ) %>%
      mutate(
        log10P = -log10(adj.P.Val)
      )
    
    # Add DotSize 
    final_results$DotSize <- sqrt(final_results$MeanLFQ)
    final_results$DotSize[is.na(final_results$DotSize)] <- 0
    
    
    # Add Significance Column
    fc_thresh <- input$fc_threshold
    p_thresh <- input$pval_threshold
    
    final_results$Significance <- "Not significant"
    final_results$Significance[final_results$logFC > fc_thresh & final_results$log10P > p_thresh] <- "Up"
    final_results$Significance[final_results$logFC < -fc_thresh & final_results$log10P > p_thresh] <- "Down"
    
    # Add Highlighting Columns
    
    # PE/PPE Highlighting
    final_results$Is_PE_PPE <- FALSE
    if (input$highlight_pe_ppe && input$pe_ppe_col %in% colnames(final_results)) {
      pe_ppe_col_values <- final_results[[input$pe_ppe_col]]
      # Add !is.na() check
      final_results$Is_PE_PPE <- !is.na(pe_ppe_col_values) & grepl("^(PE|PPE)", pe_ppe_col_values, ignore.case = TRUE)
    }
    
    # Manual Highlighting
    final_results$Is_Manual_Highlight <- FALSE
    if (!is.null(input$manual_genes) && input$manual_col %in% colnames(final_results)) {
      manual_col_values <- final_results[[input$manual_col]]
      # Add !is.na() check
      final_results$Is_Manual_Highlight <- !is.na(manual_col_values) & (manual_col_values %in% input$manual_genes)
    }
    
    # Create Plotting Groups
    # Create new columns:
    # 1. 'Fill_Group' controls the dot's fill color
    # 2. 'Edge_Group' controls the dot's edge color
    
    # 1. Create the new Fill_Group (4 levels)
    final_results$Fill_Group <- final_results$Significance # ("Up", "Down", "Not significant")
    # Overwrite with "Manual Highlight" if criteria are met (and significant)
    final_results$Fill_Group[final_results$Is_Manual_Highlight & final_results$Significance != "Not significant"] <- "Manual Highlight"
    
    # 2. Create the (renamed) Edge_Group (2 levels)
    final_results$Edge_Group <- "Default"
    # MODIFICATION: Only apply PE/PPE edge if the protein is ALSO significant
    final_results$Edge_Group[final_results$Is_PE_PPE & final_results$Significance != "Not significant"] <- "PE/PPE"
    
    
    # Create Sort Key for plotting order, higher numbers are plotted on top
    final_results$Sort_Key <- 1 # Default: Not significant
    
    # Significant (Up or Down)
    final_results$Sort_Key[final_results$Significance != "Not significant"] <- 2
    
    # PE/PPE (and significant)
    final_results$Sort_Key[final_results$Edge_Group == "PE/PPE"] <- 3 # Edge_Group already includes significance check
    
    # Manual Highlight (and significant)
    final_results$Sort_Key[final_results$Fill_Group == "Manual Highlight"] <- 4 # Fill_Group already includes significance check
    
    
    return(final_results)
  })
  
  # Reactive: Build the Plot
  # renderPlot and downloadHandler
  plotInput <- reactive({
    req(de_results(), input$label_col, input$manual_col)
    
    df <- de_results()
    
    # Sort data for plotting
    # Create a factor for Fill_Group *for the legend order*.
    plot_levels <- c("Not significant", "Down", "Up", "Manual Highlight")
    
    # Re-order the data frame based on the new Sort_Key for plotting z-index
    # Mutate Fill_Group factor for correct legend order
    df_for_plot <- df %>%
      mutate(Fill_Group = factor(.data$Fill_Group, levels = plot_levels)) %>%
      arrange(.data$Sort_Key) 
    
    # Get thresholds for lines
    fc_thresh <- input$fc_threshold
    p_thresh <- input$pval_threshold
    
    # Get plot parameters
    size_mult <- input$dot_size_multiplier
    dot_edge_wid <- input$dot_edge_width 
    
    # Get colors
    default_edge_col <- input$default_edge_color
    pe_ppe_edge_col <- input$pe_ppe_edge_color
    manual_fill_col <- input$manual_highlight_color # This is a fill color now
    
    # Get Up/Down Colors from input
    up_reg_col <- input$up_color
    down_reg_col <- input$down_color
    
    # Define the colors for the edges
    edge_colors <- c(
      "Default" = default_edge_col,
      "PE/PPE" = pe_ppe_edge_col
    )
    
    # Define the colors for the fill using input variables
    fill_colors <- c(
      "Up" = up_reg_col, # Use input
      "Down" = down_reg_col, # Use input
      "Not significant" = "grey",
      "Manual Highlight" = manual_fill_col
    )
    
    # Prepare data for labels
    
    # Get top N significant points to label
    df_top_n <- df %>% 
      filter(Significance != "Not significant") %>%
      arrange(desc(abs(logFC)), desc(log10P)) %>%
      head(input$top_n_labels)
    
    # Get manually highlighted points to label (must be significant)
    df_manual <- df %>% 
      filter(Is_Manual_Highlight == TRUE & Significance != "Not significant")
    
    # Combine them and remove duplicates
    df_for_labels <- bind_rows(df_top_n, df_manual) %>%
      distinct() 
    
    # Filter out any rows where the selected label column is NA
    if (nrow(df_for_labels) > 0 && input$label_col %in% colnames(df_for_labels)) {
      df_for_labels <- df_for_labels %>%
        filter(!is.na(!!sym(input$label_col)))
    }
    
    # Simplify Labels Logic (All labels)
    legend_text <- NULL
    
    if (input$simplify_duplicates && nrow(df_for_labels) > 0) {
      
      # Get ALL unique labels that are about to be shown
      unique_labels <- sort(unique(df_for_labels[[input$label_col]]))
      
      # Create a mapping: Label -> Numeric ID
      # seq_along gives 1, 2, 3... corresponding to the sorted labels
      label_map <- setNames(seq_along(unique_labels), unique_labels)
      
      # Create legend string (e.g., "1 = GeneA\n2 = GeneB")
      legend_items <- paste(label_map, "=", names(label_map))
      legend_text <- paste(legend_items, collapse = "\n")
      
      # Create a new column with the numeric IDs
      # Apply the mapping to the original labels
      df_for_labels$Label_To_Show <- as.character(label_map[df_for_labels[[input$label_col]]])
      
      # Use this new column for the plot
      label_column_name <- "Label_To_Show" 
      
    } else {
      # Feature off, use original column
      label_column_name <- input$label_col
    }
    
    
    # Create the plot
    p <- ggplot(df_for_plot, aes(x = logFC, y = log10P, fill = Fill_Group, size = DotSize, color = Edge_Group)) +
      # Updated geom_point with stroke
      geom_point(shape = 21, alpha = 0.7, stroke = dot_edge_wid) +
      
      # Add labels using ggrepel 
      geom_text_repel(data = df_for_labels, 
                      aes(label = !!sym(label_column_name)), # Use the dynamic column name
                      size = 4, 
                      color = "black", # Label color is always black
                      max.overlaps = input$max_overlaps, # Use slider input
                      
                      # Ensure lines connect to points
                      min.segment.length = 0, # Always draw line
                      segment.color = "black",
                      segment.size = 0.5,
                      
                      box.padding = 0.5,
                      point.padding = 0.5) +
      
      # Add threshold lines (hardcoded style)
      geom_hline(yintercept = p_thresh, linetype = "dashed", 
                 color = "grey50", linewidth = 0.5) +
      geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", 
                 color = "grey50", linewidth = 0.5) +
      
      # Scales and theme
      # Use the fill_colors vector
      scale_fill_manual(name = "Significance / Highlight",
                        values = fill_colors) +
      
      scale_size(name = "Mean LFQ (sqrt)", 
                 range = c(1 * size_mult, 10 * size_mult)) +
      
      # Use the edge_colors vector
      scale_color_manual(name = "PE/PPE Edge",
                         values = edge_colors) +
      
      # Add Coordinate Limits
      # ggplot's coord_cartesian handles NA values gracefully,
      # so it will only apply limits that are not NA.
      coord_cartesian(
        xlim = c(input$xmin, input$xmax),
        ylim = c(input$ymin, input$ymax) # 
      ) +
      
      labs(
        title = paste("Group 2 vs Group 1"),
        x = "Log2 Fold Change",
        y = "-log10 (Adjusted P-Value)"
      ) +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "bottom",
        # Remove gridlines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      # Ensure legends are combined where possible
      guides(fill = guide_legend(title = "Significance / Highlight", override.aes = list(size = 5)),
             color = guide_legend(title = "PE/PPE Edge", override.aes = list(size = 5)),
             size = "none") # Hide the size legend
    
    # Add Annotation for Duplicate Key
    if (!is.null(legend_text)) {
      # Add the text annotation to the plot
      # We place it in the top-left or top-right, or outside. 
      # Here, we'll try adding it as a caption-like element using annotate or labs(caption=...)
      # labs(caption) is cleaner.
      p <- p + labs(caption = paste("Label Key:\n", legend_text)) +
        theme(plot.caption = element_text(hjust = 0, size = 10, face = "italic"))
    }
    
    return(p)
  })
  
  # Render: Volcano Plot
  output$volcano_plot <- renderPlot({
    # Call the reactive plot
    plotInput()
  })
  
  # Render: Results Table
  output$results_table <- DT::renderDataTable({
    req(de_results())
    
    df_for_table <- de_results() %>%
      # Select relevant columns for display
      select(
        !!sym(input$label_col), # Use the selected label column
        !!sym(input$manual_col), # Also show the manual highlight column
        logFC, 
        P.Value, 
        adj.P.Val, 
        log10P, 
        MeanLFQ, 
        Significance,
        Fill_Group,
        Edge_Group,
        Sort_Key, # Show the sort key
        everything() # Put the rest of the columns after
      ) %>%
      # De-duplicate
      distinct(across(everything()), .keep_all = TRUE) %>%
      # Round numeric columns for cleaner display
      mutate(across(where(is.numeric), ~ round(., 4)))
    
    DT::datatable(df_for_table, 
                  extensions = 'Buttons',
                  options = list(pageLength = 10, scrollX = TRUE, 
                                 dom = 'Bfrtip',
                                 buttons = c('copy', 'csv', 'excel')))
  })
  
  # Download Handler: Plot 
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("volcano_plot-", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(plotInput())
      # Save the plot to a PDF file
      ggsave(file, plot = plotInput(), device = "pdf", width = 10, height = 12, units = "in", dpi = 300)
    }
  )
  
  # Download Handler: Results
  output$download_results <- downloadHandler(
    filename = function() {
      paste("differential_expression_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(de_results())
      # Write the full results table to CSV
      write.csv(de_results(), file, row.names = FALSE)
    }
  )
  
  # Render UI: Download Plot Button
  output$download_plot_button <- renderUI({
    req(de_results())
    downloadButton("download_plot", "Download Plot (PDF)", style = "width: 100%;")
  })
  
  # Render UI: Download Results Button
  output$download_button <- renderUI({
    req(de_results())
    # Note: The 'Download Table' button is now part of the DT::datatable itself (dom = 'Bfrtip')
    # This button provides a separate download for the *full* unfiltered dataset.
    downloadButton("download_results", "Download Full Results (CSV)", style = "width: 100%;")
  })
  
}

# Run the Application
shinyApp(ui = ui, server = server)
