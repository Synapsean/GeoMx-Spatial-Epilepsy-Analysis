# =============================================================================
# GeoMx Chronic KA Analysis - Interactive Dashboard
# =============================================================================
# Launch: shiny::runApp("shiny_app")

library(shiny)
library(DT)
library(plotly)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(grDevices)

# =============================================================================
# LOAD DATA
# =============================================================================
results_path <- file.path("..", "results", "de_results_all.rds")
if (!file.exists(results_path)) {
  # Try alternative path (if launched from project root)
  results_path <- file.path("results", "de_results_all.rds")
}

de_results <- readRDS(results_path)

# Parse comparison names into metadata
parse_comparison <- function(name) {
  # Naming: Type_[Treatment]_[Region]_[Side]_CellType
  # Trt = Treatment effect (KA vs PBS)
  # Reg = Regional effect (CA3 vs CA1)
  # Side = Side effect (Ipsi vs Contra)
  # Int = Interaction (Region x Treatment)
  
  parts <- strsplit(name, "_")[[1]]
  type <- parts[1]
  
  contrast <- switch(type,
    "Trt" = "KA vs PBS",
    "Reg" = "CA3 vs CA1",
    "Side" = "Ipsi vs Contra",
    "Int" = "Region x Treatment",
    "Unknown"
  )
  
  # Extract cell type (always last part)
  cell_type <- parts[length(parts)]
  cell_type <- switch(cell_type,
    "Astro" = "Astrocyte",
    "Micro" = "Microglia",
    "Neuron" = "Neuron",
    cell_type
  )
  
  # Extract other info from middle parts
  middle <- parts[2:(length(parts)-1)]
  
  region <- ifelse(any(middle %in% c("CA3", "CA1")), 
                   middle[middle %in% c("CA3", "CA1")][1], "All")
  side <- ifelse(any(middle %in% c("Ipsi", "Contra")), 
                 middle[middle %in% c("Ipsi", "Contra")][1], "All")
  treatment <- ifelse(any(middle %in% c("KA", "PBS")), 
                      middle[middle %in% c("KA", "PBS")][1], "All")
  
  data.frame(
    name = name,
    contrast = contrast,
    cell_type = cell_type,
    region = region,
    side = side,
    treatment = treatment,
    stringsAsFactors = FALSE
  )
}

comparison_meta <- do.call(rbind, lapply(names(de_results), parse_comparison))

# Build summary stats
comparison_meta$n_fdr_sig <- sapply(de_results, function(x) sum(x$adj.P.Val < 0.05))
comparison_meta$n_nominal_sig <- sapply(de_results, function(x) sum(x$P.Value < 0.05))
comparison_meta$top_gene <- sapply(de_results, function(x) rownames(x)[1])
comparison_meta$top_logFC <- sapply(de_results, function(x) round(x$logFC[1], 2))
comparison_meta$top_fdr <- sapply(de_results, function(x) signif(x$adj.P.Val[1], 3))
comparison_meta$n_genes <- sapply(de_results, nrow)

# Get all gene names
all_genes <- sort(unique(unlist(lapply(de_results, rownames))))

# Build a master gene table (gene x comparison)
build_gene_summary <- function(gene_name) {
  results <- lapply(names(de_results), function(comp) {
    de <- de_results[[comp]]
    if (gene_name %in% rownames(de)) {
      row <- de[gene_name, ]
      data.frame(
        Comparison = comp,
        logFC = round(row$logFC, 3),
        P.Value = signif(row$P.Value, 3),
        FDR = signif(row$adj.P.Val, 3),
        AveExpr = round(row$AveExpr, 2),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
  do.call(rbind, results)
}

# Readable labels
readable_label <- function(name) {
  meta <- comparison_meta[comparison_meta$name == name, ]
  parts <- c()
  
  if (meta$treatment != "All") parts <- c(parts, meta$treatment)
  if (meta$region != "All") parts <- c(parts, meta$region)
  if (meta$side != "All") parts <- c(parts, meta$side)
  parts <- c(parts, meta$cell_type)
  
  prefix <- switch(meta$contrast,
    "KA vs PBS" = "KA vs PBS:",
    "CA3 vs CA1" = "CA3 vs CA1:",
    "Ipsi vs Contra" = "Ipsi vs Contra:",
    "Region x Treatment" = "Region x Trt:",
    ""
  )
  
  paste(prefix, paste(parts, collapse = " "))
}

comparison_meta$label <- sapply(comparison_meta$name, readable_label)

# =============================================================================
# Custom value box (simple version without shinydashboard)
# =============================================================================
valueBoxUI <- function(id) {
  uiOutput(NS(id, "box"))
}

valueBoxServer <- function(id, value, subtitle, color = "#3c8dbc") {
  moduleServer(id, function(input, output, session) {
    output$box <- renderUI({
      div(style = paste0("background:", color, 
                         "; color:white; padding:15px; border-radius:5px; text-align:center; margin:10px 0;"),
        h2(value, style = "margin:0;"),
        p(subtitle, style = "margin:0; font-size:14px;")
      )
    })
  })
}

# =============================================================================
# UI
# =============================================================================
ui <- navbarPage(
  title = "GeoMx Chronic KA Analysis",
  theme = NULL,
  
  # --- TAB 1: Overview ---
  tabPanel("Overview",
    fluidRow(
      column(12,
        h3("Analysis Summary"),
        p("60 differential expression comparisons using mixed-effects models (dream) with Mouse as random effect."),
        p("Fallback to limma-trend when insufficient observations for random effects.")
      )
    ),
    fluidRow(
      column(4, valueBoxUI("n_comparisons")),
      column(4, valueBoxUI("n_genes_tested")),
      column(4, valueBoxUI("n_with_hits"))
    ),
    fluidRow(
      column(12,
        h4("Significant Genes per Comparison (FDR < 0.05)"),
        plotlyOutput("overview_barplot", height = "600px")
      )
    ),
    fluidRow(
      column(12,
        h4("Full Comparison Summary"),
        DTOutput("overview_table")
      )
    )
  ),
  
  # --- TAB 2: Comparison Explorer ---
  tabPanel("Comparison Explorer",
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Filter Comparisons"),
        selectInput("explore_contrast", "Contrast Type:",
                    choices = c("All", "KA vs PBS", "CA3 vs CA1", 
                               "Ipsi vs Contra", "Region x Treatment"),
                    selected = "All"),
        selectInput("explore_cell", "Cell Type:",
                    choices = c("All", "Astrocyte", "Microglia", "Neuron"),
                    selected = "All"),
        selectInput("explore_region", "Region:",
                    choices = c("All", "CA3", "CA1"),
                    selected = "All"),
        selectInput("explore_side", "Side:",
                    choices = c("All", "Ipsi", "Contra"),
                    selected = "All"),
        hr(),
        selectInput("selected_comparison", "Select Comparison:",
                    choices = NULL),
        hr(),
        sliderInput("fdr_threshold", "FDR Threshold:",
                    min = 0.01, max = 0.25, value = 0.05, step = 0.01),
        sliderInput("lfc_threshold", "Log2FC Threshold:",
                    min = 0, max = 3, value = 0.5, step = 0.1),
        hr(),
        downloadButton("download_de_table", "Download DE Table (CSV)")
      ),
      mainPanel(width = 9,
        fluidRow(
          column(6,
            h4("Volcano Plot"),
            plotlyOutput("volcano_plot", height = "500px")
          ),
          column(6,
            h4("MA Plot"),
            plotlyOutput("ma_plot", height = "500px")
          )
        ),
        fluidRow(
          column(12,
            h4("Differentially Expressed Genes"),
            DTOutput("de_table")
          )
        )
      )
    )
  ),
  
  # --- TAB 3: Gene Search ---
  tabPanel("Gene Search",
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Search for a Gene"),
        selectizeInput("gene_search", "Gene Name:",
                       choices = NULL,
                       options = list(
                         placeholder = "Type gene name...",
                         maxOptions = 50
                       )),
        helpText("Start typing to search. Shows the gene's results across all 60 comparisons."),
        hr(),
        h4("Or Paste a Gene List"),
        textAreaInput("gene_list", "Gene List (one per line):",
                      rows = 6,
                      placeholder = "Gfap\nMbp\nBdnf\nIl1b"),
        actionButton("search_list", "Search Gene List", class = "btn-primary"),
        hr(),
        downloadButton("download_gene_results", "Download Results (CSV)")
      ),
      mainPanel(width = 9,
        conditionalPanel(
          condition = "input.gene_search != '' && input.gene_search != null",
          fluidRow(
            column(12,
              h4(textOutput("gene_title")),
              plotlyOutput("gene_logfc_plot", height = "400px")
            )
          ),
          fluidRow(
            column(12,
              h4("Results Across All Comparisons"),
              DTOutput("gene_table")
            )
          )
        ),
        conditionalPanel(
          condition = "input.search_list > 0",
          fluidRow(
            column(12,
              h4("Gene List Heatmap"),
              plotOutput("gene_list_heatmap", height = "600px")
            )
          ),
          fluidRow(
            column(12,
              h4("Gene List Results"),
              DTOutput("gene_list_table")
            )
          )
        )
      )
    )
  ),
  
  # --- TAB 4: Compare Across ---
  tabPanel("Cross-Comparison",
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Compare Cell Types"),
        selectInput("cross_contrast", "Contrast Type:",
                    choices = c("KA vs PBS", "CA3 vs CA1", "Ipsi vs Contra"),
                    selected = "KA vs PBS"),
        selectInput("cross_region", "Region:",
                    choices = c("All", "CA3", "CA1"),
                    selected = "All"),
        selectInput("cross_side", "Side:",
                    choices = c("All", "Ipsi", "Contra"),
                    selected = "All"),
        selectInput("cross_treatment", "Treatment:",
                    choices = c("All", "KA", "PBS"),
                    selected = "All"),
        hr(),
        numericInput("cross_top_n", "Top N genes per comparison:", 
                     value = 20, min = 5, max = 100)
      ),
      mainPanel(width = 9,
        fluidRow(
          column(12,
            h4("Overlap of Significant Genes Across Cell Types"),
            plotOutput("venn_plot", height = "400px")
          )
        ),
        fluidRow(
          column(12,
            h4("Top Gene Heatmap"),
            plotOutput("cross_heatmap", height = "500px")
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {
  
  # --- Value Boxes ---
  valueBoxServer("n_comparisons", 
                 length(de_results), 
                 "Comparisons Run", "#3c8dbc")
  valueBoxServer("n_genes_tested", 
                 comparison_meta$n_genes[1], 
                 "Genes Tested", "#00a65a")
  valueBoxServer("n_with_hits",
                 sum(comparison_meta$n_fdr_sig > 0),
                 "Comparisons with FDR < 0.05 Hits", "#f39c12")
  
  # --- Overview Bar Plot ---
  output$overview_barplot <- renderPlotly({
    df <- comparison_meta %>%
      arrange(contrast, desc(n_fdr_sig)) %>%
      mutate(label = factor(label, levels = rev(label)))
    
    color_map <- c("KA vs PBS" = "#e74c3c", "CA3 vs CA1" = "#3498db",
                   "Ipsi vs Contra" = "#2ecc71", "Region x Treatment" = "#9b59b6")
    bar_colors <- color_map[df$contrast]
    
    plot_ly(df, x = ~n_fdr_sig, y = ~label, type = "bar", orientation = "h",
            marker = list(color = bar_colors),
            text = ~paste0(label, "<br>FDR sig: ", n_fdr_sig,
                           "<br>Nominal sig: ", n_nominal_sig,
                           "<br>Top gene: ", top_gene),
            hoverinfo = "text") %>%
      layout(
        xaxis = list(title = "Number of Significant Genes (FDR < 0.05)"),
        yaxis = list(title = "", tickfont = list(size = 9)),
        height = 600,
        showlegend = FALSE,
        margin = list(l = 250)
      )
  })
  
  # --- Overview Table ---
  output$overview_table <- renderDT({
    df <- comparison_meta %>%
      select(Label = label, Contrast = contrast, `Cell Type` = cell_type,
             Region = region, Side = side, Treatment = treatment,
             `FDR Sig` = n_fdr_sig, `Nominal Sig` = n_nominal_sig,
             `Top Gene` = top_gene, `Top logFC` = top_logFC, `Top FDR` = top_fdr)
    
    datatable(df, filter = "top", rownames = FALSE,
              options = list(pageLength = 20, scrollX = TRUE)) %>%
      formatStyle("FDR Sig", 
                  backgroundColor = styleInterval(c(0, 10, 50), 
                                                   c("white", "#fff3cd", "#ffc107", "#ff5722")))
  })
  
  # --- Comparison Explorer: Filter logic ---
  filtered_comparisons <- reactive({
    df <- comparison_meta
    if (input$explore_contrast != "All") df <- df[df$contrast == input$explore_contrast, ]
    if (input$explore_cell != "All") df <- df[df$cell_type == input$explore_cell, ]
    if (input$explore_region != "All") df <- df[df$region == input$explore_region, ]
    if (input$explore_side != "All") df <- df[df$side == input$explore_side, ]
    df
  })
  
  observe({
    df <- filtered_comparisons()
    choices <- setNames(df$name, df$label)
    updateSelectInput(session, "selected_comparison", choices = choices)
  })
  
  # Update gene search choices (server-side for performance)
  updateSelectizeInput(session, "gene_search", 
                       choices = all_genes, 
                       server = TRUE)
  
  # Current DE result
  current_de <- reactive({
    req(input$selected_comparison)
    de <- de_results[[input$selected_comparison]]
    de$Gene <- rownames(de)
    de$Significant <- ifelse(de$adj.P.Val < input$fdr_threshold & 
                              abs(de$logFC) > input$lfc_threshold, 
                            "Significant", "Not Significant")
    de
  })
  
  # --- Volcano Plot ---
  output$volcano_plot <- renderPlotly({
    de <- current_de()
    
    de$neg_log10_p <- -log10(de$P.Value)
    de$neg_log10_p[is.infinite(de$neg_log10_p)] <- max(de$neg_log10_p[is.finite(de$neg_log10_p)]) + 1
    
    colors <- ifelse(de$Significant == "Significant", "#e74c3c", "#bdc3c7")
    
    plot_ly(de, x = ~logFC, y = ~neg_log10_p,
            type = "scatter", mode = "markers",
            marker = list(color = colors, size = 5, opacity = 0.6),
            text = ~paste0("Gene: ", Gene,
                           "<br>logFC: ", round(logFC, 3),
                           "<br>P.Value: ", signif(P.Value, 3),
                           "<br>FDR: ", signif(adj.P.Val, 3)),
            hoverinfo = "text") %>%
      layout(
        xaxis = list(title = "Log2 Fold Change",
                     zeroline = TRUE,
                     zerolinecolor = "#aaaaaa"),
        yaxis = list(title = "-Log10(P-value)"),
        shapes = list(
          list(type="line", x0=input$lfc_threshold,  x1=input$lfc_threshold,
               y0=0, y1=1, yref="paper", line=list(dash="dash", color="#555", width=1)),
          list(type="line", x0=-input$lfc_threshold, x1=-input$lfc_threshold,
               y0=0, y1=1, yref="paper", line=list(dash="dash", color="#555", width=1)),
          list(type="line", x0=0, x1=1, xref="paper",
               y0=-log10(input$fdr_threshold), y1=-log10(input$fdr_threshold),
               line=list(dash="dash", color="#555", width=1))
        ),
        showlegend = FALSE
      )
  })
  
  # --- MA Plot ---
  output$ma_plot <- renderPlotly({
    de <- current_de()
    
    colors <- ifelse(de$Significant == "Significant", "#e74c3c", "#bdc3c7")
    
    plot_ly(de, x = ~AveExpr, y = ~logFC,
            type = "scatter", mode = "markers",
            marker = list(color = colors, size = 5, opacity = 0.6),
            text = ~paste0("Gene: ", Gene,
                           "<br>AveExpr: ", round(AveExpr, 2),
                           "<br>logFC: ", round(logFC, 3),
                           "<br>FDR: ", signif(adj.P.Val, 3)),
            hoverinfo = "text") %>%
      layout(
        xaxis = list(title = "Average Expression (log2)"),
        yaxis = list(title = "Log2 Fold Change", zeroline = TRUE, zerolinecolor = "#aaaaaa"),
        showlegend = FALSE
      )
  })
  
  # --- DE Table ---
  output$de_table <- renderDT({
    de <- current_de() %>%
      select(Gene, logFC, AveExpr, P.Value, FDR = adj.P.Val, Significant) %>%
      mutate(logFC = round(logFC, 3),
             AveExpr = round(AveExpr, 2),
             P.Value = signif(P.Value, 3),
             FDR = signif(FDR, 3))
    
    datatable(de, filter = "top", rownames = FALSE,
              options = list(pageLength = 20, scrollX = TRUE)) %>%
      formatStyle("Significant",
                  backgroundColor = styleEqual("Significant", "#ffe0e0"))
  })
  
  # --- Download DE Table ---
  output$download_de_table <- downloadHandler(
    filename = function() {
      paste0("DE_", input$selected_comparison, ".csv")
    },
    content = function(file) {
      write.csv(current_de(), file, row.names = FALSE)
    }
  )
  
  # --- Gene Search ---
  gene_results <- reactive({
    req(input$gene_search)
    build_gene_summary(input$gene_search)
  })
  
  output$gene_title <- renderText({
    paste0("Results for: ", input$gene_search)
  })
  
  output$gene_logfc_plot <- renderPlotly({
    df <- gene_results()
    req(nrow(df) > 0)
    
    df <- merge(df, comparison_meta[, c("name", "contrast", "cell_type", "label")],
                by.x = "Comparison", by.y = "name")
    df$sig <- ifelse(df$FDR < 0.05, "FDR < 0.05", "Not Significant")
    df <- df %>% arrange(contrast, logFC)
    df$label <- factor(df$label, levels = df$label)
    
    color_map <- c("Astrocyte" = "#e74c3c", "Microglia" = "#3498db", "Neuron" = "#2ecc71")
    bar_colors <- color_map[df$cell_type]
    bar_opacity <- ifelse(df$sig == "FDR < 0.05", 1, 0.3)
    
    plot_ly(df, x = ~logFC, y = ~label, type = "bar", orientation = "h",
            marker = list(color = bar_colors, opacity = bar_opacity),
            text = ~paste0(label, "<br>logFC: ", round(logFC, 3), "<br>FDR: ", signif(FDR, 3)),
            hoverinfo = "text") %>%
      layout(
        xaxis = list(title = "Log2 Fold Change", zeroline = TRUE, zerolinecolor = "#aaa"),
        yaxis = list(title = "", tickfont = list(size = 8)),
        height = 400,
        showlegend = FALSE,
        margin = list(l = 250)
      )
  })
  
  output$gene_table <- renderDT({
    df <- gene_results()
    req(nrow(df) > 0)
    
    df <- merge(df, comparison_meta[, c("name", "label", "contrast", "cell_type")], 
                by.x = "Comparison", by.y = "name")
    df$Significant <- ifelse(df$FDR < 0.05, "Yes", "")
    df <- df %>%
      select(Label = label, Contrast = contrast, `Cell Type` = cell_type,
             logFC, P.Value, FDR, AveExpr, Significant) %>%
      arrange(P.Value)
    
    datatable(df, filter = "top", rownames = FALSE,
              options = list(pageLength = 20)) %>%
      formatStyle("Significant",
                  backgroundColor = styleEqual("Yes", "#ffe0e0"))
  })
  
  # --- Gene List Search ---
  gene_list_results <- eventReactive(input$search_list, {
    genes <- trimws(unlist(strsplit(input$gene_list, "\n")))
    genes <- genes[genes != ""]
    genes <- genes[genes %in% all_genes]
    
    if (length(genes) == 0) return(NULL)
    
    results <- lapply(genes, function(g) {
      df <- build_gene_summary(g)
      if (!is.null(df)) df$Gene <- g
      df
    })
    do.call(rbind, results)
  })
  
  output$gene_list_heatmap <- renderPlot({
    df <- gene_list_results()
    req(df)
    
    # Create logFC matrix
    mat <- df %>%
      select(Gene, Comparison, logFC) %>%
      pivot_wider(names_from = Comparison, values_from = logFC) %>%
      tibble::column_to_rownames("Gene") %>%
      as.matrix()
    
    # Limit to readable size
    if (ncol(mat) > 30) {
      # Keep only comparisons with significant results
      sig_comps <- df %>%
        filter(FDR < 0.05) %>%
        pull(Comparison) %>%
        unique()
      if (length(sig_comps) > 0) {
        mat <- mat[, colnames(mat) %in% sig_comps, drop = FALSE]
      }
    }
    
    # Use readable labels
    colnames(mat) <- sapply(colnames(mat), function(x) {
      lbl <- comparison_meta$label[comparison_meta$name == x]
      if (length(lbl) > 0) lbl else x
    })
    
    if (ncol(mat) > 1 && nrow(mat) > 1) {
      pheatmap(mat,
               color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
               breaks = seq(-max(abs(mat), na.rm = TRUE), max(abs(mat), na.rm = TRUE), length.out = 101),
               cluster_rows = nrow(mat) > 2,
               cluster_cols = ncol(mat) > 2,
               fontsize_row = 10,
               fontsize_col = 7,
               angle_col = 45,
               main = "Log2 Fold Change Heatmap")
    }
  })
  
  output$gene_list_table <- renderDT({
    df <- gene_list_results()
    req(df)
    
    df <- merge(df, comparison_meta[, c("name", "label", "contrast")],
                by.x = "Comparison", by.y = "name")
    df$Significant <- ifelse(df$FDR < 0.05, "Yes", "")
    df <- df %>%
      select(Gene, Label = label, Contrast = contrast, logFC, P.Value, FDR, Significant) %>%
      arrange(Gene, P.Value)
    
    datatable(df, filter = "top", rownames = FALSE,
              options = list(pageLength = 25))
  })
  
  # --- Download Gene Results ---
  output$download_gene_results <- downloadHandler(
    filename = function() {
      gene <- ifelse(!is.null(input$gene_search) && input$gene_search != "",
                     input$gene_search, "gene_list")
      paste0("Gene_", gene, "_results.csv")
    },
    content = function(file) {
      if (!is.null(input$gene_search) && input$gene_search != "") {
        df <- gene_results()
      } else {
        df <- gene_list_results()
      }
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # --- Cross-Comparison Tab ---
  cross_comparisons <- reactive({
    df <- comparison_meta %>%
      filter(contrast == input$cross_contrast)
    
    if (input$cross_region != "All") df <- df[df$region == input$cross_region, ]
    if (input$cross_side != "All") df <- df[df$side == input$cross_side, ]
    if (input$cross_treatment != "All") df <- df[df$treatment == input$cross_treatment, ]
    
    df
  })
  
  output$venn_plot <- renderPlot({
    df <- cross_comparisons()
    req(nrow(df) >= 2)
    
    # Get significant genes per cell type
    sig_genes <- lapply(split(df, df$cell_type), function(ct_df) {
      genes <- lapply(ct_df$name, function(comp) {
        de <- de_results[[comp]]
        rownames(de)[de$adj.P.Val < 0.05]
      })
      unique(unlist(genes))
    })
    
    # Simple overlap visualization (bar chart style since Venn for >2 is messy)
    overlap_df <- data.frame(
      CellType = names(sig_genes),
      Count = sapply(sig_genes, length)
    )
    
    colors <- c("Astrocyte" = "#e74c3c", "Microglia" = "#3498db", "Neuron" = "#2ecc71")
    
    p <- ggplot(overlap_df, aes(x = CellType, y = Count, fill = CellType)) +
      geom_col() +
      scale_fill_manual(values = colors) +
      labs(x = "", y = "Number of Significant Genes (FDR < 0.05)",
           title = paste("Significant Genes by Cell Type:", input$cross_contrast)) +
      theme_minimal() +
      theme(text = element_text(size = 14))
    
    # Add overlap info if we have multiple cell types
    if (length(sig_genes) >= 2) {
      all_sig <- Reduce(intersect, sig_genes)
      p <- p + annotate("text", x = 1.5, y = max(overlap_df$Count) * 0.9,
                        label = paste("Shared across all:", length(all_sig)),
                        size = 5, fontface = "italic")
    }
    
    print(p)
  })
  
  output$cross_heatmap <- renderPlot({
    df <- cross_comparisons()
    req(nrow(df) >= 2)
    
    # Get top genes from each comparison
    top_genes <- unique(unlist(lapply(df$name, function(comp) {
      de <- de_results[[comp]]
      head(rownames(de), input$cross_top_n)
    })))
    
    # Build logFC matrix
    mat <- sapply(df$name, function(comp) {
      de <- de_results[[comp]]
      de[top_genes, "logFC"]
    })
    rownames(mat) <- top_genes
    colnames(mat) <- df$label
    
    # Remove NA rows
    mat <- mat[complete.cases(mat), , drop = FALSE]
    
    if (nrow(mat) > 1 && ncol(mat) > 1) {
      pheatmap(mat,
               color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
               breaks = seq(-max(abs(mat), na.rm = TRUE), max(abs(mat), na.rm = TRUE), length.out = 101),
               fontsize_row = 8,
               fontsize_col = 9,
               angle_col = 45,
               main = paste("Top", input$cross_top_n, "Genes per Comparison (logFC)"))
    }
  })
}

# =============================================================================
# RUN
# =============================================================================
shinyApp(ui = ui, server = server)
