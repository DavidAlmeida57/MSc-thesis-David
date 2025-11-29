
#!/usr/bin/env Rscript
# pathway_plots_script.R
# Purpose: load G:Profiler pathway results, prepare data and produce dot/segment plots.

# ---------------------------
# 0) Setup - packages & paths
# ---------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
# optional: write outputs
# library(writexl)

# ---- Automatically use current working directory ----
base_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(base_dir)
base_dir <- file.path(getwd(),"Data")

# Define Excel file paths
gp_sa_kegg        <- file.path(base_dir, "G_profiler_SA_kegg_wiki_reactome.xlsx")
gp_se_kegg        <- file.path(base_dir, "G_profiler_SE_keeg_wiki_reactome.xlsx")
gp_total_kegg     <- file.path(base_dir, "G_profiler_total_kegg_wiki_reactome.xlsx")
selected_20       <- file.path(base_dir, "20_pathways_escolhidos.xlsx")
gp_common_sa_se_d3 <- file.path(base_dir, "Common_SA_SE_D3_G_profiler_pathways.xlsx")
gp_common_sa_se_d7 <- file.path(base_dir, "Common_SA_SE_D7_G_profiler_pathways.xlsx")

# ---------------------------
# 1) Helper: safe sheet reader
# ---------------------------
read_sheet_safe <- function(path, sheet) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  message(sprintf("Reading sheet '%s' from %s", sheet, path))
  read_excel(path, sheet = sheet)
}

# ---------------------------
# 2) Load data (sheets used by original script)
# ---------------------------
SA_D3 <- read_sheet_safe(gp_sa_kegg, sheet = "SA Day 3")
SA_D7 <- read_sheet_safe(gp_sa_kegg, sheet = "SA Day 7")

SE_D3 <- read_sheet_safe(gp_se_kegg, sheet = "SE Day 3")
SE_D7 <- read_sheet_safe(gp_se_kegg, sheet = "SE Day 7")

Total_D3 <- read_sheet_safe(gp_total_kegg, sheet = "Total Dia 3")
Total_D7 <- read_sheet_safe(gp_total_kegg, sheet = "Total Dia 7")

Total_exc_D7_20 <- read_sheet_safe(selected_20, sheet = "Total_7_exclusive")
Total_common_20  <- read_sheet_safe(selected_20, sheet = "Total_Common")
SA_exc_D7_20     <- read_sheet_safe(selected_20, sheet = "SA_7_exclusive")

SE_SA_Comm20_3 <- read_sheet_safe(gp_common_sa_se_d3, sheet = "Mais relevantes")
SE_SA_Comm20_7 <- read_sheet_safe(gp_common_sa_se_d7, sheet = "Mais relevantes")

# ---------------------------
# 3) Numeric conversion helper
# ---------------------------
numericize <- function(df) {
  # Convert expected numeric-like columns to numeric safely
  df %>% mutate(
    negative_log10_of_adjusted_p_value = as.numeric(as.character(negative_log10_of_adjusted_p_value)),
    intersection_size = as.numeric(as.character(intersection_size)),
    adjusted_p_value = as.numeric(as.character(adjusted_p_value))
  )
}

SA_D3 <- numericize(SA_D3)
SA_D7 <- numericize(SA_D7)
SE_D3 <- numericize(SE_D3)
SE_D7 <- numericize(SE_D7)
Total_D3 <- numericize(Total_D3)
Total_D7 <- numericize(Total_D7)

# ---------------------------
# 4) Compute global color/size limits (optional)
#    This allows using the same color and size range when comparing D3 <-> D7
# ---------------------------
get_global_limits <- function(df3, df7) {
  list(
    pval = c(min(c(df3$adjusted_p_value, df7$adjusted_p_value), na.rm = TRUE),
             max(c(df3$adjusted_p_value, df7$adjusted_p_value), na.rm = TRUE)),
    size = c(min(c(df3$intersection_size, df7$intersection_size), na.rm = TRUE),
             max(c(df3$intersection_size, df7$intersection_size), na.rm = TRUE))
  )
}

limits_total <- get_global_limits(Total_D3, Total_D7)

# ---------------------------
# 5) Plot function
#    Creates a horizontal segment + point plot where x = -log10(adj p), size = # genes, color = FDR
# ---------------------------
gg_shinygo <- function(df, title = "", top_n = 20, color_limits = NULL, size_limits = NULL) {
  # If df is empty, return NULL with a message
  if (nrow(df) == 0) {
    message(sprintf("Data frame is empty — nothing to plot for: %s", title))
    return(NULL)
  }
  df_to_plot <- df %>%
    arrange(desc(negative_log10_of_adjusted_p_value)) %>%
    slice_head(n = min(top_n, nrow(.))) %>%
    mutate(term_name = factor(term_name, levels = rev(unique(term_name))))
  
  p <- ggplot(df_to_plot, aes(x = negative_log10_of_adjusted_p_value, y = term_name)) +
    geom_segment(aes(x = 0, xend = negative_log10_of_adjusted_p_value,
                     y = term_name, yend = term_name, color = adjusted_p_value), size = 1.1) +
    geom_point(aes(size = intersection_size, color = adjusted_p_value)) +
    labs(x = expression(-log[10](Adj.~P~value)), y = NULL, title = title) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_text(size = 11, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right")
  
  # Color scale (FDR)
  if (is.null(color_limits)) {
    p <- p + scale_color_gradient(low = "deeppink1", high = "deepskyblue", name = "FDR")
  } else {
    p <- p + scale_color_gradient(low = "deeppink1", high = "deepskyblue", limits = color_limits, name = "FDR")
  }
  
  # Size scale (# genes)
  if (is.null(size_limits)) {
    p <- p + scale_size_continuous(range = c(3, 8), name = "# Genes")
  } else {
    p <- p + scale_size_continuous(range = c(3, 8), limits = size_limits, name = "# Genes")
  }
  
  return(p)
}

# ---------------------------
# 6) Helper: get unique term names in a sheet
# ---------------------------
get_terms <- function(df) unique(df$term_name)

# ---------------------------
# 7) Create exclusive / common sets for SA and SE
# ---------------------------
exclusive_SA_D3 <- SA_D3 %>% filter(!term_name %in% get_terms(SE_D3))
exclusive_SA_D7 <- SA_D7 %>% filter(!term_name %in% get_terms(SE_D7))

exclusive_SE_D3 <- SE_D3 %>% filter(!term_name %in% get_terms(SA_D3))
exclusive_SE_D7 <- SE_D7 %>% filter(!term_name %in% get_terms(SA_D7))

common_SA_SE_terms_D3 <- intersect(get_terms(SA_D3), get_terms(SE_D3))
common_SA_SE_terms_D7 <- intersect(get_terms(SA_D7), get_terms(SE_D7))

common_SA_SE_D3 <- SA_D3 %>% filter(term_name %in% common_SA_SE_terms_D3)
common_SA_SE_D7 <- SA_D7 %>% filter(term_name %in% common_SA_SE_terms_D7)

# ---------------------------
# 8) Plot SA vs SE (examples)
# ---------------------------
p_exclusive_SA_D3 <- gg_shinygo(exclusive_SA_D3, "Exclusive SA — Day 3", top_n = nrow(exclusive_SA_D3))
p_exclusive_SA_D7 <- gg_shinygo(exclusive_SA_D7, "Exclusive SA — Day 7 (Top 20)", top_n = 20)
p_exclusive_SE_D3 <- gg_shinygo(exclusive_SE_D3, "Exclusive SE — Day 3", top_n = nrow(exclusive_SE_D3))
p_exclusive_SE_D7 <- gg_shinygo(exclusive_SE_D7, "Exclusive SE — Day 7 (Top 20)", top_n = 20)
p_common_SA_SE_D3 <- gg_shinygo(common_SA_SE_D3, "Common SA & SE — Day 3 (Top 20)", top_n = 20)
p_common_SA_SE_D7 <- gg_shinygo(common_SA_SE_D7, "Common SA & SE — Day 7 (Top 20)", top_n = 20)

# Print plots to device (RStudio will show them)
print(p_exclusive_SA_D3)
print(p_exclusive_SA_D7)
print(p_exclusive_SE_D3)
print(p_exclusive_SE_D7)
print(p_common_SA_SE_D3)
print(p_common_SA_SE_D7)

# ---------------------------
# 9) Total group with unified scales (examples)
# ---------------------------
exclusive_Total_D3 <- Total_D3 %>% filter(!term_name %in% get_terms(Total_D7))
exclusive_Total_D7 <- Total_D7 %>% filter(!term_name %in% get_terms(Total_D3))
common_Total_terms <- intersect(get_terms(Total_D3), get_terms(Total_D7))
common_Total <- Total_D3 %>% filter(term_name %in% common_Total_terms)

p_exclusive_Total_D3 <- gg_shinygo(exclusive_Total_D3, "Total — Exclusive Day 3", top_n = nrow(exclusive_Total_D3),
                                   color_limits = limits_total$pval, size_limits = limits_total$size)
p_exclusive_Total_D7 <- gg_shinygo(exclusive_Total_D7, "Total — Exclusive Day 7 (Top 20)", top_n = 20,
                                   color_limits = limits_total$pval, size_limits = limits_total$size)
p_common_Total <- gg_shinygo(common_Total, "Total — Common to Days 3 & 7 (Top 20)", top_n = 20,
                             color_limits = limits_total$pval, size_limits = limits_total$size)

print(p_exclusive_Total_D3)
print(p_exclusive_Total_D7)
print(p_common_Total)

# ---------------------------
# 10) Plot selected "top 20" sheets (already curated)
# ---------------------------
p_total_exc_D7_20 <- gg_shinygo(Total_exc_D7_20, "Selected — Total Exclusive Day 7", top_n = nrow(Total_exc_D7_20),
                                color_limits = limits_total$pval, size_limits = limits_total$size)
p_total_common_20 <- gg_shinygo(Total_common_20, "Selected — Total Common to Both Days", top_n = nrow(Total_common_20),
                                color_limits = limits_total$pval, size_limits = limits_total$size)
p_sa_exc_D7_20 <- gg_shinygo(SA_exc_D7_20, "Selected — Exclusive SA Day 7", top_n = nrow(SA_exc_D7_20),
                             color_limits = limits_total$pval, size_limits = limits_total$size)

print(p_total_exc_D7_20)
print(p_total_common_20)
print(p_sa_exc_D7_20)

# ---------------------------
# 11) Helper to plot all intersection terms if present
# ---------------------------
plot_all_common <- function(df, title) {
  if (nrow(df) == 0) {
    message(sprintf("No terms in intersection for: %s", title))
    return(NULL)
  }
  p <- gg_shinygo(df, title = title, top_n = nrow(df))
  print(p)
  return(p)
}

p_common_D3_all <- plot_all_common(common_SA_SE_D3, "Common SA & SE — Day 3 (All Terms)")
p_common_D7_all <- plot_all_common(common_SA_SE_D7, "Common SA & SE — Day 7 (All Terms)")

# ---------------------------
# 12) Additional curated plots (if present)
# ---------------------------
p_relevant_se_sa_3 <- gg_shinygo(SE_SA_Comm20_3, "Relevant Common SE & SA — Day 3",
                                color_limits = limits_total$pval, size_limits = limits_total$size)
p_relevant_se_sa_7 <- gg_shinygo(SE_SA_Comm20_7, "Relevant Common SE & SA — Day 7",
                                color_limits = limits_total$pval, size_limits = limits_total$size)

print(p_relevant_se_sa_3)
print(p_relevant_se_sa_7)

# Optionally: save plots to files, uncomment and edit filenames as needed
# ggsave("exclusive_SA_D3.png", p_exclusive_SA_D3, width = 10, height = 7, dpi = 300)
# ggsave("exclusive_Total_D7.png", p_exclusive_Total_D7, width = 10, height = 7, dpi = 300)
# write_xlsx(list(Exclusive_SA_D7 = exclusive_SA_D7), "exclusive_SA_D7_terms.xlsx")
