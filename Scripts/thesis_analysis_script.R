
#!/usr/bin/env Rscript
# thesis_analysis_script.R
# The script follows sequential sections: load data, preprocessing, filtering,
# FC calculation, ratio calculation, set construction (I/P/D), comparison with external PPI,
# and final outputs and plots.
#Note: On the Venn Digram's it's better if you c

# ---------------------------
# 0. Setup: packages and paths
# ---------------------------
library(readxl)
library(dplyr)
library(VennDiagram)
library(ggplot2)
library(writexl)
library(grid)
library(ggfortify)
library(factoextra)
library(tibble)

# --- Automatically use current working directory ---
base_dir <- getwd()  # Gets the folder where you are running the script

# --- Input file paths ---
input_excel_path <- file.path(base_dir, "Dissertação", "Proteome_Discoverer", "Results_amostras.xlsx")
biogrid_path     <- file.path(base_dir, "Dissertação", "Interactores da APP Biogrid", 
                              "BIOGRID-GENE-106848-4.4.246.DOWNLOADS", 
                              "BIOGRID-GENE-106848-4.4.246.tab3.txt")

# --- Output file paths ---
out_pre_filtering <- file.path(base_dir, "Dissertação", "Resultados", 
                               "Dados Proteinas antes da filtragem.xlsx")
out_groups_list   <- file.path(base_dir, "Dissertação", "Resultados", 
                               "Lista de Proteinas(I,P,D)_FC_1,5.xlsx")


# ---------------------------
# 1) Load data
# ---------------------------
Results_amostras <- read_excel(input_excel_path)

# ---------------------------
# 2) Preprocessing: replace NAs and compute group means
# ---------------------------
# Identify abundance columns matching the pattern "Abundances (Grouped): F"
abundance_cols <- grep("^Abundances \\(Grouped\\): F", names(Results_amostras), value = TRUE)

# Replace NA with 0 in abundance columns
Results_amostras[abundance_cols] <- lapply(Results_amostras[abundance_cols], function(x) ifelse(is.na(x), 0, x))

# Arithmetic mean helper (rounded)
arit_mean <- function(x) {
  round(mean(x), 2)
}

# Compute group-level mean abundance columns (map F columns to groups/days)
# WT group
Results_amostras$`Abundances (Grouped): APP_Wt_day_3` <- apply(Results_amostras[, c("Abundances (Grouped): F5", "Abundances (Grouped): F6")], 1, arit_mean)
Results_amostras$`Abundances (Grouped): APP_Wt_day_6` <- apply(Results_amostras[, c("Abundances (Grouped): F17", "Abundances (Grouped): F18")], 1, arit_mean)

# SA group
Results_amostras$`Abundances (Grouped): APP_SA_day_3` <- apply(Results_amostras[, c("Abundances (Grouped): F1", "Abundances (Grouped): F2")], 1, arit_mean)
Results_amostras$`Abundances (Grouped): APP_SA_day_6` <- apply(Results_amostras[, c("Abundances (Grouped): F11", "Abundances (Grouped): F12")], 1, arit_mean)

# SE group
Results_amostras$`Abundances (Grouped): APP_SE_day_3` <- apply(Results_amostras[, c("Abundances (Grouped): F3", "Abundances (Grouped): F4", "Abundances (Grouped): F13")], 1, arit_mean)
Results_amostras$`Abundances (Grouped): APP_SE_day_6` <- apply(Results_amostras[, c("Abundances (Grouped): F16", "Abundances (Grouped): F14", "Abundances (Grouped): F15")], 1, arit_mean)

# N1 group
Results_amostras$`Abundances (Grouped): APP_N1_day_3` <- apply(Results_amostras[, c("Abundances (Grouped): F7", "Abundances (Grouped): F9")], 1, arit_mean)
Results_amostras$`Abundances (Grouped): APP_N1_day_6` <- apply(Results_amostras[, c("Abundances (Grouped): F10", "Abundances (Grouped): F8")], 1, arit_mean)

# Write intermediate file (before filtering)
write_xlsx(Results_amostras, out_pre_filtering)

# ---------------------------
# 3) Initial filtering (contaminants, FDR, coverage, unique peptides)
# ---------------------------
Results_m <- Results_amostras %>%
  filter(Master == "Master Protein", Contaminant == FALSE) %>%
  filter(!Accession %in% c("P02533", "P13645", "Q5T749", "P35527",
                           "Q92764", "P04264", "Q7Z794", "P13646", "P19013", "Q9NSB4"))

results_no_contam <- Results_m %>%
  filter(`Protein FDR Confidence: Combined` == "High") %>%
  filter(!(`MW [kDa]` > 20 & `# Unique Peptides` == 1)) %>%
  filter(!(`Coverage [%]` < 20 & `# Unique Peptides` == 1))

results_filtered <- results_no_contam

# ---------------------------
# 4) Prepare per-comparison data frames (copies)
# ---------------------------
wt3 <- results_filtered
wt7 <- results_filtered
sa3 <- results_filtered
sa7 <- results_filtered
se3 <- results_filtered
se7 <- results_filtered

# ---------------------------
# 5) Adjusted fold-change function
# ---------------------------
ajusta_FC <- function(x, y) {
  x[is.na(x)] <- 0
  y[is.na(y)] <- 0
  fc <- mapply(function(a, b) {
    if (a == 0 && b == 0) {
      0
    } else if (b == 0 && a > 0) {
      100
    } else {
      a / b
    }
  }, x, y)
  round(fc, 2)
}

# ---------------------------
# 6) Compute FC vs N1 and filter FC > 1.5
# ---------------------------
# WT vs N1
wt3$FC <- ajusta_FC(wt3$`Abundances (Grouped): APP_Wt_day_3`, wt3$`Abundances (Grouped): APP_N1_day_3`)
wt7$FC <- ajusta_FC(wt7$`Abundances (Grouped): APP_Wt_day_6`, wt7$`Abundances (Grouped): APP_N1_day_6`)

# SA vs N1
sa3$FC <- ajusta_FC(sa3$`Abundances (Grouped): APP_SA_day_3`, sa3$`Abundances (Grouped): APP_N1_day_3`)
sa7$FC <- ajusta_FC(sa7$`Abundances (Grouped): APP_SA_day_6`, sa7$`Abundances (Grouped): APP_N1_day_6`)

# SE vs N1
se3$FC <- ajusta_FC(se3$`Abundances (Grouped): APP_SE_day_3`, se3$`Abundances (Grouped): APP_N1_day_3`)
se7$FC <- ajusta_FC(se7$`Abundances (Grouped): APP_SE_day_6`, se7$`Abundances (Grouped): APP_N1_day_6`)

# Optional aesthetic reordering function to move FC to column 5
reordena_FC <- function(df) {
  df[, c(1:4, ncol(df), 5:(ncol(df)-1))]
}
wt3 <- reordena_FC(wt3)
wt7 <- reordena_FC(wt7)
sa3 <- reordena_FC(sa3)
sa7 <- reordena_FC(sa7)
se3 <- reordena_FC(se3)
se7 <- reordena_FC(se7)

# Combine all datasets (not necessary for downstream but kept for record)
total_c_fc <- bind_rows(wt3, wt7, sa3, sa7, se3, se7)

# Filter by FC > 1.5 and remove duplicates
filtra_FC11 <- function(df) {
  df2 <- df %>% filter(FC > 1.5)
  unique(df2)
}
wt3 <- filtra_FC11(wt3)
wt7 <- filtra_FC11(wt7)
sa3 <- filtra_FC11(sa3)
sa7 <- filtra_FC11(sa7)
se3 <- filtra_FC11(se3)
se7 <- filtra_FC11(se7)

# ---------------------------
# 7) Compute SE/SA and SA/SE ratios in results_filtered
# ---------------------------
results_filtered <- results_filtered %>%
  mutate(
    ratio_SE_SA_day3 = round(
      ifelse(
        `Abundances (Grouped): APP_SA_day_3` == 0 & `Abundances (Grouped): APP_SE_day_3` == 0, 
        0,
        ifelse(
          `Abundances (Grouped): APP_SA_day_3` == 0, 
          100, 
          `Abundances (Grouped): APP_SE_day_3` / `Abundances (Grouped): APP_SA_day_3`
        )
      ), 
      2
    ),
    ratio_SA_SE_day3 = round(1 / pmax(ratio_SE_SA_day3, 1e-8), 2),
    ratio_SE_SA_day7 = round(
      ifelse(
        `Abundances (Grouped): APP_SA_day_6` == 0 & `Abundances (Grouped): APP_SE_day_6` == 0,
        0,
        ifelse(
          `Abundances (Grouped): APP_SA_day_6` == 0,
          100,
          `Abundances (Grouped): APP_SE_day_6` / `Abundances (Grouped): APP_SA_day_6`
        )
      ),
      2
    ),
    ratio_SA_SE_day7 = round(1 / pmax(ratio_SE_SA_day7, 1e-8), 2)
  )

# ---------------------------
# 8) Build accession lists per condition (day 3 and day 7)
# ---------------------------
wt3_acc <- wt3$Accession
sa3_acc <- sa3$Accession
se3_acc <- se3$Accession

wt7_acc <- wt7$Accession
sa7_acc <- sa7$Accession
se7_acc <- se7$Accession

# Day 3 groups
P_3 <- setdiff(se3_acc, sa3_acc)        # SE-only -> Phosphorylated (P)
D_3 <- setdiff(sa3_acc, se3_acc)        # SA-only -> Dephosphorylated (D)
common_3 <- intersect(sa3_acc, se3_acc) # SA ∩ SE
triple_3 <- Reduce(intersect, list(wt3_acc, sa3_acc, se3_acc))
I_3 <- unique(c(common_3, triple_3))    # Independent (I)

# Day 7 groups
P_7 <- setdiff(se7_acc, sa7_acc)
D_7 <- setdiff(sa7_acc, se7_acc)
common_7 <- intersect(sa7_acc, se7_acc)
triple_7 <- Reduce(intersect, list(wt7_acc, sa7_acc, se7_acc))
I_7 <- unique(c(common_7, triple_7))

# ---------------------------
# 9) Extract full data frames for each group (all columns)
# ---------------------------
I_3_all <- results_filtered %>% filter(Accession %in% I_3)
P_3_all <- results_filtered %>% filter(Accession %in% P_3)
D_3_all <- results_filtered %>% filter(Accession %in% D_3)

I_7_all <- results_filtered %>% filter(Accession %in% I_7)
P_7_all <- results_filtered %>% filter(Accession %in% P_7)
D_7_all <- results_filtered %>% filter(Accession %in% D_7)

# ---------------------------
# Helper: create modified SA/SE sets moving strongly enriched proteins out of intersection
# ---------------------------
cria_conjuntos_modificados <- function(sa_acc, se_acc, ratio_col_se_sa, ratio_col_sa_se) {
  base_inter <- intersect(sa_acc, se_acc)

  sa_de <- results_filtered %>%
    filter(Accession %in% base_inter, !!sym(ratio_col_sa_se) >= 2) %>%
    pull(Accession)
  se_de <- results_filtered %>%
    filter(Accession %in% base_inter, !!sym(ratio_col_se_sa) >= 2) %>%
    pull(Accession)

  inter_nao_de <- setdiff(base_inter, union(sa_de, se_de))

  sa_only_orig <- setdiff(sa_acc, se_acc)
  sa_only_novo  <- unique(c(sa_only_orig, sa_de))

  se_only_orig <- setdiff(se_acc, sa_acc)
  se_only_novo  <- unique(c(se_only_orig, se_de))

  list(
    SA_only    = sa_only_novo,
    SE_only    = se_only_novo,
    SA_SE_int  = inter_nao_de,
    SE_enrich  = se_de,
    SA_enrich  = sa_de
  )
}

# Apply for Day 3
mod3 <- cria_conjuntos_modificados(
  sa_acc = sa3_acc,
  se_acc = se3_acc,
  ratio_col_se_sa = "ratio_SE_SA_day3",
  ratio_col_sa_se = "ratio_SA_SE_day3"
)
SA3_only   <- mod3$SA_only
SE3_only   <- mod3$SE_only
INT3_naoDE <- mod3$SA_SE_int
SA3_enriched <- mod3$SA_enrich
SE3_enriched <- mod3$SE_enrich

# Prepare Venn counts for Day 3 (modified sets)
area1_3 <- length(wt3_acc)
area2_3 <- length(SA3_only) + length(INT3_naoDE)
area3_3 <- length(SE3_only) + length(INT3_naoDE)

n12_3 <- length(intersect(wt3_acc, SA3_only)) + length(intersect(wt3_acc, INT3_naoDE))
n13_3 <- length(intersect(wt3_acc, SE3_only)) + length(intersect(wt3_acc, INT3_naoDE))
n23_3 <- length(INT3_naoDE)
n123_3 <- length(intersect(wt3_acc, INT3_naoDE))

grid::grid.newpage() #Cleanes the plots that were in the environmet before because the Venn Digrams
# tend to be sobreposed on the other plots, so in order to see each plot is better to run each one separately
# Draw triple Venn (Day 3)
venn3_dia3 <- draw.triple.venn(
  area1 = area1_3,
  area2 = area2_3,
  area3 = area3_3,
  n12   = n12_3,
  n13   = n13_3,
  n23   = n23_3,
  n123  = n123_3,
  category = c("SE", "SA", "WT"),
  fill     = c("orange", "#0090B9", "springgreen3"),
  alpha    = 0.75,
  cex      = 2,
  cat.cex  = 1.3
)

grid::grid.newpage()
# Apply for Day 7
mod7 <- cria_conjuntos_modificados(
  sa_acc = sa7_acc,
  se_acc = se7_acc,
  ratio_col_se_sa = "ratio_SE_SA_day7",
  ratio_col_sa_se = "ratio_SA_SE_day7"
)
SA7_only <- mod7$SA_only
SE7_only <- mod7$SE_only
INT7_naoDE <- mod7$SA_SE_int
SA7_enriched <- mod7$SA_enrich
SE7_enriched <- mod7$SE_enrich

# Prepare Venn counts for Day 7 (modified sets)
area1_7 <- length(wt7_acc)
area2_7 <- length(SA7_only) + length(INT7_naoDE)
area3_7 <- length(SE7_only) + length(INT7_naoDE)

n12_7 <- length(intersect(wt7_acc, SA7_only)) + length(intersect(wt7_acc, INT7_naoDE))
n13_7 <- length(intersect(wt7_acc, SE7_only)) + length(intersect(wt7_acc, INT7_naoDE))
n23_7 <- length(INT7_naoDE)
n123_7 <- length(intersect(wt7_acc, INT7_naoDE))

# Draw triple Venn (Day 7)
venn3_dia7 <- draw.triple.venn(
  area1 = area1_7,
  area2 = area2_7,
  area3 = area3_7,
  n12   = n12_7,
  n13   = n13_7,
  n23   = n23_7,
  n123  = n123_7,
  category = c("SE", "SA", "WT"),
  fill     = c("orange", "#0090B9", "springgreen3"),
  alpha    = 0.75,
  cex      = 2,
  cat.cex  = 1.3
)

# ---------------------------
# 11) Optional: combine Day 3 and Day 7 groups into aggregated lists
# ---------------------------
I_all <- bind_rows(I_3_all, I_7_all) %>% distinct()
P_all <- bind_rows(P_3_all, P_7_all) %>% distinct()
D_all <- bind_rows(D_3_all, D_7_all) %>% distinct()
Total <- bind_rows(I_all, P_all, D_all) %>% distinct()

# Create final lists of accessions per group/day and combined interactors
interactors_day3 <- unique(c(wt3$Accession, sa3$Accession, se3$Accession))
interactors_day7 <- unique(c(wt7$Accession, sa7$Accession, se7$Accession))
Interactors_totals <- unique(c(interactors_day3, interactors_day7))

# Prepare df_final: columns for I_3, I_7, D_3, D_7, P_3, P_7, Total_Dia3, Total_Dia7
max_len <- max(length(I_3), length(P_3), length(D_3), length(I_7), length(P_7), length(D_7),
               length(interactors_day3), length(interactors_day7))

fill_na <- function(x, max_len) {
  c(x, rep(NA, max_len - length(x)))
}

df_final <- data.frame(
  I_3 = fill_na(I_3, max_len),
  I_7 = fill_na(I_7, max_len),
  D_3 = fill_na(D_3, max_len),
  D_7 = fill_na(D_7, max_len),
  P_3 = fill_na(P_3, max_len),
  P_7 = fill_na(P_7, max_len),
  Total_Dia3 = fill_na(interactors_day3, max_len),
  Total_Dia7 = fill_na(interactors_day7, max_len)
)

# Write final Excel with grouped lists
write_xlsx(df_final, out_groups_list)

# ---------------------------
# Interactor comparisons with external PPI resources (BioGRID and APID)
# ---------------------------
BIOGRID_APP <- read.delim(biogrid_path, quote = "")

interactors_common <- Interactors_totals[
  Interactors_totals %in% BIOGRID_APP$SWISS.PROT.Accessions.Interactor.A |
  Interactors_totals %in% BIOGRID_APP$SWISS.PROT.Accessions.Interactor.B
]

interactors_common_3 <- interactors_day3[
  interactors_day3 %in% BIOGRID_APP$SWISS.PROT.Accessions.Interactor.A |
  interactors_day3 %in% BIOGRID_APP$SWISS.PROT.Accessions.Interactor.B
]

interactors_common_7 <- interactors_day7[
  interactors_day7 %in% BIOGRID_APP$SWISS.PROT.Accessions.Interactor.A |
  interactors_day7 %in% BIOGRID_APP$SWISS.PROT.Accessions.Interactor.B
]

interactors_APP_total <- results_filtered %>% filter(Accession %in% interactors_common)
interactors_APP_3     <- results_filtered %>% filter(Accession %in% interactors_common_3)
interactors_APP_7     <- results_filtered %>% filter(Accession %in% interactors_common_7)

# ---------------------------
# 12) Build compact dataset summarizing per-condition accessions
# ---------------------------
interactors_day3 <- unique(c(wt3$Accession, sa3$Accession, se3$Accession))
interactors_day7 <- unique(c(wt7$Accession, sa7$Accession, se7$Accession))

max_len2 <- max(
  length(se3_acc), length(se7_acc),
  length(sa3_acc), length(sa7_acc),
  length(wt3_acc), length(wt7_acc),
  length(interactors_day3), length(interactors_day7)
)

fill_na2 <- function(x, max_len2) {
  c(x, rep(NA, max_len2 - length(x)))
}

df_groups <- data.frame(
  SE_3       = fill_na2(se3_acc, max_len2),
  SE_7       = fill_na2(se7_acc, max_len2),
  SA_3       = fill_na2(sa3_acc, max_len2),
  SA_7       = fill_na2(sa7_acc, max_len2),
  WT_3       = fill_na2(wt3_acc, max_len2),
  WT_7       = fill_na2(wt7_acc, max_len2),
  Total_Dia3 = fill_na2(interactors_day3, max_len2),
  Total_Dia7 = fill_na2(interactors_day7, max_len2)
)

# Display first rows for quick inspection (prints to console)
print(head(df_groups))

# ---------------------------
# 13) Plots: barplots and Venn diagrams
# ---------------------------
# Barplots: counts per condition by day
counts_df <- data.frame(
  Condition = c("Wt Day 3", "SE Day 3", "SA Day 3", "Wt Day 7", "SE Day 7", "SA Day 7"),
  Count = c(nrow(wt3), nrow(se3), nrow(sa3), nrow(wt7), nrow(se7), nrow(sa7)),
  Group = c("Wt", "SE", "SA", "Wt", "SE", "SA"),
  Day = c("Day 3", "Day 3", "Day 3", "Day 7", "Day 7", "Day 7")
)
counts_df$Condition <- factor(counts_df$Condition, levels = c("Wt Day 3", "SE Day 3", "SA Day 3", "Wt Day 7", "SE Day 7", "SA Day 7"))

cores_grupos <- c("Wt" = "springgreen3", "SE" = "orange", "SA" = "#0090B9")

# Plot: Day 3
ggplot(counts_df %>% filter(Day == "Day 3"), aes(x = Condition, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 300)) +
  scale_fill_manual(values = cores_grupos) +
  theme_classic() +
  labs(y = "Number of proteins", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))

# Plot: Day 7
ggplot(counts_df %>% filter(Day == "Day 7"), aes(x = Condition, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_cartesian(ylim = c(0, 300)) +
  scale_fill_manual(values = cores_grupos) +
  theme_classic() +
  labs(y = "Number of proteins", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))

# Pairwise Venn for SA Day3 vs Day7
grid::grid.newpage()
venn_D_SA <- draw.pairwise.venn(
  area1 = length(sa3_acc),
  area2 = length(sa7_acc),
  cross.area = length(intersect(sa3_acc, sa7_acc)),
  category = c("Day 3", "Day 7"),
  fill = c("lightseagreen", "yellow"),
  alpha = 0.75,
  cat.cex = 1.5,
  cex = 2,
  inverted = TRUE,
  scaled = TRUE,
  cat.pos = c(180, 180),
  cat.dist = c(0.05, 0.05)
)

# Pairwise Venn for SE Day3 vs Day7
grid::grid.newpage()
venn_D_SE <- draw.pairwise.venn(
  area1 = length(se3_acc),
  area2 = length(se7_acc),
  cross.area = length(intersect(se3_acc, se7_acc)),
  category = c("Day 3", "Day 7"),
  fill = c("lightseagreen", "yellow"),
  alpha = 0.75,
  cat.cex = 1.5,
  cex = 2,
  inverted = TRUE,
  scaled = TRUE,
  cat.pos = c(180, 180),
  cat.dist = c(0.05, 0.05)
)

# Classic triple Venns using raw per-condition sets - the ones used in the thesis
grid::grid.newpage()
venn3_raw_day3 <- draw.triple.venn(
  area1 = length(wt3$Accession),
  area2 = length(sa3$Accession),
  area3 = length(se3$Accession),
  n12 = length(intersect(wt3$Accession, sa3$Accession)),
  n13 = length(intersect(wt3$Accession, se3$Accession)),
  n23 = length(intersect(sa3$Accession, se3$Accession)),
  n123 = length(Reduce(intersect, list(wt3$Accession, sa3$Accession, se3$Accession))),
  category = c("WT", "SA", "SE"),
  fill = c("orange", "#0090B9", "springgreen3"),
  alpha = 0.75,
  cat.cex = 1.5,
  cex = 2
)

grid::grid.newpage()
venn3_raw_day7 <- draw.triple.venn(
  area1 = length(wt7$Accession),
  area2 = length(sa7$Accession),
  area3 = length(se7$Accession),
  n12 = length(intersect(wt7$Accession, sa7$Accession)),
  n13 = length(intersect(wt7$Accession, se7$Accession)),
  n23 = length(intersect(sa7$Accession, se7$Accession)),
  n123 = length(Reduce(intersect, list(wt7$Accession, sa7$Accession, se7$Accession))),
  category = c("WT", "SA", "SE"),
  fill = c("orange", "#0090B9", "springgreen3"),
  alpha = 0.75,
  cat.cex = 1.5,
  cex = 2
)

################################################################################
# 🧩 PCA Analysis of Filtered Proteins (FC > 1.5)
################################################################################
#Creating a dataset with all the FC calculated
results_fc_all <- results_filtered %>%
  mutate(
    FC_WT_day3 = ajusta_FC(`Abundances (Grouped): APP_Wt_day_3`, `Abundances (Grouped): APP_N1_day_3`),
    FC_WT_day7 = ajusta_FC(`Abundances (Grouped): APP_Wt_day_6`, `Abundances (Grouped): APP_N1_day_6`),
    FC_SA_day3 = ajusta_FC(`Abundances (Grouped): APP_SA_day_3`, `Abundances (Grouped): APP_N1_day_3`),
    FC_SA_day7 = ajusta_FC(`Abundances (Grouped): APP_SA_day_6`, `Abundances (Grouped): APP_N1_day_6`),
    FC_SE_day3 = ajusta_FC(`Abundances (Grouped): APP_SE_day_3`, `Abundances (Grouped): APP_N1_day_3`),
    FC_SE_day7 = ajusta_FC(`Abundances (Grouped): APP_SE_day_6`, `Abundances (Grouped): APP_N1_day_6`)
  )

# ---------------------------
# 1) Prepare Day 3 and Day 7 datasets
# ---------------------------

# Unique filtered proteins per day
proteins_d3 <- unique(c(sa3$Accession, se3$Accession, wt3$Accession))
proteins_d7 <- unique(c(sa7$Accession, se7$Accession, wt7$Accession))

# Select abundance columns
samples_d3 <- results_filtered %>%
  filter(Accession %in% proteins_d3) %>%
  select(
    `Abundances (Grouped): F1`, `Abundances (Grouped): F2`,                # SA
    `Abundances (Grouped): F3`, `Abundances (Grouped): F4`, `Abundances (Grouped): F13`,  # SE
    `Abundances (Grouped): F5`, `Abundances (Grouped): F6`                 # WT
  )

samples_d7 <- results_filtered %>%
  filter(Accession %in% proteins_d7) %>%
  select(
    `Abundances (Grouped): F11`, `Abundances (Grouped): F12`,              # SA
    `Abundances (Grouped): F14`, `Abundances (Grouped): F15`, `Abundances (Grouped): F16`, # SE
    `Abundances (Grouped): F17`, `Abundances (Grouped): F18`               # WT
  )

# ---------------------------
# 2) Preprocessing function
# ---------------------------
prep_pca_data <- function(df, accessions) {
  df <- log2(df + 1)
  df <- as.data.frame(t(df))
  colnames(df) <- accessions
  df <- df[, apply(df, 2, var) != 0]  # Remove proteins with zero variance
  return(df)
}

# Apply preprocessing
pca_data_d3 <- prep_pca_data(samples_d3, proteins_d3)
pca_data_d7 <- prep_pca_data(samples_d7, proteins_d7)

# ---------------------------
# 3) Define groups and run PCA
# ---------------------------
groups_d3 <- c("SA", "SA", "SE", "SE", "SE", "WT", "WT")
groups_d7 <- c("SA", "SA", "SE", "SE", "SE", "WT", "WT")

pca_d3 <- prcomp(pca_data_d3, center = TRUE, scale. = TRUE)
pca_d7 <- prcomp(pca_data_d7, center = TRUE, scale. = TRUE)

# ---------------------------
# 4) Plot PCA (Day 3 and Day 7)
# ---------------------------
autoplot(pca_d3, data = data.frame(Group = groups_d3),
         colour = 'Group', size = 4) +
  theme_classic() +
  ggtitle("PCA (Day 3)") +
  scale_color_manual(values = c("WT" = "springgreen3", "SE" = "#0090B9", "SA" = "orange"))

autoplot(pca_d7, data = data.frame(Group = groups_d7),
         colour = 'Group', size = 4) +
  theme_classic() +
  ggtitle("PCA (Day 7)") +
  scale_color_manual(values = c("WT" = "springgreen3", "SE" = "#0090B9", "SA" = "orange"))


################################################################################
# 🔎 Top contributing proteins to PC1 and PC2
################################################################################

get_top_pca_proteins <- function(pca_obj, results_fc_all, top_n = 10, day = 3) {
  loadings <- pca_obj$rotation
  
  # Select top N by absolute loading values
  top_pc1 <- sort(abs(loadings[, 1]), decreasing = TRUE)[1:top_n]
  top_pc2 <- sort(abs(loadings[, 2]), decreasing = TRUE)[1:top_n]
  
  # Map Accession → Gene Symbol
  accession_to_gene <- results_fc_all %>%
    select(Accession, `Gene Symbol`) %>%
    distinct() %>%
    deframe()
  
  df_pc1 <- data.frame(
    Gene = accession_to_gene[names(top_pc1)],
    Accession = names(top_pc1),
    Loading = top_pc1,
    PC = "PC1"
  )
  
  df_pc2 <- data.frame(
    Gene = accession_to_gene[names(top_pc2)],
    Accession = names(top_pc2),
    Loading = top_pc2,
    PC = "PC2"
  )
  
  df_out <- bind_rows(df_pc1, df_pc2)
  
  # Identify abundance columns (fallback-safe)
  get_existing_col <- function(base, day, df_names) {
    candidates <- c(paste0(base, day),
                    paste0(base, ifelse(day == 3, 2, ifelse(day == 7, 6, day - 1))))
    found <- intersect(candidates, df_names)
    if (length(found) == 0) stop("No abundance column found for base: ", base, " day: ", day)
    return(found[1])
  }
  
  df_names <- names(results_fc_all)
  se_col <- get_existing_col("Abundances (Grouped): APP_SE_day_", day, df_names)
  sa_col <- get_existing_col("Abundances (Grouped): APP_SA_day_", day, df_names)
  wt_col <- get_existing_col("Abundances (Grouped): APP_Wt_day_", day, df_names)
  n1_col <- get_existing_col("Abundances (Grouped): APP_N1_day_", day, df_names)
  
  # Join relative abundance values
  df_out <- df_out %>%
    left_join(
      results_fc_all %>% select(Accession, all_of(c(se_col, sa_col, wt_col, n1_col))),
      by = "Accession"
    ) %>%
    mutate(
      WT_rel = pmax(.data[[wt_col]] - .data[[n1_col]], 0),
      SE_rel = pmax(.data[[se_col]] - .data[[n1_col]], 0),
      SA_rel = pmax(.data[[sa_col]] - .data[[n1_col]], 0)
    ) %>%
    select(Gene, Accession, PC, Loading, WT_rel, SE_rel, SA_rel)
  
  return(df_out)
}

# Extract top proteins
top_day3 <- get_top_pca_proteins(pca_d3, results_fc_all, top_n = 10, day = 3)
top_day7 <- get_top_pca_proteins(pca_d7, results_fc_all, top_n = 10, day = 7)

# Export to Excel
write_xlsx(list("Top_Day3" = top_day3,
                "Top_Day7" = top_day7),
           "Top_proteins_PCA_abun_rel.xlsx")


################################################################################
# 🧩 Combined PCA (Day 3 + Day 7)
################################################################################

# Combine all proteins across days
proteins_all <- unique(c(sa3$Accession, se3$Accession, wt3$Accession,
                         sa7$Accession, se7$Accession, wt7$Accession))

# Select abundance columns
samples_all <- results_filtered %>%
  filter(Accession %in% proteins_all) %>%
  select(
    # Day 3
    `Abundances (Grouped): F1`, `Abundances (Grouped): F2`,                # SA
    `Abundances (Grouped): F3`, `Abundances (Grouped): F4`, `Abundances (Grouped): F13`,  # SE
    `Abundances (Grouped): F5`, `Abundances (Grouped): F6`,                # WT
    # Day 7
    `Abundances (Grouped): F11`, `Abundances (Grouped): F12`,              # SA
    `Abundances (Grouped): F14`, `Abundances (Grouped): F15`, `Abundances (Grouped): F16`, # SE
    `Abundances (Grouped): F17`, `Abundances (Grouped): F18`               # WT
  )

# Preprocess
pca_all_data <- prep_pca_data(samples_all, proteins_all)

# Group and day vectors
group_all <- c("SA", "SA", "SE", "SE", "SE", "WT", "WT",
               "SA", "SA", "SE", "SE", "SE", "WT", "WT")
day_all <- c(rep("Day 3", 7), rep("Day 7", 7))

# PCA and plot
pca_all <- prcomp(pca_all_data, center = TRUE, scale. = TRUE)

autoplot(pca_all, data = data.frame(Group = group_all, Day = day_all),
         colour = 'Group', shape = 'Day', size = 4) +
  theme_classic() +
  ggtitle("PCA (Day 3 + Day 7)") +
  scale_color_manual(values = c("WT" = "springgreen3", "SE" = "#0090B9", "SA" = "orange")) +
  scale_shape_manual(values = c("Day 3" = 16, "Day 7" = 17))  # circles = day3, triangles = day7


# End of script
