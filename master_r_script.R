# Master R Script for MSc Thesis: TE Analysis in H.contortus
# 
# Below are both the scripts used for to obtain the results in the thesis for reproducibility:
# 1. EarlGrey Annotation Analysis (TE identification, genome composition, overlaps with features, and comparative plots across species)
# 2. PoPoolationTE2 Analysis (TE insertion frequencies, PCA of TE family presence, Upset plots, logistic regression, Manhattan plots)
#
# All functions, plotting themes, and processing steps from both scripts are included here.

# SCRIPT 1 - Earlgrey Annotation Analysis
# -------------------------------
# 1. Load Required Libraries
# -------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# -------------------------------
# 2. Define Plotting Themes & Colors
# -------------------------------
# Species-specific colors
species_colors <- c(
  "H. contortus" = "#1f78b4",
  "T. circumcincta" = "#33a02c",
  "C. elegans" = "#e31a1c"
)

# Clean minimal theme
clean_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Class colors for TE plots
class_colors <- c(
  "DNA"   = "#1f77b4",
  "LINE"  = "#ff7f0e",
  "LTR"   = "#2ca02c",
  "PLE"   = "#d62728",
  "RC"    = "#9467bd",
  "SINE"  = "#e377c2"
)

# -------------------------------
# 3. File Paths 
# -------------------------------
hc_repeat_gff_file <- "earlgrey_new_output/Hcontortus_EarlGrey/Hcontortus_summaryFiles/Hcontortus.filteredRepeats.gff"
tc_repeat_gff_file <- "Tcircumcincta_EarlGrey/Tcircumcincta_summaryFiles/Tcircumcincta.filteredRepeats.gff"
ce_repeat_gff_file <- "Celegans_summaryFiles/Celegans.filteredRepeats.gff"

hc_repeat_bed_file <- "earlgrey_new_output/Hcontortus_EarlGrey/Hcontortus_summaryFiles/Hcontortus.filteredRepeats.bed"

hc_family_counts_file <- "earlgrey_new_output/Hcontortus_EarlGrey/Hcontortus_summaryFiles/Hcontortus.familyLevelCount.txt"
tc_family_counts_file <- "Tcircumcincta_EarlGrey/Tcircumcincta_summaryFiles/Tcircumcincta.familyLevelCount.txt"
ce_family_counts_file <- "Celegans_summaryFiles/Celegans.familyLevelCount.txt"

hc_annot_file <- "earlgrey_new_output/Hcontortus_EarlGrey/Hcontortus_summaryFiles/original_data/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3"
tc_annot_file <- "Tcircumcincta_EarlGrey/Tcircumcincta_summaryFiles/original_data/teladorsagia_circumcincta_tci2_wsi3.0.annotation.gff3"
ce_annot_file <- "Celegans_summaryFiles/ce_filtered.gff3"

# -------------------------------
# 4. Functions
# -------------------------------

# 4.1 Extract repeat info from GFF
extract_repeat_info <- function(gff_df) {
  colnames(gff_df) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  gff_df <- gff_df %>%
    mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      score = as.numeric(score)
    ) %>%
    separate(type, into = c("class", "subclass"), sep = "/", fill = "right") %>%
    mutate(
      length = end - start + 1,
      family = str_extract(attributes, "ID=([^;]+)") %>% str_remove("ID="),
      kimura = str_extract(attributes, "KIMURA80=\\d+\\.\\d+") %>% str_remove("KIMURA80=") %>% as.numeric()
    )
  return(gff_df)
}

# 4.2 Get copy number per class
get_copy_number <- function(gff_df, species_name) {
  gff_df %>%
    count(class, name = "copy_number") %>%
    mutate(species = species_name, te_class = class) %>%
    select(-class)
}

# 4.3 Get family counts
get_family_counts <- function(family_df, species_name) {
  family_df %>%
    mutate(species = species_name, te_class = `TE Family`, family_count = `Copy Number`) %>%
    select(-`TE Family`, -`Copy Number`)
}

# 4.4 Family richness
get_family_richness <- function(gff_df) {
  gff_df %>%
    filter(!is.na(family)) %>%
    distinct(class, family) %>%
    group_by(class) %>%
    summarise(distinct_families = n(), .groups = "drop")
}

# 4.5 Estimate genome size
estimate_genome_size <- function(gff_df) {
  max(gff_df$end, na.rm = TRUE)
}

# 4.6 Filter TE classes
filter_te <- function(gff_df) {
  non_te_classes <- c("Simple_repeat", "Unknown", "Satellite", "Low_complexity")
  gff_df %>% filter(!class %in% non_te_classes)
}

# 4.7 Convert DF to GRanges
df_to_gr <- function(df) {
  GRanges(seqnames = df$seqid,
          ranges = IRanges(start = df$start, end = df$end),
          strand = df$strand,
          class = df$class,
          family = df$family)
}

# 4.8 Load annotation features
load_annotation_features <- function(path) {
  g <- import(path)
  list(
    cds = g[g$type == "CDS"],
    utr = g[g$type %in% c("five_prime_UTR", "three_prime_UTR")],
    intron = {
      intr <- g[g$type == "intron"]
      if (length(intr) == 0) {
        txdb <- split(g[g$type == "exon"], g$Parent)
        intronsByTranscript(txdb, use.names = TRUE) |> unlist()
      } else intr
    }
  )
}

# 4.9 Count overlaps between TE and features
count_overlaps <- function(te_gr, feat_gr, feat_name, species) {
  if (length(feat_gr) == 0) return(tibble())
  hits <- findOverlaps(te_gr, feat_gr)
  tibble(
    species = species,
    feature = feat_name,
    te_class = mcols(te_gr)$class[queryHits(hits)]
  ) %>%
    count(species, feature, te_class, name = "count")
}

# 4.10 Load BED coverage and merge with TE data
load_cov_and_merge <- function(te_df, bed_path) {
  cov_df <- read_tsv(bed_path, col_names = FALSE, show_col_types = FALSE) %>%
    rename(seqid = X1, start = X2, end = X3, coverage = X4)
  left_join(te_df, cov_df, by = c("seqid", "start", "end"))
}

# 4.11 Function to calculate TE overlap with genomic features
calc_te_feature_overlap <- function(te_gr, features_gr, genome_size) {
  overlaps <- findOverlaps(te_gr, features_gr)
  te_hits <- te_gr[queryHits(overlaps)]
  feat_hits <- features_gr[subjectHits(overlaps)]
  
  overlap_widths <- width(pintersect(ranges(te_hits), ranges(feat_hits)))
  
  tibble(
    te_class = mcols(te_hits)$class,
    feature_type = mcols(feat_hits)$feature_type,
    bp_overlap = overlap_widths
  ) %>%
    group_by(te_class, feature_type) %>%
    summarise(bp_covered = sum(bp_overlap), .groups = "drop") %>%
    mutate(percent_of_genome = 100 * bp_covered / genome_size)
}

# 4.12 Calculate TE vs non-TE genome coverage
calc_te_vs_non_te <- function(gff_df, genome_size) {
  gr <- GRanges(seqnames = gff_df$seqid, ranges = IRanges(start = gff_df$start, end = gff_df$end))
  gr_merged <- reduce(gr)
  te_bp <- sum(width(gr_merged))
  non_te_bp <- genome_size - te_bp
  tibble(class = c("TE", "Non-repetitive genome"),
         bp = c(te_bp, non_te_bp),
         percent = 100 * c(te_bp, non_te_bp) / genome_size)
}

# 4.13 Calculate TE coverage per class
calc_te_per_class <- function(gff_df, genome_size) {
  gff_gr <- GRanges(seqnames = gff_df$seqid, ranges = IRanges(gff_df$start, gff_df$end), class = gff_df$class)
  per_class <- sapply(split(gff_gr, mcols(gff_gr)$class), function(gr) sum(width(reduce(gr))))
  enframe(per_class, name = "class", value = "bp") %>% mutate(percent = 100 * bp / genome_size)
}

# -------------------------------
# 5. Load Data
# -------------------------------
hc_gff_raw <- read_tsv(hc_repeat_gff_file, col_names = FALSE, comment = "#", show_col_types = FALSE)
tc_gff_raw <- read_tsv(tc_repeat_gff_file, col_names = FALSE, comment = "#", show_col_types = FALSE)
ce_gff_raw <- read_tsv(ce_repeat_gff_file, col_names = FALSE, comment = "#", show_col_types = FALSE)

hc_gff <- extract_repeat_info(hc_gff_raw)
tc_gff <- extract_repeat_info(tc_gff_raw)
ce_gff <- extract_repeat_info(ce_gff_raw)

hc_gff <- filter_te(hc_gff)
tc_gff <- filter_te(tc_gff)
ce_gff <- filter_te(ce_gff)

# -------------------------------
# 6. Compute Summaries
# -------------------------------
hc_genome_size <- estimate_genome_size(hc_gff)
tc_genome_size <- estimate_genome_size(tc_gff)
ce_genome_size <- estimate_genome_size(ce_gff)

hc_copy <- get_copy_number(hc_gff, "H. contortus")
tc_copy <- get_copy_number(tc_gff, "T. circumcincta")
ce_copy <- get_copy_number(ce_gff, "C. elegans")

hc_family_raw <- read_tsv(hc_family_counts_file, show_col_types = FALSE)
tc_family_raw <- read_tsv(tc_family_counts_file, show_col_types = FALSE)
ce_family_raw <- read_tsv(ce_family_counts_file, show_col_types = FALSE)

hc_family <- get_family_counts(hc_family_raw, "H. contortus")
tc_family <- get_family_counts(tc_family_raw, "T. circumcincta")
ce_family <- get_family_counts(ce_family_raw, "C. elegans")

hc_richness <- get_family_richness(hc_gff) %>% mutate(species = "H. contortus")
tc_richness <- get_family_richness(tc_gff) %>% mutate(species = "T. circumcincta")
ce_richness <- get_family_richness(ce_gff) %>% mutate(species = "C. elegans")

# Combine summaries
combined_summary <- bind_rows(hc_copy, tc_copy, ce_copy)
combined_family <- bind_rows(hc_family, tc_family, ce_family)
combined_richness <- bind_rows(hc_richness, tc_richness, ce_richness)

# -------------------------------
# 7. Save Filtered GFFs
# -------------------------------
write_tsv(hc_gff, "Hcontortus_filtered_TE.gff", col_names = FALSE)
write_tsv(tc_gff, "Tcircumcincta_filtered_TE.gff", col_names = FALSE)
write_tsv(ce_gff, "Celegans_filtered_TE.gff", col_names = FALSE)

# Optional: rebuild GFF with attributes
hc_gff_fixed <- hc_gff %>%
  mutate(
    tstart = str_extract(attributes, "TSTART=[^;]+"),
    tend = str_extract(attributes, "TEND=[^;]+"),
    orig_id = str_extract(attributes, "ID=([^;]+)") %>% str_remove("ID="),
    shortte = str_extract(attributes, "SHORTTE=[^;]+"),
    id = class,
    family = orig_id,
    kimura80 = ifelse(!is.na(kimura), paste0("KIMURA80=", kimura), NA),
    attributes = paste(tstart, tend, paste0("ID=", id, "_", family), shortte, kimura80, sep = ";") %>% str_remove(";NA$")
  ) %>%
  select(seqid, source, class, start, end, score, strand, phase, attributes) %>%
  arrange(seqid, start)

write_tsv(hc_gff_fixed, "Hcontortus_te_eg.gff", col_names = FALSE)

# -------------------------------
# 8. TE Overlaps with Features
# -------------------------------

# Load genomic features
hc_feats <- load_annotation_features(hc_annot_file)
hc_te_gr <- df_to_gr(hc_gff)

# Combine features into one GRanges with feature_type
hc_feats_combined <- c(hc_feats$cds, hc_feats$utr, hc_feats$intron)
mcols(hc_feats_combined)$feature_type <- rep(
  c("CDS","UTR","Intron"),
  times = c(length(hc_feats$cds), length(hc_feats$utr), length(hc_feats$intron))
)

# -------------------------------
# 9. TE Feature Overlap Analysis
# -------------------------------

coverage_df <- calc_te_feature_overlap(hc_te_gr, hc_feats, hc_genome_size)

# Zero-fill missing TE class x feature combinations
all_te_classes <- unique(mcols(hc_te_gr)$class)
all_features <- unique(mcols(hc_feats)$feature_type)

complete_df <- expand.grid(te_class = all_te_classes, feature_type = all_features) %>%
  as_tibble() %>%
  left_join(coverage_df, by = c("te_class", "feature_type")) %>%
  mutate(
    bp_covered = replace_na(bp_covered, 0),
    percent_of_genome = replace_na(percent_of_genome, 0)
  )

# Plot TE coverage by feature and class
covp <- complete_df %>%
  filter(bp_covered > 0) %>%  # optional: only show non-zero coverage
  ggplot(aes(x = feature_type, y = percent_of_genome, fill = te_class)) +
  geom_col(position = position_stack(reverse = TRUE), color = "black") +
  geom_text(data = . %>% filter(percent_of_genome > 0.1),
            aes(label = paste0(round(percent_of_genome, 2), "%")),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3, color = "black") +
  scale_fill_brewer(palette = "Set2", name = "TE Class") +
  labs(title = "Percent of TE Classes Overlapping Genomic Features\nin Haemonchus contortus",
       x = "Genomic Feature",
       y = "TE Coverage (% of genome)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  clean_theme()

ggsave("hc_te_coverage_by_feature_and_class.png", covp, width = 9, height = 6, dpi = 300)

# -------------------------------
# 10. Plotting the data for visualization
# -------------------------------

# Genome Composition Pie Chart (TE classes + non-repetitive)

hc_te_non_te <- calc_te_vs_non_te(hc_gff, genome_size)
hc_per_class <- calc_te_per_class(hc_gff, genome_size)

hc_pie_data <- bind_rows(
  hc_per_class %>% select(class, bp, percent),
  hc_te_non_te %>% filter(class == "Non-repetitive genome") %>% select(class, bp, percent)
) %>%
  arrange(desc(class)) %>%
  mutate(ypos = cumsum(percent) - 0.5 * percent)

pie_colors <- c(
  "DNA" = "#6baed6", "LINE"="#fdae6b", "LTR"="#74c476", "PLE"="#fb6a4a", "RC"="#9e9ac8",
  "Non-repetitive genome"="grey70"
)


ggplot(hc_pie_data, aes(x = "", y = percent, fill = class)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = pie_colors) +
  geom_label_repel(
    aes(y = ypos, label = paste0(round(percent,2), "%")),
    size = 2, nudge_x = 1, show.legend = FALSE,
    segment.color = "grey50", segment.size = 0.5,
    box.padding = 0.1, point.padding = 0.3
  ) +
  labs(title = "Genome Composition of Haemonchus contortus",
       fill = "Category",
       y = "Percent (%)") +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) -> p
ggsave("hc_genome_composition_pie_labels.png", plot = p, width = 7, height = 5, dpi = 300)


# -------------------------------
# 11. Comparative Plots Across Species
# -------------------------------

# TE copy number comparison
p <- ggplot(combined_copy, aes(x = reorder(te_class, copy_number), y = copy_number, fill = species)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = species_colors) +
  labs(title = "TE Copy Number Across Species", x = "TE Class", y = "Copy Number") +
  clean_theme()
ggsave("comparison_TE_copy_number.png", plot = p, width = 10, height = 6)

# Genome coverage comparison
p <- ggplot(combined_summary, aes(x = reorder(class, genome_percent), y = genome_percent, fill = species)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = species_colors) +
  labs(title = "Genome Coverage by TE Class Across Species", x = "TE Class", y = "Percent Genome (%)") +
  clean_theme()
ggsave("comparison_TE_genome_coverage.png", plot = p, width = 10, height = 6)

# TE Kimura Divergence Across Species (faceted by species)
hc_gff$species <- "H. contortus"
tc_gff$species <- "T. circumcincta"
ce_gff$species <- "C. elegans"

kimura_all <- bind_rows(hc_gff, tc_gff, ce_gff) %>%
  filter(!is.na(kimura))

p_kimura_all <- ggplot(kimura_all, aes(x = kimura, fill = species)) +
  geom_histogram(binwidth = 0.01, alpha = 0.8) +
  facet_wrap(~ species) +
  scale_fill_manual(values = species_colors) +
  labs(
    title = "TE Kimura Divergence by Species",
    x = "Kimura Divergence",
    y = "TE Count"
  ) +
  clean_theme()

ggsave("kimura_divergence_all_species.png", plot = p_kimura_all, width = 15, height = 6)

# ---------------------------------------------
# Script 2 - Analysis of TE Insertions identified by PoPoolationTE2

# -------------------------------
# 1. Load required libraries
# -------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(broom)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(scales)
library(ComplexUpset)
library(UpSetR)
library(GenomicRanges)
library(rtracklayer)
library(grid)
library(pheatmap)
library(stringr)
library(uwot)
library(ggfortify)
library(gtools)
library(tidyverse)

# -------------------------------
# 2. Define plotting theme and colors
# -------------------------------

# Custom clean theme for all plots
clean_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Color palette for TE classes
te_colors <- c(
  "DNA" = "#6baed6", 
  "LINE" = "#fdae6b",
  "LTR" = "#74c476",
  "PLE" = "#fb6a4a",
  "RC" = "#9e9ac8",
  "Non-repetitive genome" = "grey70"
)

# -------------------------------
# 3. Functions for parsing TE classes
# -------------------------------

# Parse TE classes from EarlGrey output
parse_earlgrey_classes <- function(dt) {
  dt[, c("main_class", "raw_subclass") := tstrsplit(te_classes, "/", fixed = TRUE)]
  dt[, c("te_class", "te_subclass") := {
    sub_parts <- tstrsplit(raw_subclass, "-", fixed = TRUE)
    class_part <- sub_parts[[1]]
    subclass_part <- ifelse(is.na(sub_parts[[2]]), class_part, sub_parts[[2]])
    .(class_part, subclass_part)
  }]
  dt[, te_class := fifelse(is.na(te_class), "Unknown", te_class)]
  return(dt)
}

# -------------------------------
# 4. PoPoolationTE2 processing
# -------------------------------

# Sample names
sample_names <- c(
  "1-14031_postCTL3", "2-14036_postCTL3", "3-14949_postCTL3",
  "4-14907_postIVM3", "5-14916_postIVM3", "6-14917_postIVM3",
  "7-14029_postMOX3", "8-14953_postMOX3", "9-14954_postMOX3"
)

# Load joint PoPoolationTE2 results
popte <- fread("joint.tes.finalresult.txt", header = FALSE)
colnames(popte)[1:8] <- c("GroupID", "Chrom", "Position", "Strand", "TE_ID", "TE_Class", "Signature", "Comments")
colnames(popte)[9:17] <- sample_names

# Pivot to long format
popte_long <- popte %>%
  pivot_longer(cols = all_of(sample_names), names_to = "SampleName", values_to = "Frequency") %>%
  mutate(
    TE_suffix = str_extract(TE_ID, "_\\d+$") %>% str_remove("_"),
    TE_full = gsub("_[0-9]+$", "", TE_ID),
    TE_subclass = str_split_fixed(TE_full, "-", 2)[, 1],
    TE_family = str_split_fixed(TE_full, "-", 2)[, 2],
    TE_family = ifelse(is.na(TE_family) | TE_family == "", TE_subclass, TE_family)
  )

# Assign sample groups
sample_groups <- tibble(
  SampleName = sample_names,
  Group = rep(c("control", "ivermectin", "moxidectin"), each = 3)
)
popte_long <- popte_long %>% left_join(sample_groups, by = "SampleName")

# Calculate presence/absence and categories
haplotypes_per_sample <- 400
popte_long <- popte_long %>%
  mutate(
    Start = pmax(Position - 25, 1),
    End = Position + 25,
    coord = paste(Chrom, Start, End, sep = "_"),
    present = round(Frequency * haplotypes_per_sample),
    absent = haplotypes_per_sample - present,
    present_bin = ifelse(Frequency > 0.1, 1, 0),
    sample_type = case_when(
      Group == "control" ~ "control",
      Group %in% c("ivermectin", "moxidectin") ~ "treated",
      TRUE ~ NA_character_
    ),
    Category = case_when(
      Frequency < 0.10 ~ "LowFreq",
      Frequency > 0.90 ~ "HighFreq",
      TRUE ~ "Polymorphic"
    )
  )

# -------------------------------
# 5. Per-sample TE Insertion Counts
# -------------------------------

te_counts <- popte_long %>%
  group_by(SampleName, Group) %>%
  summarise(NumInsertions = n(), .groups = "drop")

ggplot(te_counts, aes(x = SampleName, y = NumInsertions, fill = Group)) +
  geom_col() +
  labs(
    title = "Number of TE insertions per sample",
    y = "Number of TE insertions"
  ) +
  scale_fill_manual(values = c(
    "control" = "#1b9e77",
    "ivermectin" = "#d95f02",
    "moxidectin" = "#7570b3"
  )) +
  clean_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("TE_insertions_per_sample.png", width = 8, height = 5)

# -------------------------------
# 6. Compare TE insertions between samples and reference genome
# -------------------------------
# EarlGrey TE dataset

earlgrey <- fread("Hc_TE_annotation.1kb.bed", header = FALSE)
setnames(earlgrey, c("Chrom", "Start", "End", "te_classes", "Score", "Strand"))

# Apply the parsing function
earlgrey <- parse_earlgrey_classes(as.data.table(earlgrey))

# Add metadata columns
earlgrey[, Length := End - Start]
earlgrey[, Source := "EarlGrey"]
earlgrey[, coord := paste(Chrom, Start, End, sep = "_")]

# Count TE insertions per class for plotting
eg_class_group <- earlgrey %>%
  count(te_class) %>%
  rename(TE_Class = te_class, n = n) %>%
  mutate(Group = "EarlGrey")

# Count TE insertions per class for PoPoolationTE2 by Group
pt_class_group <- popte_long %>%
  group_by(Group, te_subclass) %>%
  summarise(n = n(), .groups = "drop") %>%
  rename(TE_Class = te_subclass) %>%
  mutate(TE_Class = ifelse(is.na(TE_Class) | TE_Class == "", "Unknown", TE_Class))

# Combine both datasets
combined_class_group <- bind_rows(pt_class_group, eg_class_group)

# -------------------------------
# Plot Barplot of TE Insertions per Class by Group (PoPoolationTE2 + EarlGrey)
# -------------------------------

ggplot(combined_class_group, aes(x = fct_reorder(TE_Class, n, .desc = TRUE), y = n, fill = Group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "TE Insertions per Class by Group (PoPoolationTE2 + EarlGrey)",
    x = "TE Class",
    y = "Number of Insertions"
  ) +
  scale_fill_manual(values = c(
    "control" = "#1b9e77",
    "ivermectin" = "#d95f02",
    "moxidectin" = "#7570b3",
    "EarlGrey" = "lightblue"
  )) +
  clean_theme()

ggsave("insertions_per_class_by_group_with_earlgrey.png", width = 10, height = 7)

# ----------------------------------------
# 7. PCA on TE Family Presence (All Insertions)
# ----------------------------------------

# Create presence/absence matrix of TE families
te_matrix_all <- popte_long %>%
  group_by(SampleName, TE_family) %>%
  summarise(Present = 1, .groups = "drop") %>%
  pivot_wider(names_from = TE_family, values_from = Present, values_fill = 0)

# Filter out columns with zero variance
te_matrix_filtered_all <- te_matrix_all %>%
  select(-SampleName) %>%
  select(where(~ var(.) > 0))

# Run PCA only if matrix has more than one TE family
if (ncol(te_matrix_filtered_all) > 1) {
  te_pca_all <- prcomp(te_matrix_filtered_all, scale. = TRUE)
  
  autoplot(te_pca_all, data = te_matrix_all, label = TRUE, label.size = 4) +
    clean_theme() +
    labs(title = "PCA: TE Family Presence (All Insertions)")
  
  ggsave("pca_te_family_all.png", width = 10, height = 6)
}

# ----------------------------------------
# 8. Histogram of TE insertion frequencies (10–90%)
# ----------------------------------------

ggplot(popte_long %>% filter(Frequency >= 0.10, Frequency <= 0.90), 
       aes(x = Frequency, fill = Group)) +
  geom_histogram(binwidth = 0.05, color = "black", alpha = 0.7, position = "identity") +
  facet_wrap(~Group, ncol = 3) +
  labs(title = "TE Insertion Frequency Distribution by Group (PoPoolationTE2, 10–90%)",
       x = "Population Frequency", y = "Count") +
  scale_fill_manual(values = c("control" = "#1b9e77", "ivermectin" = "#d95f02", "moxidectin" = "#7570b3")) +
  clean_theme()

ggsave("PoPoolationTE2_Frequency_Distribution_10to90.png", width = 10, height = 5, dpi = 300)

# ----------------------------------------
# 9. Prepare binary presence/absence matrix for TE insertions for Complex UpSet plots
# ----------------------------------------
binary_matrix_long <- popte_long %>%
  filter(Frequency > 0.10) %>%  # filter insertions with frequency > 10%
  distinct(SampleName, coord) %>%  # remove duplicates
  mutate(present = 1) %>%  # mark presence
  pivot_wider(
    names_from = SampleName, 
    values_from = present, 
    values_fill = 0  # absent = 0
  )

# Define sample groups and colors
sample_groups <- data.frame(
  sample = colnames(binary_matrix_long)[-1],  # exclude coord column
  group  = c(rep("Control", 3), rep("Ivermectin", 3), rep("Moxidectin", 3))
)
rownames(sample_groups) <- sample_groups$sample

group_colors <- c(
  Control = "#1b9e77",
  Ivermectin = "#d95f02",
  Moxidectin = "#7570b3"
)

# Create color-coded queries for each sample
query_list <- lapply(1:nrow(sample_groups), function(i) {
  upset_query(
    set  = sample_groups$sample[i],
    fill = group_colors[sample_groups$group[i]]
  )
})

# Define function to generate Upset plots
make_upset <- function(data, samples, min_size, max_size = 16333, title_suffix = "") {
  ComplexUpset::upset(
    data,
    intersect = samples,
    name = "TE Insertion",
    min_size = min_size,
    max_size = max_size,
    base_annotations = list(
      'Intersection size' = intersection_size(text = list(size = 3)) +
        labs(title = paste0("Shared TE Insertions ", title_suffix))
    ),
    queries = query_list,
    themes = upset_modify_themes(list(
      'intersections_matrix' = theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 3)
      ),
      'overall_sizes' = theme(
        axis.text.y = element_text(size = 7)
      )
    ))
  )
}

# Generate Upset plots for different frequency bins
plots <- list(
  low  = make_upset(binary_matrix_long, sample_groups$sample, 1, 15, "(1–15)"),
  mid1 = make_upset(binary_matrix_long, sample_groups$sample, 16, 100, "(16–100)"),
  high = make_upset(binary_matrix_long, sample_groups$sample, 101, 16333, "(>100)"),
  
  bin_16_25  = make_upset(binary_matrix_long, sample_groups$sample, 16, 25, "(16-25)"),
  bin_26_50  = make_upset(binary_matrix_long, sample_groups$sample, 26, 50, "(25–50)"),
  bin_51_70 = make_upset(binary_matrix_long, sample_groups$sample, 51, 70, "(51–70)"),
  bin_71_500 = make_upset(binary_matrix_long, sample_groups$sample, 71, 500, "(71–500)"),
  bin_over500 = make_upset(binary_matrix_long, sample_groups$sample, 501, 16333, "(>500)")
)

# Create an output folder
dir.create("upset_plots", showWarnings = FALSE)

# Loop over the list and save each plot separately
for (nm in names(plots)) {
  ggsave(
    filename = paste0("upset_plots/upset_", nm, ".png"),
    plot = plots[[nm]],
    width = 14,   # make big enough for labels
    height = 6,
    dpi = 300
  )
}

# ----------------------------------------
# 10. Logistic Regression
# ----------------------------------------

# Aggregate TE insertions by sample type and group
agg_sample_type <- popte_long %>%
  filter(!is.na(sample_type)) %>%
  group_by(coord, sample_type) %>%
  summarise(
    present = sum(present),
    absent  = sum(absent),
    .groups = "drop"
  )

agg_group <- popte_long %>%
  filter(!is.na(Group)) %>%
  group_by(coord, Group) %>%
  summarise(
    present = sum(present),
    absent  = sum(absent),
    .groups = "drop"
  )


# Filter only coordinates with variation between categories
filtered_sample_type <- agg_sample_type %>%
  group_by(coord) %>%
  filter(n_distinct(sample_type) > 1) %>%
  nest()

filtered_group <- agg_group %>%
  group_by(coord) %>%
  filter(n_distinct(Group) > 1) %>%
  nest()


# Fit logistic regression models and extract p-values
model_sample_type <- filtered_sample_type %>%
  mutate(
    model = map(data, ~ tryCatch(glm(cbind(present, absent) ~ sample_type, binomial, data = .x),
                                 error = function(e) NULL)),
    result = map(model, ~ if (!is.null(.x)) tidy(anova(.x, test = "LRT")) else NULL),
    p_sample_type = map_dbl(result, ~ if (!is.null(.x) && any(.x$term == "sample_type")) 
      .x$p.value[.x$term == "sample_type"] else NA_real_),
    padj_sample_type = p.adjust(p_sample_type, method = "BH")
  ) %>%
  select(coord, p_sample_type, padj_sample_type)

model_group <- filtered_group %>%
  mutate(
    model = map(data, ~ tryCatch(glm(cbind(present, absent) ~ Group, binomial, data = .x),
                                 error = function(e) NULL)),
    result = map(model, ~ if (!is.null(.x)) tidy(anova(.x, test = "LRT")) else NULL),
    p_group = map_dbl(result, ~ if (!is.null(.x) && any(.x$term == "Group")) 
      .x$p.value[.x$term == "Group"] else NA_real_),
    padj_group = p.adjust(p_group, method = "BH")
  ) %>%
  select(coord, p_group, padj_group)

# Merge sample_type and group results
final_results <- full_join(model_sample_type, model_group, by = "coord") %>%
  arrange(padj_sample_type)

# Separate nuclear and mitochondrial variants
nuclear_results <- final_results %>%
  filter(grepl("^hcontortus_chr", coord))

mito_results <- final_results %>%
  filter(grepl("^mitochondrion_", coord))

# Prepare genomic positions for Manhattan plots
merged_results <- nuclear_results %>%
  extract(
    col = coord,
    into = c("chr", "start", "end"),
    regex = "hcontortus_(chr[^_]+)_.*_(\\d+)_(\\d+)"
  ) %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end),
    chr   = factor(chr, levels = mixedsort(unique(chr))),
    genomic_pos = (start + end) / 2
  )

chr_offsets <- merged_results %>%
  group_by(chr) %>%
  summarise(chr_len = max(genomic_pos, na.rm = TRUE)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))

merged_results <- merged_results %>%
  left_join(chr_offsets, by = "chr") %>%
  mutate(genomic_pos = genomic_pos + chr_start)

# Avoid log10(0)
merged_results <- merged_results %>%
  mutate(
    p_sample_type = ifelse(p_sample_type == 0, 1e-300, p_sample_type),
    p_group       = ifelse(p_group == 0, 1e-300, p_group)
  )

saveRDS(merged_results, file = "merged_results.rds")

# Generate Manhattan plots
manhattan_sample_type <- ggplot(merged_results, aes(x = genomic_pos, y = -log10(p_sample_type))) +
  geom_point(aes(color = chr), alpha = 0.75, size = 1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = length(unique(merged_results$chr)), name = "Paired")) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = chr_offsets$chr_start + chr_offsets$chr_len / 2,
    labels = chr_offsets$chr
  ) +
  labs(y = expression(-log[10](italic(p))), title = "Manhattan Plot: sample_type") +
  ylim(0, 350) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  clean_theme()

ggsave("manhattan_sample_type.png", manhattan_sample_type, width = 10, height = 6, dpi = 300)

manhattan_group <- ggplot(merged_results, aes(x = genomic_pos, y = -log10(p_group))) +
  geom_point(aes(color = chr), alpha = 0.75, size = 1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = length(unique(merged_results$chr)), name = "Paired")) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = chr_offsets$chr_start + chr_offsets$chr_len / 2,
    labels = chr_offsets$chr
  ) +
  labs(y = expression(-log[10](italic(p))), title = "Manhattan Plot: group") +
  ylim(0, 350) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  clean_theme()

ggsave("manhattan_group.png", manhattan_group, width = 10, height = 6, dpi = 300)

