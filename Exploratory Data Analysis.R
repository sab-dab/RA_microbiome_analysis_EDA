setwd("C:/Users/sabir/Downloads/gut microbiota")

library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(ggplot2)
library(vegan)

# -----------------------------
# 1) Load data
# -----------------------------
asv_profile <- readRDS("1.ASV.profile.rds")
tax_info    <- readRDS("1.taxonomy.info.rds")

asv_profile_df <- asv_profile %>%
  as.data.frame() %>%
  rownames_to_column("ASV")

tax_info_df <- tax_info %>%
  as.data.frame() %>%
  rownames_to_column("ASV")

sample_cols <- setdiff(colnames(asv_profile_df), "ASV")

hc_cols <- grep("^HC", sample_cols, value = TRUE)
ra_cols <- grep("^RA", sample_cols, value = TRUE)

cat("HC samples:", length(hc_cols), "\n")
cat("RA samples:", length(ra_cols), "\n")
cat("Total ASVs:", nrow(asv_profile_df), "\n\n")

# -----------------------------
# 2) Objective 1: Dataset suitability
# -----------------------------
cat("Taxonomy rows:", nrow(tax_info_df), "\n")
cat("ASV overlap between profile and taxonomy:",
    sum(asv_profile_df$ASV %in% tax_info_df$ASV), "\n\n")

# -----------------------------
# 3) Objective 2: Data quality checks
# -----------------------------
cat("Missing values in ASV profile:", sum(is.na(asv_profile_df)), "\n")
cat("Missing values in taxonomy:", sum(is.na(tax_info_df)), "\n\n")

# Total sequencing depth per sample
asv_mat <- asv_profile_df %>% select(-ASV) %>% as.matrix()
seq_depth <- colSums(asv_mat, na.rm = TRUE)

depth_df <- data.frame(
  Sample = names(seq_depth),
  Depth = as.numeric(seq_depth),
  Group = ifelse(grepl("^HC", names(seq_depth)), "HC", "RA")
)

ggplot(depth_df, aes(x = Group, y = Depth)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_minimal() +
  labs(title = "Sequencing Depth per Sample", y = "Total Reads (Sum of ASVs)")

# Identify very low depth samples
low_depth_cutoff <- quantile(depth_df$Depth, 0.05)
cat("Low-depth cutoff (5th percentile):", low_depth_cutoff, "\n")
cat("Number of low-depth samples:", sum(depth_df$Depth < low_depth_cutoff), "\n\n")

# -----------------------------
# 4) Objective 3: Diversity + separability
# -----------------------------

# Richness (ASV count per sample)
richness <- colSums(asv_mat > 0, na.rm = TRUE)

rich_df <- data.frame(
  Sample = names(richness),
  Richness = as.numeric(richness),
  Group = ifelse(grepl("^HC", names(richness)), "HC", "RA")
)

ggplot(rich_df, aes(x = Group, y = Richness)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_minimal() +
  labs(title = "ASV Richness (Alpha Diversity)", y = "Number of ASVs")

# Shannon diversity
shannon <- vegan::diversity(t(asv_mat), index = "shannon")

div_df <- data.frame(
  Sample = names(shannon),
  Shannon = as.numeric(shannon),
  Group = ifelse(grepl("^HC", names(shannon)), "HC", "RA")
)

ggplot(div_df, aes(x = Group, y = Shannon)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.4) +
  theme_minimal() +
  labs(title = "Shannon Diversity (Alpha Diversity)", y = "Shannon Index")

# -----------------------------
# 5) Genus-level EDA
# -----------------------------
tax_abund <- tax_info_df %>%
  left_join(asv_profile_df, by = "ASV") %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+")) %>%
  mutate(Genus = ifelse(is.na(Genus), "g__Unknown", Genus))

genus_abund <- tax_abund %>%
  select(Genus, all_of(sample_cols)) %>%
  group_by(Genus) %>%
  summarise(across(all_of(sample_cols), \(x) sum(x, na.rm = TRUE)), .groups = "drop")

# Relative abundance
genus_rel <- genus_abund
genus_rel[, sample_cols] <- sweep(genus_rel[, sample_cols], 2,
                                  colSums(genus_rel[, sample_cols]), "/")

# Mean relative abundance per group
genus_summary <- genus_rel %>%
  mutate(
    mean_HC = rowMeans(across(all_of(hc_cols)), na.rm = TRUE),
    mean_RA = rowMeans(across(all_of(ra_cols)), na.rm = TRUE),
    diff_RA_minus_HC = mean_RA - mean_HC
  ) %>%
  arrange(desc(abs(diff_RA_minus_HC)))

print(head(genus_summary, 20))

# -----------------------------
# 6) Beta diversity (PCoA)
# -----------------------------
bray <- vegan::vegdist(t(asv_mat), method = "bray")
pcoa <- cmdscale(bray, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  Sample = rownames(pcoa$points),
  PC1 = pcoa$points[,1],
  PC2 = pcoa$points[,2],
  Group = ifelse(grepl("^HC", rownames(pcoa$points)), "HC", "RA")
)

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(title = "PCoA (Bray-Curtis Distance)", x = "PC1", y = "PC2")

# PERMANOVA test (group difference)
meta <- data.frame(Group = pcoa_df$Group)
adonis_res <- vegan::adonis2(bray ~ Group, data = meta)
print(adonis_res)
adonis_table <- as.data.frame(adonis_res)
adonis_table$Term <- rownames(adonis_table)
adonis_table <- adonis_table %>% select(Term, everything())
adonis_table
write.csv(adonis_table, "PERMANOVA_results.csv", row.names = FALSE)

disp <- betadisper(bray, pcoa_df$Group)
anova(disp)
permutest(disp, permutations = 999)


# -----------------------------
# ASV prevalence across samples
# -----------------------------

# asv_mat is already: asv_profile_df %>% select(-ASV) %>% as.matrix()

# prevalence = fraction of samples where ASV is present (>0)
asv_prev_all <- rowMeans(asv_mat > 0, na.rm = TRUE)

prev_summary <- data.frame(
  Threshold = c(">= 0.5%", ">= 1%", ">= 5%", ">= 10%", ">= 20%"),
  Num_ASVs = c(
    sum(asv_prev_all >= 0.005),
    sum(asv_prev_all >= 0.01),
    sum(asv_prev_all >= 0.05),
    sum(asv_prev_all >= 0.10),
    sum(asv_prev_all >= 0.20)
  )
)

print(prev_summary)

# Optional: prevalence distribution plot
prev_df <- data.frame(Prevalence = asv_prev_all)

ggplot(prev_df, aes(x = Prevalence)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(title = "Distribution of ASV Prevalence (All Samples)",
       x = "Fraction of samples where ASV is present",
       y = "Number of ASVs")

# -----------------------------
# Top 20 genera barplot (HC vs RA)
# -----------------------------

# genus_rel already created earlier (relative abundance genus x samples)

genus_bar <- genus_rel %>%
  mutate(
    mean_HC = rowMeans(across(all_of(hc_cols)), na.rm = TRUE),
    mean_RA = rowMeans(across(all_of(ra_cols)), na.rm = TRUE)
  ) %>%
  select(Genus, mean_HC, mean_RA) %>%
  pivot_longer(cols = c(mean_HC, mean_RA),
               names_to = "Group",
               values_to = "MeanRelAbund") %>%
  mutate(Group = recode(Group, mean_HC = "HC", mean_RA = "RA"))

# Choose top genera by overall mean across both groups
top_genera <- genus_bar %>%
  group_by(Genus) %>%
  summarise(overall_mean = mean(MeanRelAbund), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  slice(1:20) %>%
  pull(Genus)

genus_bar_top <- genus_bar %>%
  filter(Genus %in% top_genera) %>%
  mutate(Genus = factor(Genus, levels = rev(top_genera)))

ggplot(genus_bar_top, aes(x = Genus, y = MeanRelAbund, fill = Group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Genera (Mean Relative Abundance): HC vs RA",
       x = "Genus",
       y = "Mean Relative Abundance")

p_top20 <- ggplot(genus_bar_top, aes(x = Genus, y = MeanRelAbund, fill = Group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Genera (Mean Relative Abundance): HC vs RA",
       x = "Genus",
       y = "Mean Relative Abundance")
ggsave("Figure_Top20_Genera.png", plot = p_top20,
       width = 10, height = 7, dpi = 300)
