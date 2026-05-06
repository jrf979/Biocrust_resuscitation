#title: "16S data analysis, including stats and plots"
#author: "Raul Roman"
#date: "2025-10-08"
#note: This script work with phyloseq_resuscitation.rds, a phyloseq object containing the ASV table, taxonomy table and associated metadata used for 16S data analysis.

rm(list = ls())
gc()

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggh4x)
library(grid)

# ---- Load phyloseq_resuscitation.rds object ----

setwd()

ps <- readRDS("phyloseq_resuscitation.rds")

ps

# remove ASVs belonging to organelles
ps_filt <- subset_taxa(
  ps,
  (is.na(Family)  | Family  != "Mitochondria") &
    (is.na(Order) | Order   != "Chloroplast")
)

ps_filt <- subset_samples(ps_filt, Name2 != "LSA4") # remove outlier

ps_filt <- prune_taxa(taxa_sums(ps_filt) > 0, ps_filt)

ps_filt

nsamples(ps_filt)
ntaxa(ps_filt)

table(sample_data(ps_filt)$Name2) # Two samples per microcosms, one for each Fraction: Boncat-active or Bulk DNA



# ---- Relative abundance by Phylum (Fig S1) ----

# Collapse to phylum

ps_phylum <- tax_glom(ps_filt, taxrank = "Phylum", NArm = FALSE)

ps_rel <- transform_sample_counts(
  ps_phylum,
  function(x) if (sum(x) > 0) x / sum(x) else x
)

df_rel <- psmelt(ps_rel)

df_rel <- df_rel %>%
  mutate(
    FractionLab = recode(
      Fraction,
      "boncat" = "BONCAT-active",
      "bulk"   = "Bulk DNA"
    ),
    Biocrust = factor(Biocrust, levels = c("L-BSC", "D-BSC")),
    Treatment = factor(Treatment, levels = c("Sun", "Night")),
    FractionLab = factor(FractionLab, levels = c("BONCAT-active", "Bulk DNA")),
    Panel = factor(
      paste(Biocrust, "•", FractionLab),
      levels = c(
        "L-BSC • BONCAT-active",
        "L-BSC • Bulk DNA",
        "D-BSC • BONCAT-active",
        "D-BSC • Bulk DNA"
      )
    )
  )


df_rel <- df_rel %>%
  arrange(Name2, Fraction, Treatment) %>%
  mutate(
    Sample = factor(Sample, levels = unique(Sample))
  )


topN <- df_rel %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 9) %>%
  pull(Phylum)

df_rel <- df_rel %>%
  mutate(
    Phylum_plot = ifelse(Phylum %in% topN, as.character(Phylum), "Other"),
    Phylum_plot = factor(Phylum_plot, levels = c(topN, "Other"))
  )

legend_breaks <- c(topN, "Other")

pal <- ggsci::pal_npg("nrc")(length(legend_breaks))
names(pal) <- legend_breaks


# Plot

p_rel_combined <- ggplot(
  df_rel,
  aes(x = Sample, y = Abundance, fill = Phylum_plot)
) +
  geom_col(width = 0.9) +
  
  ggh4x::facet_nested(
    . ~ Panel + Treatment,
    scales = "free_x",
    independent = "x"
  ) +
  
  scale_x_discrete(
    drop = TRUE,
    labels = function(x) {
      df_rel %>%
        distinct(Sample, Name2) %>%
        filter(Sample %in% x) %>%
        arrange(match(Sample, x)) %>%
        pull(Name2)
    }
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  scale_fill_manual(
    values = pal,
    breaks = legend_breaks,
    name = "Phylum"
  ) +
  
  labs(
    x = NULL,
    y = "Relative abundance"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 6
    ),
    axis.ticks.x = element_blank(),
    
    strip.text = element_text(face = "bold"),
    strip.background = element_blank(),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.height = grid::unit(4, "mm"),
    legend.key.width  = grid::unit(8, "mm")
  )

p_rel_combined


# ---- Remove Cyanobacteria. This object will be used from now on ----

ps_no_cyano <- subset_taxa(
  ps_filt,
  !(Phylum %in% c("Cyanobacteria", "Cyanobacteriota"))
)

ps_no_cyano <- prune_samples(sample_sums(ps_no_cyano) > 0, ps_no_cyano)
ps_no_cyano <- prune_taxa(taxa_sums(ps_no_cyano) > 0, ps_no_cyano)


# ---- Fig S2 and S3 -> TOP6 Phyla per fraction and biocrust -----

PAL_FRACTION <- c(
  "BONCAT-active" = "#1b9e77",
  "Bulk DNA"      = "#d95f02"
)

prep_top6_phyla <- function(ps_in, biocrust_value) {
  
  ps_b <- prune_samples(sample_data(ps_in)$Biocrust == biocrust_value, ps_in)
  ps_b <- prune_taxa(taxa_sums(ps_b) > 0, ps_b)
  
  ps_phylum <- tax_glom(ps_b, taxrank = "Phylum", NArm = TRUE)
  
  ps_rel <- transform_sample_counts(
    ps_phylum,
    function(x) if (sum(x) > 0) x / sum(x) else x
  )
  
  df <- psmelt(ps_rel) %>%
    mutate(
      Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", as.character(Phylum)),
      FractionLab = recode(
        Fraction,
        "boncat" = "BONCAT-active",
        "bulk"   = "Bulk DNA"
      ),
      FractionLab = factor(FractionLab, levels = c("BONCAT-active", "Bulk DNA"))
    )
  
  top6 <- df %>%
    group_by(Phylum) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = 6) %>%
    pull(Phylum)
  
  df_top <- df %>%
    filter(Phylum %in% top6) %>%
    mutate(Phylum = factor(Phylum, levels = top6))
  
  df_top
}


add_pvals <- function(df_in) {
  
  df_in %>%
    group_by(Phylum) %>%
    summarise(
      p = if (n_distinct(FractionLab) > 1) {
        wilcox.test(Abundance ~ FractionLab, exact = FALSE, correct = TRUE)$p.value
      } else {
        NA_real_
      },
      y = max(Abundance, na.rm = TRUE) * 1.08,
      .groups = "drop"
    ) %>%
    mutate(
      label = ifelse(is.na(p), "n/a", paste0("p = ", signif(p, 3)))
    )
}

df_L <- prep_top6_phyla(ps_no_cyano, "L-BSC")
df_D <- prep_top6_phyla(ps_no_cyano, "D-BSC")

pv_L <- add_pvals(df_L)
pv_D <- add_pvals(df_D)

# plot L-BSC (Figure S2)

p_L_top6 <- ggplot(
  df_L,
  aes(x = FractionLab, y = Abundance, fill = FractionLab)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.7,
    color = "black"
  ) +
  geom_point(
    position = position_jitter(width = 0.12),
    size = 1.6,
    alpha = 0.75
  ) +
  geom_text(
    data = pv_L,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~ Phylum, ncol = 3, scales = "free_y") +
  scale_y_continuous(
    breaks = seq(0, 0.4, 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 0.45)) +
  scale_fill_manual(values = PAL_FRACTION, guide = "none") +
  labs(
    x = NULL,
    y = "Relative abundance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)
  )

p_L_top6


# plot D-BSC (Figure S3)

p_D_top6 <- ggplot(
  df_D,
  aes(x = FractionLab, y = Abundance, fill = FractionLab)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.7,
    color = "black"
  ) +
  geom_point(
    position = position_jitter(width = 0.12),
    size = 1.6,
    alpha = 0.75
  ) +
  geom_text(
    data = pv_D,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~ Phylum, ncol = 3, scales = "free_y") +
  scale_y_continuous(
    breaks = seq(0, 0.4, 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 0.45)) +
  scale_fill_manual(values = PAL_FRACTION, guide = "none") +
  labs(
    x = NULL,
    y = "Relative abundance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9)
  )

p_D_top6


# ---- Figure 3A - NMDS + Table1 - PERMANOVA ----

library(vegan)
library(scales)


PAL_TREAT <- c(
  "Sun"   = "#b39600",
  "Night" = "#66CCFF"
)


add_label_cols <- function(df_like) {
  df_like %>%
    mutate(
      FractionLab = recode(
        as.character(Fraction),
        "boncat" = "BONCAT-active",
        "bulk"   = "Bulk DNA"
      ),
      BiocrustLab = factor(Biocrust, levels = c("L-BSC", "D-BSC")),
      TreatmentLab = factor(Treatment, levels = c("Sun", "Night")),
      FractionLab = factor(FractionLab, levels = c("BONCAT-active", "Bulk DNA"))
    )
}


# Bray-Curtis Distance + NMDS

otu_mat <- as(otu_table(ps_no_cyano), "matrix")

if (taxa_are_rows(ps_no_cyano)) {
  otu_mat <- t(otu_mat)
}

dist_bray <- vegan::vegdist(otu_mat, method = "bray")

set.seed(123)
nmds_bray <- vegan::metaMDS(
  dist_bray,
  k = 2,
  trymax = 100,
  autotransform = FALSE,
  trace = FALSE
)

scores_unrare <- as.data.frame(vegan::scores(nmds_bray, display = "sites"))
scores_unrare$Sample <- rownames(scores_unrare)

meta <- as(sample_data(ps_no_cyano), "data.frame") %>%
  add_label_cols()

meta$Sample <- rownames(meta)

scores_unrare <- scores_unrare %>%
  left_join(meta, by = "Sample")


# Plot function

plot_nmds <- function(scores) {
  
  ggplot(scores, aes(NMDS1, NMDS2)) +
    
    stat_ellipse(
      aes(
        group = interaction(BiocrustLab, FractionLab),
        linetype = BiocrustLab
      ),
      color = "black",
      type = "t",
      level = 0.95,
      linewidth = 0.6,
      alpha = 0.9
    ) +
    
    geom_point(
      aes(color = TreatmentLab, shape = FractionLab),
      size = 3.2,
      alpha = 0.9,
      stroke = 0.5
    ) +
    
    scale_color_manual(
      values = PAL_TREAT,
      name = "Incubation"
    ) +
    
    scale_shape_manual(
      values = c(
        "BONCAT-active" = 15,
        "Bulk DNA" = 17
      ),
      name = "Fraction"
    ) +
    
    scale_linetype_manual(
      values = c(
        "L-BSC" = "solid",
        "D-BSC" = "33"
      ),
      name = "Biocrust"
    ) +
    
    scale_x_continuous(breaks = scales::breaks_width(0.5)) +
    scale_y_continuous(breaks = scales::breaks_width(0.5)) +
    
    labs(
      x = "NMDS1",
      y = "NMDS2"
    ) +
    
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      legend.position = "right"
    )
}


# NMDS plot

p_nmds <- plot_nmds(scores_unrare)

p_nmds

nmds_bray$stress


# PERMANOVA 


meta_perm <- meta[labels(dist_bray), , drop = FALSE]

set.seed(123)
perm_bray_unrare <- vegan::adonis2(
  dist_bray ~ (TreatmentLab + BiocrustLab + FractionLab)^2,
  data = meta_perm,
  permutations = 999,
  by = "terms"
)

perm_bray_unrare


# ---- Figure 3B - Distance to centroid + EMMs ----

library(emmeans)
library(multcomp)
library(tibble)

meta_perm$G_All <- interaction(
  meta_perm$TreatmentLab,
  meta_perm$BiocrustLab,
  meta_perm$FractionLab,
  drop = TRUE
)


grp <- factor(meta_perm$G_All)

bd_all <- vegan::betadisper(
  dist_bray,
  grp,
  type = "centroid"
)

# extract distance-to-centroid table
df_disp <- data.frame(
  Sample = names(bd_all$distances),
  Dist2Centroid = unname(bd_all$distances)
) %>%
  left_join(
    meta_perm %>%
      dplyr::select(Sample, BiocrustLab, FractionLab, TreatmentLab),
    by = "Sample"
  ) %>%
  mutate(
    BiocrustLab  = factor(BiocrustLab, levels = c("L-BSC", "D-BSC")),
    FractionLab  = factor(FractionLab, levels = c("BONCAT-active", "Bulk DNA")),
    TreatmentLab = factor(TreatmentLab, levels = c("Sun", "Night"))
  )

# LM + EMMs on distance-to-centroid
lm_disp <- lm(
  Dist2Centroid ~ (BiocrustLab + TreatmentLab + FractionLab)^2,
  data = df_disp
)

anova(lm_disp)

emm_bt_byF <- emmeans(
  lm_disp,
  ~ BiocrustLab * TreatmentLab | FractionLab
)

pairs(emm_bt_byF, adjust = "tukey")

cld_bt_byF <- multcomp::cld(
  emm_bt_byF,
  adjust = "tukey",
  Letters = letters
)

df_emm <- as.data.frame(cld_bt_byF) %>%
  mutate(
    Letters = gsub("\\s+", "", .group),
    BiocrustLab  = factor(BiocrustLab, levels = c("L-BSC", "D-BSC")),
    FractionLab  = factor(FractionLab, levels = c("BONCAT-active", "Bulk DNA")),
    TreatmentLab = factor(TreatmentLab, levels = c("Sun", "Night"))
  ) %>%
  group_by(FractionLab) %>%
  mutate(y_lab = max(emmean + SE, na.rm = TRUE) * 1.06) %>%
  ungroup()

# EMM plot
pd <- position_dodge2(width = 0.65, preserve = "single")

p_emm_dist_centroid <- ggplot(
  df_emm,
  aes(x = BiocrustLab, y = emmean, fill = TreatmentLab, color = TreatmentLab)
) +
  geom_linerange(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    color = "black",
    linewidth = 0.3,
    position = pd
  ) +
  geom_point(size = 3.6, color = "black", position = pd) +
  geom_point(size = 2.8, position = pd) +
  geom_text(
    aes(label = Letters, y = y_lab, group = TreatmentLab),
    color = "black",
    size = 4,
    vjust = 0,
    position = pd
  ) +
  facet_wrap(~ FractionLab, nrow = 1) +
  scale_fill_manual(values = PAL_TREAT, name = "Incubation") +
  scale_color_manual(values = PAL_TREAT, name = "Incubation") +
  labs(
    x = NULL,
    y = "β-dispersion (EMM distance to centroid)"
  ) +
  scale_y_continuous(
    limits = c(0.2, 0.5),
    breaks = seq(0.2, 0.5, 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

p_emm_dist_centroid


# ---- Combine Figure 3A and 3B ----

library(cowplot)

p_nmds_clean <- p_nmds +
  labs(
    color = NULL,
    shape = NULL,
    linetype = NULL
  ) +
  guides(
    color = guide_legend(title = NULL),
    shape = guide_legend(title = NULL),
    linetype = guide_legend(title = NULL)
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

leg_shared <- cowplot::get_legend(p_nmds_clean)

pA <- p_nmds_clean +
  theme(legend.position = "none")

pB <- p_emm_dist_centroid +
  theme(legend.position = "none") +
  labs(
    y = expression(beta*"-dispersion (EMM distance to centroid)")
  )

pA_labeled <- cowplot::ggdraw(pA) +
  cowplot::draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
                      fontface = "bold", size = 14)

pB_labeled <- cowplot::ggdraw(pB) +
  cowplot::draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
                      fontface = "bold", size = 14)

fig3_top <- cowplot::plot_grid(
  pA_labeled,
  pB_labeled,
  nrow = 1,
  rel_widths = c(1.55, 1)
)

fig3 <- cowplot::plot_grid(
  fig3_top,
  leg_shared,
  ncol = 1,
  rel_heights = c(1, 0.14)
)

fig3
