#title: "ANCOM-BC2 data analysis and plots"
#author: "Raul Roman"
#date: "2025-10-08"
#note: This script works with Tables S3 and S4.

rm(list = ls())
gc()

## ---- libraries ----

library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggvenn)

## ---- paths ----

setwd()

input_dir  <- "data"
output_dir <- "figures"

if (!dir.exists(output_dir)) dir.create(output_dir)

bio_all <- readr::read_csv(
  file.path(input_dir, "TableS3.csv"),
  show_col_types = FALSE
)

frac_all <- readr::read_csv(
  file.path(input_dir, "TableS4.csv"),
  show_col_types = FALSE
)

## palette

alpha <- 0.05

PHYLUM_PALETTE <- c(
  "Pseudomonadota"    = "#1F77B4",
  "Actinomycetota"    = "#FF7F0E",
  "Bacteroidota"      = "#2CA02C",
  "Chloroflexota"     = "#D62728",
  "Firmicutes"        = "#9467BD",
  "Bacillota"         = "#9467BD",
  "Acidobacteriota"   = "#E377C2",
  "Planctomycetota"   = "#7F7F7F",
  "Verrucomicrobiota" = "#BCBD22",
  "Gemmatimonadota"   = "#17BECF",
  "Myxococcota"       = "#AEC7E8",
  "Bdellovibrionota"  = "#FFBB78",
  "Cyanobacteriota"   = "#98DF8A",
  "Armatimonadota"    = "#FF9896",
  "Abditibacteriota"  = "#C5B0D5",
  "Deinococcota"      = "#CCCCCC",
  "Unknown"           = "#999999"
)

## tables

bio_all <- bio_all %>%
  mutate(
    lfc = as.numeric(lfc),
    q = as.numeric(q),
    se = as.numeric(se),
    diff_robust = as.logical(diff_robust),
    Phylum = ifelse(is.na(Phylum) | !nzchar(Phylum), "Unknown", as.character(Phylum))
  )

frac_all <- frac_all %>%
  mutate(
    lfc = as.numeric(lfc),
    q = as.numeric(q),
    se = as.numeric(se),
    diff_robust = as.logical(diff_robust),
    Phylum = ifelse(is.na(Phylum) | !nzchar(Phylum), "Unknown", as.character(Phylum))
  )

# filter significant taxa

sig_biocrust <- bio_all %>%
  mutate(enriched_in = ifelse(lfc > 0, "BONCAT-active", "Bulk DNA")) %>%
  filter(!is.na(q), q < alpha, !is.na(diff_robust), diff_robust)

## ---- Figure 4 —---

## top 10 taxa enriched in each direction

d_fig4 <- sig_biocrust %>%
  filter(Biocrust %in% c("L-BSC", "D-BSC"), !is.na(lfc)) %>%
  filter(!grepl("_Unclassified$", label)) %>%
  group_by(Biocrust, enriched_in) %>%
  slice_max(abs(lfc), n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    xmin = ifelse(lfc >= 0, lfc, lfc - coalesce(se, 0)),
    xmax = ifelse(lfc >= 0, lfc + coalesce(se, 0), lfc)
  )

d_light_fig4 <- d_fig4 %>%
  filter(Biocrust == "L-BSC") %>%
  arrange(lfc) %>%
  mutate(label_f = factor(label, levels = label))

d_dark_fig4 <- d_fig4 %>%
  filter(Biocrust == "D-BSC") %>%
  arrange(lfc) %>%
  mutate(label_f = factor(label, levels = label))


phyla_fig4 <- sort(unique(d_fig4$Phylum))

PAL_FIG4 <- PHYLUM_PALETTE[phyla_fig4]
PAL_FIG4[is.na(PAL_FIG4)] <- PHYLUM_PALETTE["Unknown"]


xlim_fig4 <- 8
xbreaks_fig4 <- c(-8, -6, -4, -2, 0, 2, 4, 6, 8)

#Panel A : L-BSC

pA_fig4 <- ggplot(d_light_fig4, aes(x = lfc, y = label_f, fill = Phylum)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_col(width = 0.82, colour = "black", linewidth = 0.15) +
  geom_errorbarh(
    data = d_light_fig4 %>% filter(!is.na(se), is.finite(se), se > 0),
    aes(
      y = label_f,
      xmin = pmax(xmin, -xlim_fig4),
      xmax = pmin(xmax, xlim_fig4)
    ),
    height = 0.22,
    linewidth = 0.35,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    limits = names(PAL_FIG4),
    values = PAL_FIG4,
    drop = TRUE,
    name = "Phylum"
  ) +
  scale_x_continuous(
    limits = c(-xlim_fig4, xlim_fig4),
    breaks = xbreaks_fig4,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.25, 1.5))) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Log fold-change (BONCAT-active - Bulk DNA)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    legend.position = "none",
    plot.margin = margin(t = 22, r = 5, b = 5, l = 5)
  )

##Panel B: D-BSC

pB_fig4 <- ggplot(d_dark_fig4, aes(x = lfc, y = label_f, fill = Phylum)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_col(width = 0.82, colour = "black", linewidth = 0.15) +
  geom_errorbarh(
    data = d_dark_fig4 %>% filter(!is.na(se), is.finite(se), se > 0),
    aes(
      y = label_f,
      xmin = pmax(xmin, -xlim_fig4),
      xmax = pmin(xmax, xlim_fig4)
    ),
    height = 0.22,
    linewidth = 0.35,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    limits = names(PAL_FIG4),
    values = PAL_FIG4,
    drop = TRUE,
    name = "Phylum"
  ) +
  scale_x_continuous(
    limits = c(-xlim_fig4, xlim_fig4),
    breaks = xbreaks_fig4,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.25, 1.5))) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Log fold-change (BONCAT-active - Bulk DNA)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    legend.position = "none",
    plot.margin = margin(t = 22, r = 5, b = 5, l = 5)
  )

#legend

legend_fig4 <- cowplot::get_legend(
  ggplot(d_fig4, aes(x = lfc, y = label, fill = Phylum)) +
    geom_col() +
    scale_fill_manual(
      limits = names(PAL_FIG4),
      values = PAL_FIG4,
      drop = TRUE,
      name = "Phylum"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right")
)

# combine panels

pA_labeled <- cowplot::ggdraw(pA_fig4) +
  cowplot::draw_label(
    "A",
    x = 0.01,
    y = 0.99,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 14
  )

pB_labeled <- cowplot::ggdraw(pB_fig4) +
  cowplot::draw_label(
    "B",
    x = 0.01,
    y = 0.99,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 14
  )

fig4_panels <- cowplot::plot_grid(
  pA_labeled,
  pB_labeled,
  nrow = 1,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

fig4 <- cowplot::plot_grid(
  fig4_panels,
  legend_fig4,
  nrow = 1,
  rel_widths = c(1, 0.18)
)

fig4 #figure was further edited in .ppt to add enrichment arrows and labels



## ---- Figure 5 ----

## Figure 5A - venn

bio_light <- bio_all %>%
  filter(Biocrust == "L-BSC", !is.na(lfc)) %>%
  dplyr::select(taxon, lfc, q, label, Phylum) %>%
  rename(
    lfc_LBSC = lfc,
    q_LBSC = q,
    label_LBSC = label,
    Phylum_LBSC = Phylum
  )

bio_dark <- bio_all %>%
  filter(Biocrust == "D-BSC", !is.na(lfc)) %>%
  dplyr::select(taxon, lfc, q, label, Phylum) %>%
  rename(
    lfc_DBSC = lfc,
    q_DBSC = q,
    label_DBSC = label,
    Phylum_DBSC = Phylum
  )

bio_w <- inner_join(bio_light, bio_dark, by = "taxon") %>%
  filter(is.finite(lfc_LBSC), is.finite(lfc_DBSC)) %>%
  mutate(
    Phylum_plot = coalesce(Phylum_LBSC, Phylum_DBSC),
    label_plot = coalesce(label_LBSC, label_DBSC),
    sig_LBSC = !is.na(q_LBSC) & q_LBSC < alpha,
    sig_DBSC = !is.na(q_DBSC) & q_DBSC < alpha,
    sig_any = sig_LBSC | sig_DBSC,
    sig_both = sig_LBSC & sig_DBSC,
    alpha_pt = ifelse(sig_both, 1, 0.25)
  ) %>%
  filter(sig_any)

bio_sig <- bio_w %>%
  filter(sig_both)

cor_bio <- cor.test(
  bio_sig$lfc_LBSC,
  bio_sig$lfc_DBSC,
  method = "spearman",
  exact = FALSE
)

p_bio <- ifelse(
  cor_bio$p.value < 1e-4,
  "p < 1e-4",
  sprintf("p = %.3g", cor_bio$p.value)
)

ann_bio <- sprintf(
  "Spearman ρ = %.3f, %s (n = %d)",
  cor_bio$estimate,
  p_bio,
  nrow(bio_sig)
)

PAL_BIO <- PHYLUM_PALETTE[sort(unique(bio_w$Phylum_plot))]
PAL_BIO[is.na(PAL_BIO)] <- PHYLUM_PALETTE["Unknown"]

bio_lab_x <- bio_w %>%
  filter(sig_both) %>%
  slice_max(abs(lfc_LBSC), n = 3, with_ties = FALSE) %>%
  pull(taxon)

bio_lab_y <- bio_w %>%
  filter(sig_both) %>%
  slice_max(abs(lfc_DBSC), n = 3, with_ties = FALSE) %>%
  pull(taxon)

bio_lab <- bio_w %>%
  filter(taxon %in% unique(c(bio_lab_x, bio_lab_y)))

bio_light_ids <- unique(bio_light$taxon[!is.na(bio_light$q_LBSC) & bio_light$q_LBSC < alpha])
bio_dark_ids  <- unique(bio_dark$taxon[!is.na(bio_dark$q_DBSC) & bio_dark$q_DBSC < alpha])

venn_bio <- ggvenn(
  list(`L-BSC` = bio_light_ids, `D-BSC` = bio_dark_ids),
  fill_color = c("#F2E7D5", "#ADA9A0"),
  stroke_size = 0.4,
  set_name_size = 3.2,
  text_size = 3.5
)

## Figure 5B - venn

frac_boncat <- frac_all %>%
  filter(Fraction == "BONCAT", !is.na(lfc)) %>%
  dplyr::select(taxon, lfc, q, label, Phylum) %>%
  rename(
    lfc_BONCAT = lfc,
    q_BONCAT = q,
    label_BONCAT = label,
    Phylum_BONCAT = Phylum
  )

frac_bulk <- frac_all %>%
  filter(Fraction == "Bulk", !is.na(lfc)) %>%
  dplyr::select(taxon, lfc, q, label, Phylum) %>%
  rename(
    lfc_Bulk = lfc,
    q_Bulk = q,
    label_Bulk = label,
    Phylum_Bulk = Phylum
  )

fr_w <- inner_join(frac_boncat, frac_bulk, by = "taxon") %>%
  filter(is.finite(lfc_BONCAT), is.finite(lfc_Bulk)) %>%
  mutate(
    Phylum_plot = coalesce(Phylum_BONCAT, Phylum_Bulk),
    label_plot = coalesce(label_BONCAT, label_Bulk),
    sig_BONCAT = !is.na(q_BONCAT) & q_BONCAT < alpha,
    sig_Bulk = !is.na(q_Bulk) & q_Bulk < alpha,
    sig_any = sig_BONCAT | sig_Bulk,
    sig_both = sig_BONCAT & sig_Bulk,
    alpha_pt = ifelse(sig_both, 1, 0.25)
  ) %>%
  filter(sig_any)

fr_sig <- fr_w %>%
  filter(sig_both)

cor_fr <- cor.test(
  fr_sig$lfc_BONCAT,
  fr_sig$lfc_Bulk,
  method = "spearman",
  exact = FALSE
)

p_fr <- ifelse(
  cor_fr$p.value < 1e-4,
  "p < 1e-4",
  sprintf("p = %.3g", cor_fr$p.value)
)

ann_fr <- sprintf(
  "Spearman ρ = %.3f, %s (n = %d)",
  cor_fr$estimate,
  p_fr,
  nrow(fr_sig)
)

PAL_FR <- PHYLUM_PALETTE[sort(unique(fr_w$Phylum_plot))]
PAL_FR[is.na(PAL_FR)] <- PHYLUM_PALETTE["Unknown"]

fr_lab_x <- fr_w %>%
  filter(sig_both) %>%
  slice_max(abs(lfc_BONCAT), n = 3, with_ties = FALSE) %>%
  pull(taxon)

fr_lab_y <- fr_w %>%
  filter(sig_both) %>%
  slice_max(abs(lfc_Bulk), n = 3, with_ties = FALSE) %>%
  pull(taxon)

fr_lab <- fr_w %>%
  filter(taxon %in% unique(c(fr_lab_x, fr_lab_y)))

frac_boncat_ids <- unique(frac_boncat$taxon[!is.na(frac_boncat$q_BONCAT) & frac_boncat$q_BONCAT < alpha])
frac_bulk_ids   <- unique(frac_bulk$taxon[!is.na(frac_bulk$q_Bulk) & frac_bulk$q_Bulk < alpha])

venn_frac <- ggvenn(
  list(`BONCAT-active` = frac_boncat_ids, `Bulk DNA` = frac_bulk_ids),
  fill_color = c("#A9DACC", "#F1C3A0"),
  stroke_size = 0.4,
  set_name_size = 3.2,
  text_size = 3.5
)

# scatterplots
lim_all <- max(
  abs(c(
    bio_w$lfc_LBSC,
    bio_w$lfc_DBSC,
    fr_w$lfc_BONCAT,
    fr_w$lfc_Bulk
  )),
  na.rm = TRUE
)

lim_all <- ceiling(lim_all * 10) / 10

# Scatterplot 5A

p_sc1 <- ggplot(bio_w, aes(x = lfc_LBSC, y = lfc_DBSC)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(color = Phylum_plot, alpha = alpha_pt), size = 1.8) +
  scale_alpha_identity() +
  scale_color_manual(values = PAL_BIO, drop = FALSE, name = "Phylum") +
  coord_cartesian(xlim = c(-lim_all, lim_all), ylim = c(-lim_all, lim_all)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  ) +
  ggrepel::geom_text_repel(
    data = bio_lab,
    aes(label = label_plot),
    size = 3,
    max.overlaps = 100,
    na.rm = TRUE
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = ann_bio,
    hjust = -0.05,
    vjust = 1.1,
    size = 3.2,
    label.size = 0,
    alpha = 0.9
  )

# Scatterplot 5B

p_sc2 <- ggplot(fr_w, aes(x = lfc_BONCAT, y = lfc_Bulk)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(color = Phylum_plot, alpha = alpha_pt), size = 1.8) +
  scale_alpha_identity() +
  scale_color_manual(values = PAL_FR, drop = FALSE, name = "Phylum") +
  coord_cartesian(xlim = c(-lim_all, lim_all), ylim = c(-lim_all, lim_all)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  ) +
  ggrepel::geom_text_repel(
    data = fr_lab,
    aes(label = label_plot),
    size = 3,
    max.overlaps = 100,
    na.rm = TRUE
  ) +
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = ann_fr,
    hjust = -0.05,
    vjust = 1.1,
    size = 3.2,
    label.size = 0,
    alpha = 0.9
  )

#combine

fig_sc1 <- cowplot::plot_grid(
  venn_bio,
  p_sc1,
  ncol = 1,
  rel_heights = c(0.45, 1)
)

fig_sc2 <- cowplot::plot_grid(
  venn_frac,
  p_sc2,
  ncol = 1,
  rel_heights = c(0.45, 1)
)

fig_sc1_labeled <- cowplot::ggdraw(fig_sc1) +
  cowplot::draw_label(
    "A",
    x = 0.01,
    y = 0.99,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 14
  )

fig_sc2_labeled <- cowplot::ggdraw(fig_sc2) +
  cowplot::draw_label(
    "B",
    x = 0.01,
    y = 0.99,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 14
  )

fig5 <- cowplot::plot_grid(
  fig_sc1_labeled,
  fig_sc2_labeled,
  nrow = 1,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

fig5 # Fig5 was later modified in .ppt to add concordance arrows, modify axis titles and improve visualization.

