library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(gridExtra)
library(grid)

social_colors_bin <- c("Social" = "

repeat_cols <- c("SINEs","L2/CR1/Rex","R2/R4/NeSL","RTE/Bov-B","L1/CIN4",
                 "Gypsy/DIRS1","Retroviral","hobo-Activator","Tc1-IS630-Pogo",
                 "MULE-MuDR","PiggyBac","Tourist/Harbinger","Rolling-circles",
                 "Small RNA","Satellites","Simple repeats","Low complexity")

sig_TE <- wilcox_results %>% filter(significant) %>% pull(TE)

te_labels <- wilcox_results %>%
  filter(significant) %>%
  mutate(
    arrow      = if_else(direction == "Higher in Social", "higher in Social",
                         if_else(direction == "Higher in Non-social", "higher in Non-social",
                                 "no direction")),
    label_full = paste0(TE, "\n(", arrow, ", p_adj=",
                        formatC(p_adj, format = "e", digits = 1), ")")
  ) %>%
  select(TE, label_full, direction, p_adj)

plot_data_sig <- tree_data %>%
  select(tip_label, Social_binary, all_of(sig_TE)) %>%
  pivot_longer(cols = all_of(sig_TE), names_to = "TE", values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%
  left_join(te_labels, by = "TE") %>%
  mutate(
    TE_factor     = factor(label_full,
                           levels = te_labels %>% arrange(p_adj) %>% pull(label_full)),
    Social_binary = factor(Social_binary, levels = c("Social", "Non-social"))
  )

p_box_sig <- ggplot(plot_data_sig,
                    aes(x = Social_binary, y = value, fill = Social_binary)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.2,
               outlier.alpha = 0.6, linewidth = 0.45, width = 0.55) +
  geom_jitter(aes(color = Social_binary),
              width = 0.15, size = 0.7, alpha = 0.35) +
  facet_wrap(~ TE_factor, scales = "free_y", ncol = 4) +
  scale_fill_manual(values  = social_colors_bin, name = NULL) +
  scale_color_manual(values = social_colors_bin, guide = "none") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(
    title    = "Significant differences in transposable element content",
    subtitle = paste0(
      "Wilcoxon rank-sum test, Benjamini-Hochberg correction  |  ",
      "n(Social)=", sum(tree_data$Social_binary == "Social"),
      ", n(Non-social)=", sum(tree_data$Social_binary == "Non-social")
    ),
    x = NULL,
    y = "% of genome"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position    = "top",
    legend.key.size    = unit(0.9, "lines"),
    strip.text         = element_text(size = 8.5, face = "bold"),
    strip.background   = element_rect(fill = "grey95", color = "grey70"),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave("wilcoxon_boxplots_significant.pdf",
       plot = p_box_sig, width = 14, height = 10, device = cairo_pdf)
ggsave("wilcoxon_boxplots_significant.png",
       plot = p_box_sig, width = 14, height = 10, dpi = 300)

wilcox_table <- wilcox_results %>%
  mutate(
    p_value   = formatC(p_value, format = "e", digits = 2),
    p_adj_fmt = formatC(p_adj,   format = "e", digits = 2),
    sig_label = case_when(
      significant & as.numeric(p_adj) < 0.001 ~ "***",
      significant & as.numeric(p_adj) < 0.01  ~ "**",
      significant                              ~ "*",
      TRUE                                     ~ "ns"
    ),
    dir_label = case_when(
      direction == "Higher in Social"     ~ "up Social",
      direction == "Higher in Non-social" ~ "up Non-social",
      TRUE                                ~ "-"
    ),
    median_social    = round(median_social,    3),
    median_nonsocial = round(median_nonsocial, 3)
  ) %>%
  select(TE, median_social, median_nonsocial, p_value, p_adj_fmt, sig_label, dir_label) %>%
  rename(
    "TE class"            = TE,
    "Median Social"       = median_social,
    "Median Non-social"   = median_nonsocial,
    "p-value"             = p_value,
    "p_adj (BH)"          = p_adj_fmt,
    "Sig."                = sig_label,
    "Direction"           = dir_label
  )

tbl_grob <- tableGrob(
  wilcox_table,
  rows  = NULL,
  theme = ttheme_minimal(
    core    = list(
      fg_params = list(fontsize = 9),
      bg_params = list(fill = rep(c("white", "grey96"), length.out = nrow(wilcox_table)))
    ),
    colhead = list(
      fg_params = list(fontsize = 10, fontface = "bold"),
      bg_params = list(fill = "grey85")
    )
  )
)

pdf("wilcoxon_table.pdf", width = 12, height = 8)
grid.draw(tbl_grob)
dev.off()

message("wilcoxon_boxplots_significant.pdf / .png")
message("wilcoxon_table.pdf")
