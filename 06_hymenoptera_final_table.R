library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
library(grid)

final_table <- summary_table %>%
  mutate(
    sig_W = case_when(
      p_adj_wilcox < 0.001 ~ "***",
      p_adj_wilcox < 0.01  ~ "**",
      p_adj_wilcox < 0.05  ~ "*",
      TRUE                  ~ "ns"
    ),
    sig_pgls_sym = if_else(sig_pgls, "*", "ns"),
    sig_K_sym    = if_else(sig_K,    "*", "ns"),
    dir_arrow = case_when(
      direction == "Higher in Social"     ~ "↑ Social",
      direction == "Higher in Non-social" ~ "↑ Non-social",
      TRUE                                ~ "—"
    ),
    K      = round(K,      3),
    lambda = round(lambda, 3),
    p_adj_wilcox = formatC(p_adj_wilcox, format = "e", digits = 2),
    p_adj_pgls   = formatC(p_adj_pgls,   format = "e", digits = 2),
    p_adj_K      = formatC(p_adj_K,      format = "e", digits = 2)
  ) %>%
  dplyr::select(
    "TE class"      = TE,
    "Direction"     = dir_arrow,
    "Wilcoxon p_adj"= p_adj_wilcox,
    "W sig."        = sig_W,
    "pGLS p_adj"    = p_adj_pgls,
    "pGLS sig."     = sig_pgls_sym,
    "lambda (pGLS)" = lambda,
    "Blomberg K"    = K,
    "K p_adj"       = p_adj_K,
    "K sig."        = sig_K_sym
  )

write.csv(final_table, "final_summary_table.csv", row.names = FALSE)

bubble_data <- summary_table %>%
  mutate(
    neg_log_wilcox = pmin(-log10(p_adj_wilcox), 12),
    K_size         = pmax(K, 0.05),
    dir_color      = case_when(
      direction == "Higher in Social"     ~ "Higher in Social",
      direction == "Higher in Non-social" ~ "Higher in Non-social",
      TRUE                                ~ "Equal"
    ),
    sig_wilcox_label = if_else(sig_wilcox, TE, "")
  )

dir_colors <- c(
  "Higher in Social"     = "
  "Higher in Non-social" = "
  "Equal"                = "
)

p_bubble <- ggplot(bubble_data,
                   aes(x = neg_log_wilcox, y = lambda,
                       size = K_size, color = dir_color,
                       label = sig_wilcox_label)) +
  annotate("rect", xmin = -log10(0.05), xmax = 13,
           ymin = 0.5, ymax = 1.05,
           fill = "
  annotate("text", x = 11, y = 0.95,
           label = "Wilcoxon sig.\nhigh lambda\n→ phylogenetic\nconfounding",
           size = 2.8, color = "
  annotate("rect", xmin = -log10(0.05), xmax = 13,
           ymin = -0.05, ymax = 0.5,
           fill = "
  annotate("text", x = 11, y = 0.25,
           label = "Wilcoxon sig.\nlow lambda\n→ real signal",
           size = 2.8, color = "
  geom_vline(xintercept = -log10(0.05),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted", color = "grey50", linewidth = 0.5) +
  geom_point(alpha = 0.85) +
  ggrepel::geom_text_repel(
    size = 3.2, max.overlaps = 20,
    box.padding = 0.4, show.legend = FALSE
  ) +
  scale_color_manual(values = dir_colors, name = "Direction") +
  scale_size_continuous(range = c(3, 10), name = "Blomberg's K") +
  scale_x_continuous(
    name   = "-log10(Wilcoxon p_adj)",
    limits = c(0, 13),
    breaks = c(0, -log10(0.05), 2, 4, 6, 8, 10),
    labels = c("0", "p=0.05", "2", "4", "6", "8", "10")
  ) +
  scale_y_continuous(
    name   = "Pagel's lambda (pGLS)",
    limits = c(-0.05, 1.05)
  ) +
  labs(
    title    = "TE content: Wilcoxon significance vs phylogenetic signal strength",
    subtitle = paste0(
      "Point size = Blomberg's K  |  ",
      "Red zone: significant Wilcoxon but high lambda (phylogenetic confounding)\n",
      "Blue zone: significant Wilcoxon with low lambda (potentially real signal)"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 8.5, color = "grey40"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("final_bubble_plot.pdf",
       plot = p_bubble, width = 12, height = 8, device = cairo_pdf)
ggsave("final_bubble_plot.png",
       plot = p_bubble, width = 12, height = 8, dpi = 300)

tbl_grob <- tableGrob(
  final_table,
  rows  = NULL,
  theme = ttheme_minimal(
    core    = list(
      fg_params = list(fontsize = 8),
      bg_params = list(
        fill = rep(c("white","grey96"), length.out = nrow(final_table))
      )
    ),
    colhead = list(
      fg_params = list(fontsize = 9, fontface = "bold"),
      bg_params = list(fill = "grey85")
    )
  )
)

pdf("final_summary_table.pdf", width = 16, height = 6)
grid.draw(tbl_grob)
dev.off()

message("Готово!")
message("  final_summary_table.csv / .pdf")
message("  final_bubble_plot.pdf / .png")
