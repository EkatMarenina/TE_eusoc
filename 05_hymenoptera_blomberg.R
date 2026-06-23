library(phytools)
library(dplyr)
library(ggplot2)
library(tidyr)

cat("=== Blomberg's K (phylosignal) ===\n")
cat("Пермутаций: 1000, метод: Blomberg et al. 2003\n\n")

blomberg_results <- lapply(te_variable, function(te) {
  vals <- setNames(
    log1p(as.numeric(tree_data[[te]])),
    tree_data$label_safe
  )
  vals <- vals[!is.na(vals)]
  vals <- vals[names(vals) %in% phylo_tree$tip.label]

  tryCatch({
    res <- phylosig(
      tree   = phylo_tree,
      x      = vals,
      method = "K",
      test   = TRUE,
      nsim   = 1000
    )
    data.frame(
      TE          = te,
      K           = round(res$K, 4),
      p_value     = res$P,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  Ошибка %s: %s\n", te, e$message))
    data.frame(TE = te, K = NA, p_value = NA, stringsAsFactors = FALSE)
  })
}) %>% bind_rows() %>%
  filter(!is.na(K)) %>%
  mutate(
    p_adj       = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05,
    signal      = case_when(
      !significant    ~ "No signal",
      K >= 1          ~ "Strong (K≥1)",
      K >= 0.5        ~ "Moderate (K≥0.5)",
      TRUE            ~ "Weak (K<0.5)"
    )
  ) %>%
  arrange(p_adj)

cat("\nРезультаты Blomberg's K:\n")
print(blomberg_results, digits = 4, row.names = FALSE)

write.csv(blomberg_results, "blomberg_K_results.csv", row.names = FALSE)

summary_table <- blomberg_results %>%
  dplyr::select(TE, K, p_adj_K = p_adj, sig_K = significant) %>%
  left_join(
    pgls_results %>% dplyr::select(TE, lambda, p_adj_pgls = p_adj, sig_pgls = significant),
    by = "TE"
  ) %>%
  left_join(
    wilcox_results %>% dplyr::select(TE, p_adj_wilcox = p_adj, sig_wilcox = significant, direction),
    by = "TE"
  ) %>%
  mutate(
    interpretation = case_when(
      sig_K & sig_wilcox & !sig_pgls ~
        "Wilcoxon significant but driven by phylogeny",
      sig_K & sig_pgls ~
        "Robust signal (survives phylo correction)",
      !sig_K & sig_wilcox ~
        "Wilcoxon significant, low phylo signal",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(p_adj_K)

cat("\n=== Сводная таблица: K + lambda + Wilcoxon ===\n")
print(summary_table, row.names = FALSE)

write.csv(summary_table, "phylosignal_summary.csv", row.names = FALSE)

signal_colors <- c(
  "Strong (K≥1)"     = "
  "Moderate (K≥0.5)" = "
  "Weak (K<0.5)"     = "
  "No signal"        = "
)

p_K <- blomberg_results %>%
  mutate(
    TE        = reorder(TE, K),
    neg_log_p = -log10(p_adj)
  ) %>%
  ggplot(aes(x = K, y = neg_log_p, color = signal, label = TE)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3.2, max.overlaps = 20) +
  scale_color_manual(values = signal_colors, name = "Signal strength") +
  annotate("text", x = 1.02, y = 0.1,
           label = "K = 1\n(Brownian)", size = 3, color = "grey40", hjust = 0) +
  annotate("text", x = 0.05, y = -log10(0.05) + 0.1,
           label = "p_adj = 0.05", size = 3, color = "grey40") +
  labs(
    title    = "Phylogenetic signal in TE content (Blomberg's K)",
    subtitle = "Points above dashed horizontal line: significant signal (p_adj < 0.05, BH)\nPoints right of K=1: stronger signal than Brownian motion",
    x        = "Blomberg's K",
    y        = "-log10(p_adj)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 9, color = "grey40")
  )

te_order <- summary_table %>%
  arrange(p_adj_K) %>%
  pull(TE)

heat_methods <- summary_table %>%
  dplyr::select(TE, K, lambda, p_adj_K, p_adj_pgls, p_adj_wilcox) %>%
  pivot_longer(cols = starts_with("p_adj"),
               names_to = "method", values_to = "p_adj") %>%
  mutate(
    method    = recode(method,
                       "p_adj_K"      = "Blomberg K",
                       "p_adj_pgls"   = "pGLS",
                       "p_adj_wilcox" = "Wilcoxon"),
    neg_log_p = pmin(-log10(p_adj), 10),
    sig       = p_adj < 0.05,
    TE        = factor(TE, levels = te_order)
  )

p_heat_methods <- ggplot(heat_methods,
                         aes(x = method, y = TE, fill = neg_log_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = if_else(sig, "*", "")),
            size = 5, color = "white", vjust = 0.7) +
  scale_fill_gradient(low = "
                      name = "-log10(p_adj)",
                      na.value = "grey90") +
  labs(
    title    = "Significance across methods (* p_adj < 0.05)",
    subtitle = "Ordered by Blomberg's K p_adj",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(face = "bold", size = 10),
    plot.title   = element_text(face = "bold"),
    panel.grid   = element_blank()
  )

library(cowplot)
final <- plot_grid(p_K, p_heat_methods,
                   ncol = 2, labels = c("A", "B"),
                   rel_widths = c(1.4, 1))

ggsave("blomberg_phylosignal.pdf",
       plot = final, width = 16, height = 8, device = cairo_pdf)
ggsave("blomberg_phylosignal.png",
       plot = final, width = 16, height = 8, dpi = 300)

message("Готово!")
message("  blomberg_K_results.csv")
message("  phylosignal_summary.csv")
message("  blomberg_phylosignal.pdf / .png")
