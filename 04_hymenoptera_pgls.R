library(caper)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

te_variable <- wilcox_results %>%
  filter(!(median_social == 0 & median_nonsocial == 0)) %>%
  filter(!is.na(p_adj)) %>%
  pull(TE)

cat("TE включены в анализ:", paste(te_variable, collapse = ", "), "\n\n")

pgls_data <- tree_data %>%
  mutate(
    Social_num = if_else(Social_binary == "Social", 1L, 0L),
    across(all_of(te_variable), ~ log1p(as.numeric(.x)))
  ) %>%
  dplyr::select(label_safe, Social_binary, Social_num, all_of(te_variable)) %>%
  as.data.frame()

phylo_tree$node.label <- NULL

comp_data <- comparative.data(
  phy          = phylo_tree,
  data         = pgls_data,
  names.col    = "label_safe",
  vcv          = TRUE,
  warn.dropped = TRUE
)

cat(sprintf("Видов в анализе: %d\n\n", nrow(comp_data$data)))

cat("=== 3.2 phylo.d ===\n")

d_result <- phylo.d(
  data      = comp_data,
  binvar    = Social_num,
  permut    = 1000
)

print(d_result)

d_value <- d_result$DEstimate
p_random <- d_result$Pval1
p_brownian <- d_result$Pval0

cat(sprintf("\nD = %.4f\n", d_value))
cat(sprintf("p (vs random/D=1):    %.4f\n", p_random))
cat(sprintf("p (vs Brownian/D=0):  %.4f\n", p_brownian))

if (d_value < 0) {
  cat("Интерпретация: признак сильно консервативен (клинальный)\n")
} else if (d_value < 0.5) {
  cat("Интерпретация: признак умеренно консервативен\n")
} else if (d_value < 1) {
  cat("Интерпретация: признак слабо кластеризован\n")
} else {
  cat("Интерпретация: признак близок к случайному распределению\n")
}

cat("\n=== 3.3 pGLS (Pagel's lambda) ===\n")
cat("Модель: log1p(TE) ~ Social_binary\n")
cat("Ковариационная структура: Pagel's lambda (ML)\n\n")

pgls_results <- lapply(te_variable, function(te) {
  formula_str <- paste0("`", te, "` ~ Social_binary")
  
  tryCatch({
    mod <- pgls(
      as.formula(formula_str),
      data   = comp_data,
      lambda = "ML"
    )
    s <- summary(mod)
    coef_social <- coef(mod)["Social_binarySocial"]
    p_val       <- s$coefficients["Social_binarySocial", "Pr(>|t|)"]
    lambda_est  <- mod$param["lambda"]
    r2          <- s$r.squared
    
    data.frame(
      TE         = te,
      beta       = round(coef_social, 4),
      lambda     = round(lambda_est,  4),
      R2         = round(r2,          4),
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  Ошибка для %s: %s\n", te, e$message))
    data.frame(TE = te, beta = NA, lambda = NA, R2 = NA, p_value = NA,
               stringsAsFactors = FALSE)
  })
}) %>% bind_rows() %>%
  filter(!is.na(p_value)) %>%
  mutate(
    p_adj       = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05,
    direction   = case_when(
      beta > 0  ~ "Higher in Social",
      beta < 0  ~ "Higher in Non-social",
      TRUE      ~ "Equal"
    )
  ) %>%
  arrange(p_adj)

cat("\nРезультаты pGLS:\n")
print(pgls_results %>%
        dplyr::select(TE, beta, lambda, R2, p_value, p_adj, significant, direction),
      digits = 4, row.names = FALSE)

write.csv(pgls_results, "pgls_results.csv", row.names = FALSE)

cat("\n=== 3.1 BRUNCH ===\n")
cat("Метод: контрасты между кладами разной социальности\n")
cat("Каждый контраст = точка перехода Social <-> Non-social\n\n")

brunch_results <- lapply(te_variable, function(te) {
  formula_str <- paste0("`", te, "` ~ Social_binary")
  
  tryCatch({
    mod <- brunch(
      as.formula(formula_str),
      data = comp_data
    )
    s <- summary(mod)
    t_val <- s$coefficients[1, "t value"]
    p_val <- s$coefficients[1, "Pr(>|t|)"]
    n_contrasts <- nrow(mod$contrasts)
    
    data.frame(
      TE           = te,
      n_contrasts  = n_contrasts,
      t_value      = round(t_val, 4),
      p_value      = p_val,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  Ошибка для %s: %s\n", te, e$message))
    data.frame(TE = te, n_contrasts = NA, t_value = NA, p_value = NA,
               stringsAsFactors = FALSE)
  })
}) %>% bind_rows() %>%
  filter(!is.na(p_value)) %>%
  mutate(
    p_adj       = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05
  ) %>%
  arrange(p_adj)

cat("\nРезультаты BRUNCH:\n")
print(brunch_results, digits = 4, row.names = FALSE)

write.csv(brunch_results, "brunch_results.csv", row.names = FALSE)

p_lambda <- pgls_results %>%
  filter(!is.na(lambda)) %>%
  mutate(TE = reorder(TE, lambda)) %>%
  ggplot(aes(x = TE, y = lambda, fill = significant)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  scale_fill_manual(values = c("TRUE" = "
                    name = "pGLS significant\n(p_adj < 0.05)") +
  coord_flip() +
  labs(
    title    = "Pagel's lambda per TE class (pGLS)",
    subtitle = "lambda=1: strong phylogenetic signal; lambda=0: no signal",
    x = NULL, y = "lambda"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

compare_df <- wilcox_results %>%
  dplyr::select(TE, p_adj_wilcox = p_adj) %>%
  left_join(pgls_results   %>% dplyr::select(TE, p_adj_pgls   = p_adj), by = "TE") %>%
  left_join(brunch_results %>% dplyr::select(TE, p_adj_brunch = p_adj), by = "TE") %>%
  filter(!is.na(p_adj_pgls)) %>%
  pivot_longer(cols = starts_with("p_adj"),
               names_to = "method", values_to = "p_adj") %>%
  mutate(
    method    = recode(method,
                       "p_adj_wilcox" = "Wilcoxon",
                       "p_adj_pgls"   = "pGLS",
                       "p_adj_brunch" = "BRUNCH"),
    neg_log_p = -log10(p_adj),
    TE        = reorder(TE, neg_log_p)
  )

p_compare <- ggplot(compare_df,
                    aes(x = TE, y = neg_log_p, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  annotate("text", x = 1, y = -log10(0.05) + 0.1,
           label = "p_adj = 0.05", size = 3, color = "grey30", hjust = 0) +
  scale_fill_manual(values = c("Wilcoxon" = "
                               "pGLS"     = "
                               "BRUNCH"   = "
                    name = "Method") +
  coord_flip() +
  labs(
    title    = "Comparison of methods: -log10(p_adj) per TE class",
    subtitle = "Dashed line = significance threshold (p_adj = 0.05)",
    x = NULL, y = "-log10(p_adj)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

d_plot_df <- data.frame(
  label = c("Observed D", "Random (D=1)", "Brownian (D=0)"),
  value = c(d_value, 1, 0),
  type  = c("observed", "reference", "reference")
)

p_d <- ggplot(d_plot_df, aes(x = value, y = 1, color = type, size = type)) +
  geom_segment(aes(x = 0, xend = 1, y = 1, yend = 1),
               color = "grey80", linewidth = 1) +
  geom_point(shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = c("observed"  = "
                                "reference" = "grey50")) +
  scale_size_manual(values  = c("observed"  = 5,
                                "reference" = 3)) +
  scale_x_continuous(limits = c(-0.5, 1.5),
                     breaks = c(0, d_value, 1),
                     labels = c("0\n(Brownian)",
                                sprintf("%.2f\n(Observed)", d_value),
                                "1\n(Random)")) +
  labs(
    title    = paste0("phylo.d: D = ", round(d_value, 3)),
    subtitle = sprintf("p(vs random) = %.4f  |  p(vs Brownian) = %.4f",
                       p_random, p_brownian),
    x = "D statistic", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold")
  )

final <- plot_grid(
  p_d,
  plot_grid(p_lambda, p_compare, ncol = 2, labels = c("B", "C")),
  ncol   = 1,
  rel_heights = c(0.25, 0.75),
  labels = c("A", "")
)

ggsave("pgls_summary.pdf",
       plot = final, width = 16, height = 14, device = cairo_pdf)
ggsave("pgls_summary.png",
       plot = final, width = 16, height = 14, dpi = 300)

message("Готово!")
message("  pgls_results.csv")
message("  brunch_results.csv")
message("  pgls_summary.pdf / .png")
