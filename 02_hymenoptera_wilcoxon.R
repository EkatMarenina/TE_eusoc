library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(phytools)

xlsx_path <- "/cloud/project/Hymenoptera sociality.xlsx"

ext  <- read_xlsx(xlsx_path, sheet = "Extended results")
ag   <- read_xlsx(xlsx_path, sheet = "All genomes")
leg  <- read_xlsx(xlsx_path, sheet = "Legend")
tax  <- read_xlsx(xlsx_path, sheet = "Genomes taxonomy")

social_map <- c(
  "YES"                  = "Eusocial",
  "NO"                   = "Solitary",
  "PRIMITIVELY EUSOCIAL" = "Primitive",
  "FACULTATIVE"          = "Partially social",
  "Primitive social"     = "Primitive",
  "Partialy social"      = "Partially social",
  "Klepto"               = "Kleptoparasite"
)

meta <- ext %>%
  left_join(ag  %>% select(`Assembly Accession`, `Clear name`, Social, Family),
            by = c("Name" = "Assembly Accession")) %>%
  left_join(leg %>% select(sp, social_ext),
            by = c("Clear name" = "sp")) %>%
  left_join(tax %>% select(species, family, superfamily, suborder, subfamily, genus),
            by = c("Clear name" = "species")) %>%
  mutate(
    Social_final = coalesce(Social, social_ext),
    Social_final = recode(Social_final, !!!social_map),
    Social_final = case_when(
      is.na(Social_final) & !is.na(`Clear name`) ~ "Solitary",
      is.na(Social_final) &  is.na(`Clear name`) ~ "Unknown",
      TRUE ~ Social_final
    ),
    Social_binary = case_when(
      Social_final %in% c("Eusocial","Primitive","Partially social") ~ "Social",
      Social_final %in% c("Solitary","Kleptoparasite")              ~ "Non-social",
      TRUE ~ NA_character_
    ),
    tip_label  = str_squish(`Clear name`),
    label_safe = str_replace_all(tip_label, " ", "_")
  )

repeat_cols <- c("SINEs","L2/CR1/Rex","R2/R4/NeSL","RTE/Bov-B","L1/CIN4",
                 "Gypsy/DIRS1","Retroviral","hobo-Activator","Tc1-IS630-Pogo",
                 "MULE-MuDR","PiggyBac","Tourist/Harbinger","Rolling-circles",
                 "Small RNA","Satellites","Simple repeats","Low complexity")

tree_data <- meta %>%
  filter(!is.na(tip_label), !is.na(Social_binary)) %>%
  distinct(tip_label, .keep_all = TRUE)

tree_data_tax <- tree_data %>%
  mutate(
    suborder2    = if_else(is.na(suborder),    "Hymenoptera_unk",           suborder),
    superfamily2 = if_else(is.na(superfamily), paste0(suborder2,    "_sf"),  superfamily),
    family2      = if_else(is.na(family),       paste0(superfamily2, "_fam"), family),
    subfamily2   = if_else(is.na(subfamily),    paste0(family2,      "_sub"), subfamily),
    genus2       = if_else(is.na(genus),        paste0(family2,      "_gen"), genus),
    across(c(suborder2, superfamily2, family2, subfamily2, genus2, label_safe), as.factor)
  )

phylo_tree <- as.phylo(
  ~suborder2/superfamily2/family2/subfamily2/genus2/label_safe,
  data     = tree_data_tax,
  collapse = FALSE
)

cat(sprintf("До multi2di: %d tips, %d nodes, binary=%s\n",
    length(phylo_tree$tip.label), phylo_tree$Nnode, is.binary(phylo_tree)))

phylo_tree <- collapse.singles(phylo_tree)
cat(sprintf("После collapse.singles: %d tips, %d nodes\n",
    length(phylo_tree$tip.label), phylo_tree$Nnode))

set.seed(42)
phylo_tree <- multi2di(phylo_tree, random = TRUE)

if (!is.rooted(phylo_tree)) {
  phylo_tree <- root(phylo_tree, outgroup = phylo_tree$tip.label[1], resolve.root = TRUE)
  phylo_tree <- multi2di(phylo_tree, random = TRUE)
}

phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))

cat(sprintf("Итог: %d tips, %d nodes, rooted=%s, binary=%s\n",
    length(phylo_tree$tip.label), phylo_tree$Nnode,
    is.rooted(phylo_tree), is.binary(phylo_tree)))

expected_nodes <- length(phylo_tree$tip.label) - 1L
if (phylo_tree$Nnode != expected_nodes) {
  stop(sprintf("Всё ещё проблема: %d nodes, ожидалось %d", phylo_tree$Nnode, expected_nodes))
}
message("✓ Дерево готово для ACE")

tip_social <- setNames(
  tree_data$Social_binary,
  tree_data$label_safe
)
tip_social <- tip_social[phylo_tree$tip.label]

ace_result <- ace(
  tip_social,
  phylo_tree,
  type  = "discrete",
  model = "ER"
)

node_probs <- ace_result$lik.anc
colnames(node_probs) <- colnames(ace_result$lik.anc)

node_state <- apply(node_probs, 1, function(x) colnames(node_probs)[which.max(x)])

n_tips  <- length(phylo_tree$tip.label)
n_nodes <- phylo_tree$Nnode
node_ids <- (n_tips + 1):(n_tips + n_nodes)

transition_nodes <- c()
transition_type  <- c()

for (i in seq_along(node_ids)) {
  node <- node_ids[i]
  parent <- which(phylo_tree$edge[, 2] == node)
  if (length(parent) == 0) next

  parent_node <- phylo_tree$edge[parent, 1]
  parent_state <- if (parent_node > n_tips) {
    node_state[parent_node - n_tips]
  } else {
    tip_social[parent_node]
  }
  child_state <- node_state[i]

  if (!is.na(parent_state) && !is.na(child_state) && parent_state != child_state) {
    transition_nodes <- c(transition_nodes, node)
    transition_type  <- c(transition_type,
                          paste0(parent_state, " → ", child_state))
  }
}

transition_df <- data.frame(
  node       = transition_nodes,
  transition = transition_type
)

cat("=== Найденные узлы переходов ===\n")
print(transition_df)
cat(sprintf("\nВсего переходов: %d\n", nrow(transition_df)))
cat(sprintf("  Non-social → Social: %d\n",
            sum(grepl("Non-social → Social", transition_df$transition))))
cat(sprintf("  Social → Non-social: %d\n",
            sum(grepl("Social → Non-social", transition_df$transition))))

wilcox_results <- lapply(repeat_cols, function(col) {
  social_vals     <- tree_data %>%
    filter(Social_binary == "Social")     %>% pull(!!sym(col)) %>% as.numeric()
  nonsocial_vals  <- tree_data %>%
    filter(Social_binary == "Non-social") %>% pull(!!sym(col)) %>% as.numeric()

  wt <- wilcox.test(social_vals, nonsocial_vals, exact = FALSE)

  data.frame(
    TE          = col,
    median_social    = median(social_vals,    na.rm = TRUE),
    median_nonsocial = median(nonsocial_vals, na.rm = TRUE),
    W           = wt$statistic,
    p_value     = wt$p.value,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(
    p_adj      = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05,
    direction  = case_when(
      median_social > median_nonsocial ~ "Higher in Social",
      median_social < median_nonsocial ~ "Higher in Non-social",
      TRUE ~ "Equal"
    )
  ) %>%
  arrange(p_adj)

cat("\n=== Результаты теста Уилкоксона (поправка BH) ===\n")
print(wilcox_results %>%
        select(TE, median_social, median_nonsocial, W, p_value, p_adj, significant, direction),
      digits = 4, row.names = FALSE)

write.csv(wilcox_results, "wilcoxon_TE_social_vs_nonsocial.csv", row.names = FALSE)

social_colors_bin <- c("Social" = "

plot_data <- tree_data %>%
  select(tip_label, Social_binary, all_of(repeat_cols)) %>%
  pivot_longer(cols = all_of(repeat_cols), names_to = "TE", values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%
  left_join(wilcox_results %>% select(TE, p_adj, significant, direction),
            by = "TE") %>%
  mutate(
    TE_label = if_else(significant,
                       paste0(TE, "*"),
                       TE),
    TE = factor(TE, levels = wilcox_results$TE)
  )

p_box <- ggplot(plot_data, aes(x = Social_binary, y = value, fill = Social_binary)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.5, linewidth = 0.4) +
  geom_jitter(aes(color = Social_binary), width = 0.15, size = 0.4, alpha = 0.4) +
  facet_wrap(~ TE, scales = "free_y", ncol = 6) +
  scale_fill_manual(values  = social_colors_bin, name = NULL) +
  scale_color_manual(values = social_colors_bin, guide = "none") +
  geom_text(
    data = plot_data %>%
      filter(significant) %>%
      group_by(TE) %>%
      summarise(y_pos = max(value, na.rm = TRUE) * 1.1, .groups = "drop") %>%
      mutate(label = "*", Social_binary = "Social"),
    aes(y = y_pos, label = label),
    size = 5, color = "black", vjust = 0
  ) +
  labs(
    title    = "Transposable element content: Social vs Non-social Hymenoptera",
    subtitle = "Wilcoxon rank-sum test, BH correction; * p_adj < 0.05; ordered by p_adj",
    x        = NULL,
    y        = "% of genome"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "top",
    strip.text       = element_text(size = 8, face = "bold"),
    axis.text.x      = element_text(size = 8, angle = 20, hjust = 1),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave("wilcoxon_boxplots.pdf",
       plot = p_box, width = 18, height = 12, device = cairo_pdf)
ggsave("wilcoxon_boxplots.png",
       plot = p_box, width = 18, height = 12, dpi = 300)

tip_meta_bin <- tree_data %>%
  select(label = label_safe, tip_label, Social_binary, Social_final, family)

node_data <- data.frame(
  node         = node_ids,
  Social_binary = node_state,
  stringsAsFactors = FALSE
)

transition_df$node_label <- ifelse(
  grepl("Non-social → Social", transition_df$transition), "↑Social", "↓Non-social"
)

p_tree <- ggtree(phylo_tree, branch.length = "none",
                 layout = "rectangular", color = "grey40", size = 0.2) %<+%
  tip_meta_bin +
  geom_tippoint(aes(color = Social_binary), size = 1.8, alpha = 0.9) +
  scale_color_manual(
    values   = social_colors_bin,
    name     = "Sociality",
    na.value = "
    guide    = guide_legend(override.aes = list(size = 4))
  ) +
  geom_point2(
    aes(subset = (node %in% transition_df$node[
                    grepl("Non-social → Social", transition_df$transition)])),
    color = "
  ) +
  geom_point2(
    aes(subset = (node %in% transition_df$node[
                    grepl("Social → Non-social", transition_df$transition)])),
    color = "
  ) +
  geom_tiplab(
    aes(label = tip_label, color = Social_binary),
    size = 1.8, align = TRUE, linesize = 0.2,
    offset = 0.3, fontface = "italic", show.legend = FALSE
  ) +
  scale_color_manual(values = social_colors_bin, na.value = "
  theme_tree2() +
  theme(
    legend.position = c(0.02, 0.92),
    legend.background = element_rect(fill = "white", color = "grey70"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text  = element_text(size = 9),
    plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.margin  = margin(5, 120, 5, 5)
  ) +
  labs(title = "Hymenoptera: inferred social transitions (ACE, ER model)\n◆ yellow = Non-social→Social,  ◆ orange = Social→Non-social")

ggsave("transition_tree.pdf",
       plot = p_tree, width = 16, height = 22, device = cairo_pdf)
ggsave("transition_tree.png",
       plot = p_tree, width = 16, height = 22, dpi = 300)

message("✓ Готово!")
message("  wilcoxon_TE_social_vs_nonsocial.csv")
message("  wilcoxon_boxplots.pdf / .png")
message("  transition_tree.pdf / .png")
