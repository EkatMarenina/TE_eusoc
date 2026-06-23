library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(cowplot)

xlsx_path <- "/cloud/project/Hymenoptera sociality.xlsx"

ext  <- read_xlsx(xlsx_path, sheet = "Extended results")
ag   <- read_xlsx(xlsx_path, sheet = "All genomes")
leg  <- read_xlsx(xlsx_path, sheet = "Legend")
tax  <- read_xlsx(xlsx_path, sheet = "Genomes taxonomy")

meta <- ext %>%
  left_join(ag  %>% select(`Assembly Accession`, `Clear name`, Social, Family),
            by = c("Name" = "Assembly Accession")) %>%
  left_join(leg %>% select(sp, social_ext),
            by = c("Clear name" = "sp")) %>%
  left_join(tax %>% select(species, family, superfamily, suborder, subfamily, genus),
            by = c("Clear name" = "species"))

social_map <- c(
  "YES"                  = "Eusocial",
  "NO"                   = "Solitary",
  "PRIMITIVELY EUSOCIAL" = "Primitive",
  "FACULTATIVE"          = "Partially social",
  "Primitive social"     = "Primitive",
  "Partialy social"      = "Partially social",
  "Klepto"               = "Kleptoparasite"
)

meta <- meta %>%
  mutate(
    Social_final = coalesce(Social, social_ext),
    Social_final = recode(Social_final, !!!social_map),
    Social_final = case_when(
      is.na(Social_final) & !is.na(`Clear name`) ~ "Solitary",
      is.na(Social_final) &  is.na(`Clear name`) ~ "Unknown",
      TRUE ~ Social_final
    ),
    tip_label  = str_squish(`Clear name`),
    label_safe = str_replace_all(tip_label, " ", "_")
  )

tree_data <- meta %>%
  filter(!is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE)

repeat_cols <- c("SINEs","L2/CR1/Rex","R2/R4/NeSL","RTE/Bov-B","L1/CIN4",
                 "Gypsy/DIRS1","Retroviral","hobo-Activator","Tc1-IS630-Pogo",
                 "MULE-MuDR","PiggyBac","Tourist/Harbinger","Rolling-circles",
                 "Unclassified","Small RNA","Satellites","Simple repeats","Low complexity")

heat_mat <- tree_data %>%
  select(label_safe, any_of(repeat_cols)) %>%
  tibble::column_to_rownames("label_safe") %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(Unclassified = 0) %>%
  mutate(across(everything(), ~ {
    mn <- min(.x, na.rm = TRUE)
    mx <- max(.x, na.rm = TRUE)
    if (mx > mn) (.x - mn) / (mx - mn) else rep(0, length(.x))
  }))

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

tip_meta <- tree_data %>%
  select(label = label_safe, tip_label, Social_final, family, suborder)

social_colors <- c(
  "Eusocial"         = "
  "Primitive"        = "
  "Partially social" = "
  "Kleptoparasite"   = "
  "Solitary"         = "
  "Unknown"          = "
)

base_tree <- ggtree(phylo_tree, branch.length = "none",
                    layout = "rectangular", color = "grey30", size = 0.25) %<+%
  tip_meta +
  geom_tippoint(aes(color = Social_final), size = 2, alpha = 0.95) +
  scale_color_manual(
    values   = social_colors,
    name     = "Sociality",
    na.value = "
    guide    = guide_legend(
      title.position = "top",
      override.aes   = list(size = 6),
      keyheight      = unit(1.4, "lines")
    )
  ) +
  xlim(0, 25)

heat_offset <- 4
heat_width  <- 1.5
n_col       <- ncol(heat_mat)
max_x_tree  <- max(base_tree$data$x, na.rm = TRUE)
col_step    <- (max_x_tree * heat_width) / n_col
col_start   <- max_x_tree + heat_offset + col_step * 0.5
col_xs      <- col_start + col_step * seq(0, n_col - 1)

y_bottom    <- 0.5
tick_len    <- 0.8

hline_df <- data.frame(
  x    = c(col_xs[1] - col_step * 0.5, col_xs[n_col] + col_step * 0.5),
  y    = y_bottom,
  type = "line"
)

ticks_df <- data.frame(
  x     = col_xs,
  y_top = y_bottom,
  y_bot = y_bottom - tick_len
)

labels_df <- data.frame(
  x     = col_xs,
  y     = y_bottom - tick_len - 0.2,
  label = colnames(heat_mat)
)

p_heat <- gheatmap(
  base_tree,
  heat_mat,
  offset   = heat_offset,
  width    = heat_width,
  font.size = 0,
  colnames  = FALSE
) +
  scale_fill_gradient(
    low      = "lightgreen",
    high     = "tomato",
    na.value = "transparent",
    name     = "Relative\ncontent"
  ) +
  geom_segment(
    data        = hline_df,
    aes(x = x[1], xend = x[2], y = y_bottom, yend = y_bottom),
    color       = "grey20",
    linewidth   = 0.5,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data        = ticks_df,
    aes(x = x, xend = x, y = y_top, yend = y_bot),
    color       = "grey20",
    linewidth   = 0.4,
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = labels_df,
    aes(x = x, y = y, label = label),
    angle       = -30,
    hjust       = 1,
    vjust       = 0.5,
    size        = 2.8,
    inherit.aes = FALSE
  ) +
  new_scale_color() +
  geom_tiplab(
    aes(label = tip_label, color = Social_final),
    size        = 2.2,
    align       = TRUE,
    linesize    = 0.4,
    linetype    = "dotted",
    offset      = 2,
    fontface    = "italic",
    show.legend = FALSE
  ) +
  scale_color_manual(values = social_colors, na.value = "
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 1), limits = c(-10, NA)) +
  coord_cartesian(clip = "off") +
  theme_tree2() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.line.x      = element_blank(),
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.margin      = margin(5, 5, 70, 5)
  ) +
  labs(title = "Taxonomic tree of Hymenoptera (n = 186 genomes) with repeat element content")

leg_social <- get_legend(
  base_tree +
    theme(
      legend.position      = "right",
      legend.justification = "center",
      legend.background    = element_rect(fill = "white", color = "grey70", linewidth = 0.6),
      legend.title         = element_text(face = "bold", size = 13),
      legend.text          = element_text(size = 11),
      legend.margin        = margin(10, 14, 10, 14),
      legend.key.size      = unit(1.2, "lines")
    )
)

leg_fill <- get_legend(
  p_heat +
    theme(
      legend.position      = "right",
      legend.justification = "center",
      legend.background    = element_rect(fill = "white", color = "grey70", linewidth = 0.6),
      legend.title         = element_text(face = "bold", size = 13),
      legend.text          = element_text(size = 11),
      legend.margin        = margin(16, 20, 16, 20),
      legend.key.height    = unit(2.5, "cm"),
      legend.key.width     = unit(0.9, "cm")
    )
)

right_col <- plot_grid(
  leg_social,
  leg_fill,
  ncol        = 1,
  align       = "v",
  rel_heights = c(1, 1)
)

final_plot <- plot_grid(
  p_heat,
  right_col,
  ncol       = 2,
  rel_widths = c(6.5, 0.8)
)

ggsave("hymenoptera_tree_heatmap.pdf",
       plot = final_plot, width = 30, height = 14, device = cairo_pdf)
ggsave("hymenoptera_tree_heatmap.png",
       plot = final_plot, width = 30, height = 14, dpi = 300)

message("✓ Готово: hymenoptera_tree_heatmap.pdf / .png")
