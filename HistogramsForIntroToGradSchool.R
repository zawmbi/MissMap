library(tidyverse)
setwd("/Users/Linda/Coding/MissMap/")
# 👉 EDIT THESE to your actual file names
baseline_path <- "missmap_results/data_availability_matrix_20251119-021143.csv"
ai_path       <- "missmap_results/data_availability_matrix_20251119-021737.csv"

baseline <- read_csv(baseline_path, show_col_types = FALSE)
ai       <- read_csv(ai_path,       show_col_types = FALSE)

# Helper to count synonyms in the 'synonyms' column
count_synonyms <- function(x) {
  x <- ifelse(is.na(x), "", x)
  ifelse(
    str_trim(x) == "",
    0L,
    lengths(strsplit(x, ","))
  )
}
syn_df <- baseline %>%
  select(Species, synonyms) %>%
  mutate(baseline_syn = count_synonyms(synonyms)) %>%
  select(Species, baseline_syn) %>%
  full_join(
    ai %>%
      select(Species, synonyms) %>%
      mutate(ai_syn = count_synonyms(synonyms)) %>%
      select(Species, ai_syn),
    by = "Species"
  ) %>%
  mutate(
    baseline_syn = replace_na(baseline_syn, 0L),
    ai_syn       = replace_na(ai_syn, 0L),
    ai_added     = pmax(ai_syn - baseline_syn, 0L)
  )
syn_long <- syn_df %>%
  # order species by total synonyms (AI run)
  mutate(Species = forcats::fct_reorder(Species, ai_syn, .desc = TRUE)) %>%
  select(Species, baseline_syn, ai_added) %>%
  pivot_longer(
    cols = c(baseline_syn, ai_added),
    names_to = "component",
    values_to = "count"
  ) %>%
  mutate(
    component = recode(
      component,
      baseline_syn = "NCBI synonyms",
      ai_added     = "AI-added synonyms"
    )
  )
library(tidyverse)

# syn_df has columns: Species, baseline_syn, ai_syn, ai_added
# syn_long is the long version you already made
# (Species, component = "NCBI synonyms"/"AI-added synonyms", count)

# make sure species order matches between syn_long and totals
totals <- syn_df %>%
  mutate(Species = forcats::fct_reorder(Species, ai_syn, .desc = TRUE))

syn_long_plot <- syn_long %>%
  mutate(Species = forcats::fct_reorder(Species, ai_syn, .desc = TRUE))

ggplot(syn_long_plot, aes(x = Species, y = count, fill = component)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "NCBI synonyms"     = "#1f78b4",  # blue
      "AI-added synonyms" = "#e31a1c"   # magenta/pink
    )
  ) +
  scale_y_continuous(
    breaks = 0:8,        # show 0,1,2,...,8
    limits = c(0, 8)     # buffer above max (6)
  ) +
  # 👉 numbers at the end of each bar (total synonyms per species)
  geom_text(
    data = totals,
    aes(
      x     = Species,
      y     = ai_syn + 0.2,   # a bit past the bar end
      label = ai_syn
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  labs(
    title = "Synonym sources per species in MissMap",
    x = "Species",
    y = "Number of synonyms",
    fill = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic"),
    plot.title  = element_text(hjust = 0.5)
  )
