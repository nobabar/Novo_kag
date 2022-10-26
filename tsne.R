library(tidyverse)
library(magrittr)
train_data <- read_csv("./train_tot_tsne.csv")
train_data %<>%
  mutate(group = as.factor(group))
shapes <- unique(train_data$group)
shapes_df <- tibble("group" = shapes, "shape" = rep(1:5, len = 156))
train_data <- inner_join(train_data, shapes_df, by = "group")

train_data %>%
  ggplot(aes(x = tsne_x, y = tsne_y)) +
  geom_point(aes(color = group, shape = group)) +
  scale_shape_manual(values = rep(15:18, len = 156))
  theme_bw() +
  theme(panel.grid.minor = element_blank())
