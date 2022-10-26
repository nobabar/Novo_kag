library(tidyverse)
library(magrittr)
train_data <- read_csv("./train_tot_tsne.csv")
## shapes <- unique(train_data$group)
## shapes_df <- tibble("group" = shapes, "shape" = rep(1:5, len = 156))
## train_data <- inner_join(train_data, shapes_df, by = "group")

p <- train_data %>%
  mutate(gid = as.factor(gid)) %>%
  ggplot(aes(x = tsne_x, y = tsne_y)) +
  geom_point(aes(color = gid, shape = gid)) +
  scale_shape_manual(values = rep(15:18, len = 156))
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("tsne.png", width = 21, height = 9)


p2 <- train_data %>%
  filter(gid %in% c(147, 114)) %>%
  mutate(gid = as.factor(gid)) %>%
  ggplot(aes(x = tsne_x, y = tsne_y)) +
  geom_point(aes(color = gid, shape = gid)) +
  scale_shape_manual(values = rep(15:18, len = 156))
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pb_data <- train_data %>%
  filter(gid %in% c(0, 85))

pb_data_all <- train_data %>%
  group_by(tsne_x, tsne_y) %>%
  filter(n()>1) %>%
  ungroup() %>%
  group_by(gid) %>%
  summarise(n = n())

pb_data_all %>%
  mutate(gid = as.factor(gid)) %>%
  ggplot(aes(x = gid, y = n)) +
  geom_bar(stat = "identity")

ggsave("pb_data.png")
