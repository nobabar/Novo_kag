library(tidyverse)
library(magrittr)
library(ggrepel)

train_data <- read_csv("./train_tot_tsne.csv")
wildtype_data <- read_csv("./wildtype_tsne.csv")
train_data_2 <- train_data %>%
  select(gid, tsne_x, tsne_y) %>%
  mutate(origin = "mutated")
wildtype_data_2 <- wildtype_data %>%
  select(gid, tsne_x, tsne_y) %>%
  mutate(origin = "wildtype")

total_data <- rbind(train_data_2, wildtype_data_2)

## shapes <- unique(train_data$group)
## shapes_df <- tibble("group" = shapes, "shape" = rep(1:5, len = 156))
## train_data <- inner_join(train_data, shapes_df, by = "group")
p <- total_data %>%
  ggplot(aes(x = tsne_x, y = tsne_y)) +
  geom_point(aes(color = gid, shape = gid), show.legend = FALSE) +
  scale_shape_manual(values = rep(15:19, len = length(unique(total_data$gid)))) +
  geom_label_repel(
    data = total_data %>% filter(origin == "wildtype"),
    aes(label = gid, color = NULL),
    show.legend = FALSE,
    size = 2
    ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("tsne.png", width = 16, height = 9, dpi = 300)
