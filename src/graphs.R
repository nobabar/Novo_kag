library(tidyverse)

train_data <- read_csv("./train.csv")
train_data_upd <- read_csv("./train_updates_20220929.csv")
ids_to_delete <- train_data_upd %>%
  filter(across(c(protein_sequence, pH, data_source, tm), ~is.na(.x))) %>%
  pull(seq_id)
ids_to_switch <- train_data_upd %>%
  filter(across(c(protein_sequence, pH, tm), ~!is.na(.x))) %>%
  pull(seq_id)

train_data[train_data$seq_id %in% ids_to_switch, c("tm", "pH")] <-
  train_data[train_data$seq_id %in% ids_to_switch, c("pH", "tm")]

train_data_updated <- train_data[!(train_data$seq_id %in% ids_to_delete), ]

train_data_updated %>%
  ggplot(aes(x = pH)) +
  geom_histogram()  +
  theme_bw()

ggsave("ph_plot.png")

train_data_updated %>%
  ggplot(aes(x = tm)) +
  geom_histogram()  +
  theme_bw()

ggsave("tm_plot.png")

train_data_updated$len_prot <- map_int(train_data_updated$protein_sequence, nchar)

train_data_updated %>%
  ggplot(aes(x = len_prot)) +
  geom_histogram(binwidth = 30)  +
  theme_bw()


ggsave("len_plot.png")

train_data_updated %>%
  filter(len_prot < 3000) %>%
  ggplot(aes(x = len_prot, y = tm)) +
    geom_point(alpha = .1) +
    theme_bw()


ggsave("cor_plot.png")

lm(tm ~ len_prot, data = train_data_updated) %>% summary
