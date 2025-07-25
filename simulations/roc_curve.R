library(tidyverse)

setwd("path/spectre-paper/simulations")

get_pct <- function(n, m) {
  checks <- c(1, 5, 10, 25)
  checks_ <- c("01", "05", "10", "25")
  N <- 100*round(n*m, 2)
  diffs <- abs(checks - N)
  checks_[which(diffs == min(diffs))]
}
get_pct_V <- Vectorize(get_pct)

dat_phantom <- read.csv("results-phantom-epistasis/aggregated_phantom_epistasis_results_relate.csv") %>%
  filter(!is.na(test2_overall_min_b))
dat_true <- read.csv("results-true-interactions/aggregated_true_interactions_results_relate.csv") %>%
  filter(!is.na(freq_snp1)) %>%
  mutate(frequency_category = paste0(get_pct_V(AF_A, AF_B), "_pct")) %>%
  mutate(dist_category = round(distance/1000000, 1)) %>%
  mutate(dist_category = ifelse(dist_category > 3.0, 5.0, dist_category))
dat_phantom_true <- read.csv("results-phantom-epistasis/aggregated_phantom_epistasis_results_true-trees.csv") %>%
  filter(!is.na(test2_overall_min_b))

dat_phantom %>%
  group_by(frequency_category) %>%
  summarise(n = n())
dat_true %>%
  group_by(frequency_category, dist_category) %>%
  summarise(n = n())
dat_phantom_true %>%
  group_by(frequency_category) %>%
  summarise(n = n())

g1 <- dat_phantom_true %>%
  mutate(b1a_true = test1a_b_clade_min_p0.5,
         b1b_true = test1b_b_min) %>%
  select(run_id, b1a_true, b1b_true)
g2 <- dat_phantom %>%
  mutate(b1a_relate = test1a_b_clade_min_p0.5,
         b1b_relate = test1b_b_min) %>%
  select(run_id, b1a_relate, b1b_relate)
g <- full_join(g1, g2)

ggplot(data=g %>% filter(b1a_relate != 0)) +
  geom_point(aes(x=b1a_true, y=b1a_relate)) +
  geom_abline(slope=1) +
  theme_bw() +
  scale_x_continuous(name="True trees", limits=c(0.1, 0.6)) +
  scale_y_continuous(name="Relate trees", limits=c(0.1,0.6)) +
  ggtitle("Test 1a")
  
ggplot(data=g %>% filter(b1b_relate != 0)) +
  geom_point(aes(x=b1b_true, y=b1b_relate)) +
  geom_abline(slope=1) +
  theme_bw() +
  scale_x_continuous(name="True trees", limits=c(0.0, 0.2)) +
  scale_y_continuous(name="Relate trees", limits=c(0.0,0.2)) +
  ggtitle("Test 1b")

b_range <- seq(0,2,0.1)
roc_curve <- data.frame(
  true_test1a=numeric(0), phantom_test1a=numeric(0), 
  true_test1b=numeric(0), phantom_test1b=numeric(0), 
  true_test2=numeric(0), phantom_test2=numeric(0),
  frequency_category=character(0),
  dist_category=character(0),
  b=numeric(0)
  )
fc <- c("01_pct", "05_pct", "10_pct", "25_pct")
fc_labels = c("AF=0.01", "AF=0.05", "AF=0.10", "AF=0.25")
dist <- c(0.1, 0.5, 5)
dist_labels <- c("100kb (Relate trees)", "500kb (Relate trees)", "5Mb (Relate trees)")
for (i in 1:length(b_range)) {
  b <- b_range[i]
  for (j in 1:length(fc)) {
    f <- fc[j]
    for (k in 1:length(dist)) {
      d <- dist[k]
      true_t1a <- dat_true$test1a_b_clade_min_p0.5[dat_true$frequency_category == f & dat_true$dist_category == d]
      phantom_t1a <- dat_phantom$test1a_b_clade_min_p0.5[dat_phantom$frequency_category == f]
      true_t1b <- dat_true$test1b_b_min[dat_true$frequency_category == f & dat_true$dist_category == d]
      phantom_t1b <- dat_phantom$test1b_b_min[dat_phantom$frequency_category == f]
      true_t2 <- dat_true$test2_overall_min_b[dat_true$frequency_category == f & dat_true$dist_category == d]
      phantom_t2 <- dat_phantom$test2_overall_min_b[dat_phantom$frequency_category == f]
      roc_curve <- rbind(roc_curve, data.frame(
        true_test1a = sum(true_t1a < b)/length(true_t1a), 
        phantom_test1a = sum(phantom_t1a < b)/length(phantom_t1a), 
        true_test1b = sum(true_t1b < b)/length(true_t1b), 
        phantom_test1b = sum(phantom_t1b < b)/length(phantom_t1b), 
        true_test2 = sum(true_t2 < b)/length(true_t2), 
        phantom_test2 = sum(phantom_t2 < b)/length(phantom_t2), 
        frequency_category=fc_labels[j],
        dist_category=dist_labels[k],
        b=b))
    }
  }
}

ggplot(data=roc_curve) +
  geom_line(aes(x=true_test1a, y=phantom_test1a, color=frequency_category), linewidth=1, alpha=0.8) +
  scale_x_continuous(name="False positive rate", limits=c(0, 1)) +
  scale_y_continuous(name="True positive rate", limits=c(0,1)) +
  geom_abline(slope=1) +
  facet_wrap(~dist_category) +
  theme_bw() +
  theme(legend.title=element_blank())

ggplot(data=roc_curve) +
  geom_line(aes(x=true_test1b, y=phantom_test1b, color=frequency_category), linewidth=1, alpha=0.8) +
  scale_x_continuous(name="False positive rate", limits=c(0, 1)) +
  scale_y_continuous(name="True positive rate", limits=c(0,1)) +
  geom_abline(slope=1) +
  facet_wrap(~dist_category) +
  theme_bw() +
  theme(legend.title=element_blank())


