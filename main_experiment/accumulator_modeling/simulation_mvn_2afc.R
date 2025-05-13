library(doParallel)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tictoc)
library(cowplot)
library(sjPlot)
library(ggsci)
library(ggh4x)
library(matrixStats)

theme_set(theme_publication()) 

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

mu_d <-    0.0000798
sigma <-   0.0176
tau_nd <-  0.49
sigma_d <- 0.000005
beta <-    0.5
theta <-   0.156
w <-     0
th <-    1

trial <-  200000
sample <- 3501
condition <- c(92, 96, 100)
simdat <- c()

tic()
for (j in condition) {
  
  dr_scaler <- cbind(rnorm(trial, mu_d, sigma_d),
                     rnorm(trial, mu_d, sigma_d)) 
  
  a10 <- runif(trial, -beta, beta)
  a20 <- runif(trial, -beta, beta)
  
  d <- foreach(i = 1:trial, .combine = "rbind", .packages = c("dplyr", "matrixStats")) %dopar% {
    
    stimulus <- c(j, 90)
    # stimulus <- stimulus / (1 + w * sum(stimulus))
    stimulus <- stimulus * dr_scaler[i, ]
    
    accum <- matrix(0, nrow = sample + 1, ncol = 2)
    accum[, 1] <- c(a10[i], rep(stimulus[1], sample))
    accum[, 2] <- c(a20[i], rep(stimulus[2], sample))
    accum[-1, 1] <- accum[-1, 1] + rnorm(sample, 0, sigma)
    accum[-1, 2] <- accum[-1, 2] + rnorm(sample, 0, sigma)
    accum <- matrixStats::colCumsums(accum)
    
    mvn <- matrixStats::rowOrderStats(-accum, which = 2) - matrixStats::rowOrderStats(-accum, which = 1)
    simrt  <- which(mvn > th)[1]
    choice <- which.max(accum[simrt, ])
    second <- which.max(-accum[simrt, ])
    e_chosen_post <- ifelse(simrt <= 3000, accum[simrt + 350, choice], NA)
    e_second_post <- ifelse(simrt <= 3000, accum[simrt + 350, second], NA)
    as.numeric(c(simrt, choice, mvn[simrt], mvn[1], mvn[350], e_chosen_post, e_second_post))                                                  
    
  }
  
  d <- as.data.frame(d)
  colnames(d) <- c("rt", "choice", "mvn_choice", "mvn_0", "mvn_350", "e_chosen_post", "e_second_post")
  d$rt   <- d$rt + tau_nd * 1000
  d$condition <- j
  simdat <- rbind(simdat, filter(d, rt <= 3000))
  
}
toc()

# model visualization
df <- read.csv("exp2_data_2afc_processed.csv", header = T)
df$rt <- as.numeric(df$rt)
df <- na.omit(df)
df$chosenItem <- factor(df$chosenItem, levels = c("target", "distractor"))
df <- mutate(df, confidence = ifelse(conf > 2, 1, 0))

simdat <- na.omit(simdat)
simdat <- mutate(simdat, early_conf = mvn_350 - mvn_0)
simdat <- mutate(simdat, post_conf  = e_chosen_post - e_second_post - mvn_choice)
simdat <- mutate(simdat, cv = 0.4 * early_conf + 0.6 * post_conf)
simdat <- mutate(simdat, confidence = ifelse(cv > theta, "high confidence", "low confidence"))
simdat$choice <- ifelse(simdat$choice == 1, "target", "distractor")
simdat$choice <- factor(simdat$choice, levels = c("target", "distractor"))
simdat <- rename(simdat, chosenItem = choice, val1 = condition)


# RT histogram
df %>%
  tidyr::complete(chosenItem, val1) %>%
  dplyr::select(chosenItem, val1, rt) %>%
  mutate(rt = rt * 1000) -> emp1
emp1$chosenItem <- factor(emp1$chosenItem, levels = c("target", "distractor", "dud"))

simdat %>%
  complete(chosenItem, val1) %>%
  ggplot() +
  geom_histogram(emp1, mapping = (aes(x = rt, fill = chosenItem, color = chosenItem))) +
  geom_freqpoly(aes(x = rt, y = after_stat(count / nrow(simdat) * nrow(df)))) +
  scale_color_discrete(labels = c("target", "distractor")) +
  guides(fill = F, color = F) + xlab("RT (s)") + ylab("Count") +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000), labels = c("0", "1", "2", "3")) +
  facet_wrap(. ~ factor(val1) + factor(chosenItem), nrow = 3) +
  theme(axis.text = element_text(size = 6.5), axis.title = element_text(size = 7)) +
  scale_fill_npg() + scale_color_npg() -> p1
p1

# choice proportion
simdat %>%
  dplyr::select(chosenItem, val1) %>%
  group_by(chosenItem, val1) %>%
  mutate(n1 = n()) %>%
  ungroup(chosenItem, val1) %>%
  group_by(val1) %>%
  mutate(n2 = n()) %>%
  dplyr::select(n1, n2, chosenItem, val1) %>%
  distinct() %>%
  mutate(p = n1/n2) %>%
  ungroup(chosenItem, val1) %>%
  complete(chosenItem, val1) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_all(~replace(., is.na(.), 0)) -> m2
m2$chosenItem <- factor(m2$chosenItem, levels = c("target", "distractor"))

df %>%
  group_by(chosenItem, val1, subj) %>%
  mutate(n1 = n()) %>%
  ungroup(chosenItem, val1, subj) %>%
  group_by(val1, subj) %>%
  mutate(n2 = n()) %>%
  dplyr::select(n1, n2, val1, chosenItem, subj) %>%
  distinct() %>%
  mutate(p = n1/n2) %>%
  ungroup(chosenItem, val1, subj) %>%
  complete(chosenItem, val1, subj) %>%
  mutate_all(~ replace_na(., 0)) %>%
  ggplot() + 
  geom_line(m2, size = 0.4, mapping = aes(x = val1, y = p, color = chosenItem)) +
  stat_summary(fun = "mean", geom = "point", aes(x = val1, y = p, color = chosenItem)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val1, y = p, color = chosenItem)) +
  labs(color = "choice") + xlab("Target value") +
  theme(legend.position = c(.5, .6), legend.text = element_text(size = 10), 
        legend.background = element_rect(fill = NA, colour = NA), legend.key = element_rect(fill = NA)) + 
  guides(color = guide_legend(title = NULL)) + scale_color_npg() +
  coord_cartesian(ylim = c(0, 0.82)) +
  scale_x_continuous(breaks = c(92, 96, 100)) + ylab("Choice proportion") + 
  guides(color = F) + facet_wrap(chosenItem ~ .) -> p2
p2

# mean RT
simdat %>%
  dplyr::select(rt, chosenItem, val1) %>%
  group_by(chosenItem, val1) %>%
  summarise(mean_rt = mean(rt) / 1000) -> m3
m3$chosenItem <- factor(m3$chosenItem, levels = c("target", "distractor"))

df %>% group_by(val1, chosenItem, subj) %>%
  summarise(mean_rt = mean(rt)) %>%
  ggplot() + 
  geom_line(m3, size = 0.4, mapping = aes(x = val1, y = mean_rt, color = chosenItem)) +
  stat_summary(fun = "mean", geom = "point", aes(x = val1, y = mean_rt, color = chosenItem)) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val1, y = mean_rt, color = chosenItem)) +
  guides(color = F) + xlab("Target value") + ylab("Mean RT (s)") + scale_color_npg() +
  scale_y_continuous(breaks = seq(0.7, 2.1, 0.3)) + coord_cartesian(ylim = c(0.7, 2.1)) +
  scale_x_continuous(breaks = c(92, 96, 100)) + facet_wrap(chosenItem ~ .) -> p3
p3

# high confidence rate
simdat %>%
  dplyr::select(confidence, chosenItem, val1) %>%
  mutate(conf = ifelse(confidence == "high confidence", 1, 0)) %>%
  group_by(chosenItem, val1) %>%
  summarise(high_conf_rate = mean(conf)) -> m4
m4$chosenItem <- factor(m4$chosenItem, levels = c("target", "distractor"))

df %>% group_by(val1, chosenItem, subj) %>%
  summarise(mean_confidence = mean(conf)) %>%
  ggplot() + 
  stat_summary(fun = "mean", geom = "point", aes(x = val1, y = mean_confidence, color = chosenItem)) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val1, y = mean_confidence, color = chosenItem)) +
  xlab("Target value") + ylab("Mean confidence") + guides(color = F) + scale_color_npg() +
  scale_x_continuous(breaks = c(92, 96, 100)) + facet_wrap(chosenItem ~ .) -> p4
p4

df %>% group_by(val1, chosenItem, subj) %>%
  summarise(high_conf = mean(confidence)) %>%
  ggplot() + 
  geom_line(m4, size = 0.4, mapping = aes(x = val1, y = high_conf_rate, color = chosenItem)) +
  stat_summary(fun = "mean", geom = "point", aes(x = val1, y = high_conf, color = chosenItem)) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val1, y = high_conf, color = chosenItem)) +
  xlab("Target value") + ylab("High confidence rate") + guides(color = F) + scale_color_npg() +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2)) + coord_cartesian(ylim = c(0, 0.8)) +
  scale_x_continuous(breaks = c(92, 96, 100)) + facet_wrap(chosenItem ~ .) -> p4
p4

# rt for low confidence responses
simdat %>% group_by(chosenItem, val1, confidence) %>%
  summarise(mean_rt = mean(rt) / 1000) %>% 
  filter(confidence == "low confidence") -> m5
m5$chosenItem <- factor(m5$chosenItem, levels = c("target", "distractor"))

df <- mutate(df, confidence = ifelse(confidence == 1, "high confidence", "low confidence"))
df %>% group_by(val1, chosenItem, subj, confidence) %>%
  summarise(mean_rt = mean(rt)) %>%
  filter(confidence == "low confidence") %>%
  ggplot() + 
  geom_line(m5, size = 0.4, mapping = aes(x = val1, y = mean_rt, color = chosenItem)) +
  stat_summary(fun = "mean", geom = "point", aes(x = val1, y = mean_rt, color = chosenItem)) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val1, y = mean_rt, color = chosenItem)) +
  xlab("Target value") + ylab("Mean RT (s)") + guides(color = F) + scale_color_npg() +
  scale_y_continuous(breaks = seq(0.7, 2.1, 0.3)) + coord_cartesian(ylim = c(0.7, 2.1)) +
  scale_x_continuous(breaks = c(92, 96, 100)) + facet_nested(. ~ confidence + chosenItem) -> p5
p5

# rt for high confidence responses
simdat %>% group_by(chosenItem, val1, confidence) %>%
  summarise(mean_rt = mean(rt) / 1000) %>% 
  filter(confidence == "high confidence") -> m6
m6$chosenItem <- factor(m6$chosenItem, levels = c("target", "distractor"))

df %>% group_by(val1, chosenItem, subj, confidence) %>%
  summarise(mean_rt = mean(rt)) %>%
  filter(confidence == "high confidence") %>%
  ggplot() + 
  geom_line(m6, size = 0.4, mapping = aes(x = val1, y = mean_rt, color = chosenItem)) +
  stat_summary(fun = "mean", geom = "point", aes(x = val1, y = mean_rt, color = chosenItem)) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val1, y = mean_rt, color = chosenItem)) +
  xlab("Target value") + ylab("Mean RT (s)") + guides(color = F) + scale_color_npg() +
  scale_y_continuous(breaks = seq(0.7, 2.1, 0.3)) + coord_cartesian(ylim = c(0.7, 2.1)) +
  scale_x_continuous(breaks = c(92, 96, 100)) + facet_nested(. ~ confidence + chosenItem) -> p6
p6

# early confidence
simdat %>%
  group_by(chosenItem, val1) %>%
  summarise(mean_ec = mean(early_conf)) %>%
  ggplot() + geom_line(aes(x = val1, y = mean_ec, color = factor(chosenItem))) +
  scale_color_discrete(labels = c("target", "distractor")) +
  scale_x_continuous(breaks = c(92, 96, 100)) + 
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) + coord_cartesian(ylim = c(-0.45, 0.4)) +
  theme(legend.position = c(.37, .46), legend.text = element_text(size = 10), 
        legend.background = element_rect(fill = NA, colour = NA), legend.key = element_rect(fill = NA)) + 
  xlab("Target value") + ylab("Early-stage MVN") + 
  guides(color = guide_legend(title = NULL, nrow = 3, keyheight = unit(0.5, "lines"))) +
  scale_color_npg() + guides(color = F) + facet_wrap(chosenItem ~ .) -> p7

# post confidence
simdat %>%
  group_by(chosenItem, val1) %>%
  summarise(mean_pc = mean(post_conf)) %>%
  ggplot() + geom_line(aes(x = val1, y = mean_pc, color = factor(chosenItem))) +
  scale_color_discrete(labels = c("target", "distractor")) +
  scale_x_continuous(breaks = c(92, 96, 100)) + 
  xlab("Target value") + ylab("Post-decisional MVN") + 
  scale_color_npg() + guides(color = F) + facet_wrap(chosenItem ~ .) -> p8

# save figures
figure_6 <- cowplot::plot_grid(
    p2 + theme(legend.text = element_text(size = 6), 
               legend.background = element_rect(fill = NA, colour = NA),
               legend.key = element_rect(fill = NA),
               axis.text.x = element_text(angle = 20),
               plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
    p3 + theme(axis.text.x = element_text(angle = 20),
               plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
    p4 + theme(axis.text.x = element_text(angle = 20),
               plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
    p7 + theme(axis.text.x = element_text(angle = 20),
               plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
    p8 + theme(axis.text.x = element_text(angle = 20),
               plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
    ggplot() + theme_void(),
    p5 + theme(axis.text.x = element_text(angle = 20),
               plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
    p6 + theme(axis.text.x = element_text(angle = 20),
               plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
    labels = c("a", "b", "c", "d", "e", "", "f", "g"),
    label_size = 10, 
    nrow = 3
)

figure_6
sjPlot::save_plot("figure_6.jpg", figure_6, dpi = 600, height = 13, width = 14)
sjPlot::save_plot("figure_S4.jpg", p1, dpi = 600)
write.csv(simdat, "simdat_mvn_2afc.csv")