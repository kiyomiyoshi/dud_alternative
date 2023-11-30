library(doParallel)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tictoc)
library(cowplot)
library(sjPlot)
library(ggsci)
library(matrixStats)

theme_set(theme_publication()) 

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

mu_d    <- 0.000167
sigma <-   0.056
tau_nd  <- 0.32
sigma_d <- 0.000004
beta <-    0.4
theta <-   0.305
w     <-   0
th    <-  -0.1 # should be negative

trial <- 50000
sample <- 3501
condition <- c(0, 27, 45, 63, 76, 86)
simdat <- c()

tic()
for (j in condition) {
  
  dr_scaler <- cbind(rnorm(trial, mu_d, sigma_d),
                     rnorm(trial, mu_d, sigma_d), 
                     rnorm(trial, mu_d, sigma_d))
  
  a10 <- runif(trial, -beta, beta)
  a20 <- runif(trial, -beta, beta) 
  a30 <- runif(trial, -beta, beta) 
  
  d <- foreach(i = 1:trial, .combine = "rbind", .packages = c("matrixStats")) %dopar% {
    
    stimulus <- c(100, 90, j)
    # stimulus <- stimulus / (1 + w * sum(stimulus))
    stimulus <- stimulus * dr_scaler[i, ]
    nAlt <- sum(stimulus != 0)
    
    accum <- matrix(0, nrow = sample + 1, ncol = 3)
    accum[, 1] <- c(a10[i], rep(stimulus[1], sample))
    accum[, 2] <- c(a20[i], rep(stimulus[2], sample))
    accum[, 3] <- c(a30[i], rep(stimulus[3], sample))
    accum[-1, 1] <- accum[-1, 1] + rnorm(sample, 0, sigma)
    accum[-1, 2] <- accum[-1, 2] + rnorm(sample, 0, sigma)
    accum[-1, 3] <- accum[-1, 3] + rnorm(sample, 0, sigma)
    accum <- accum[, 1:nAlt]
    accum <- matrixStats::colCumsums(accum)
    
    ppa <- matrixStats::rowMaxs(accum) - log(rowSums(exp(accum)))
    simrt  <- which(ppa > th)[1]
    choice <- which.max(accum[simrt, ])
    third  <- which.max(-accum[simrt, ])
    second <- ifelse(nAlt == 3, setdiff(c(1, 2, 3), c(choice, third)), third)
    # pre <- ifelse(simrt <= 300, simrt - 1, 300)
    ppa_post <- ifelse(simrt <= 3000, accum[simrt + 350, choice] - log(sum(exp(accum[simrt + 350, ]))), NA) # choice-conditional post ppa
    as.numeric(c(simrt, choice, ppa[simrt], ppa[1], ppa[350], ppa_post))   
    
  }
  
  d <- as.data.frame(d)
  colnames(d) <- c("rt", "choice", "ppa_choice", "ppa_0", "ppa_350", "ppa_post")
  d$rt   <- d$rt + tau_nd * 1000
  d$condition <- j
  simdat <- rbind(simdat, filter(d, rt <= 3000))
  
}
toc()


# model visualization
df <- read.csv("exp1_data_behavior.csv", header = TRUE)
df <- na.omit(df)

simdat <- na.omit(simdat)
simdat <- mutate(simdat, early_conf = qnorm(exp(ppa_350)) - qnorm(exp(ppa_0)))
simdat <- mutate(simdat, post_conf  = qnorm(exp(ppa_post)) - qnorm(exp(ppa_choice)))
simdat <- mutate(simdat, cv = 0.6 * early_conf + 0.4 * post_conf)
simdat <- mutate(simdat, confidence = ifelse(cv > theta, "high", "low"))
simdat$choice <- ifelse(simdat$choice == 1, "target", 
                        ifelse(simdat$choice == 2, "distractor", "dud"))
simdat$choice <- factor(simdat$choice, levels = c("target", "distractor", "dud"))

# RT histogram
df %>%
  complete(chosenItem, val3) %>%
  dplyr::select(chosenItem, val3, rt) %>%
  rename(choice = chosenItem, condition = val3) %>%
  mutate(rt = rt * 1000) -> emp1
emp1$choice <- factor(emp1$choice, levels = c("target", "distractor", "dud"))

simdat %>%
  complete(choice, condition) %>%
  ggplot() +
  geom_histogram(emp1, mapping = (aes(x = rt, fill = choice, color = choice))) +
  geom_freqpoly(aes(x = rt, y = after_stat(count / nrow(simdat) * nrow(df)))) +
  scale_color_discrete(labels = c("target", "distractor", "dud")) +
  guides(fill = F, color = F) + xlab("RT (s)") + ylab("Count") +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000), labels = c("0", "1", "2", "3")) +
  facet_wrap(. ~ factor(condition) + factor(choice), nrow = 3) +
  theme(axis.text = element_text(size = 6.5), axis.title = element_text(size = 7)) +
  scale_fill_npg() + scale_color_npg() -> p1
p1

# choice proportion
df %>%
  dplyr::select(chosenItem, val3) %>%
  group_by(chosenItem, val3) %>%
  rename(choice = chosenItem, condition = val3) %>%
  mutate(n1 = n()) %>%
  ungroup(choice, condition) %>%
  group_by(condition) %>%
  mutate(n2 = n()) %>%
  dplyr::select(n1, n2, choice, condition) %>%
  distinct() %>%
  mutate(p = n1/n2) %>%
  ungroup(choice, condition) %>%
  complete(choice, condition) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_all(~replace(., is.na(.), 0)) -> emp2

simdat %>%
  group_by(choice, condition) %>%
  mutate(n1 = n()) %>%
  ungroup(choice, condition) %>%
  group_by(condition) %>%
  mutate(n2 = n()) %>%
  dplyr::select(n1, n2, choice, condition) %>%
  distinct() %>%
  mutate(p = n1/n2) %>%
  ungroup(choice, condition) %>%
  complete(choice, condition) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  ggplot() + geom_line(aes(x = condition, y = p, color = factor(choice))) + ylab("Choice proportion") +
  geom_point(emp2, size = 0.8, mapping = aes(x = condition, y = p, color = choice)) + xlab("Dud value") +
  scale_color_discrete(labels = c("target", "distractor", "dud")) + labs(color = "choice") +
  theme(legend.text = element_text(size = 7),
        legend.position = c(.5, .6), 
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA)) +
  labs(color = NULL) +
  scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 0.81) + scale_color_npg() -> p2

# target vs. distractor choice
df %>%
  dplyr::select(chosenItem, val3) %>%
  rename(choice = chosenItem, condition = val3) %>%
  group_by(condition) %>%
  summarise(n1 = sum(choice == "target"), n2 = sum(choice == "distractor"), n3 = sum(choice == "dud")) %>%
  mutate(r = n1/n2) -> emp3

simdat %>%
  group_by(condition) %>%
  summarise(n1 = sum(choice == "target"), n2 = sum(choice == "distractor"), n3 = sum(choice == "dud")) %>%
  mutate(r = n1/n2) %>% 
  ggplot() + geom_line(aes(x = condition, y = r)) + geom_point(emp3, size = 0.8, mapping = aes(x = condition, y = r)) +
  scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(3, 5) + xlab("Dud value") +
  theme(legend.position = c(.25, .8), legend.background = element_rect(fill = NA, colour = NA)) + 
  guides(linetype = guide_legend(title = NULL)) +
  ylab("Pc(target) / Pc(distractor)") + scale_color_npg() -> p3

# mean RT
df %>%
  dplyr::select(rt, chosenItem, val3) %>%
  rename(choice = chosenItem, condition = val3) %>%
  group_by(choice, condition) %>%
  summarise(mean_rt = mean(rt) * 1000) -> emp4
emp4$choice <- factor(emp4$choice, levels = c("target", "distractor", "dud"))

simdat %>%
  group_by(choice, condition) %>%
  summarise(mean_rt = mean(rt)) %>%
  
  ggplot() + geom_line(aes(x = condition, y = mean_rt / 1000, color = factor(choice))) +
  geom_point(emp4, size = 0.8, mapping = aes(x = condition, y = mean_rt / 1000, color = factor(choice))) +
  scale_color_discrete(labels = c("target", "distractor", "dud")) + labs(color = "choice") +
  scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + xlab("Dud value") +
  ylab("Mean RT (s)") + ylim(0.8, 1.6) + scale_color_npg() -> p4

# high confidence rate
df %>%
  dplyr::select(conf, chosenItem, val3) %>%
  rename(choice = chosenItem, condition = val3) %>%
  group_by(choice, condition) %>%
  mutate(conf = ifelse(conf > 2, 1, 0)) %>%
  summarise(high_conf_rate = mean(conf)) -> emp5
emp5$choice <- factor(emp5$choice, levels = c("target", "distractor", "dud"))

simdat %>%
  group_by(choice, condition) %>%
  mutate(confidence = ifelse(confidence == "high", 1, 0)) %>%
  summarise(high_conf_rate = mean(confidence)) %>%
  ggplot() + geom_line(aes(x = condition, y = high_conf_rate, color = factor(choice))) +
  geom_point(emp5, mapping = aes(x = condition, y = high_conf_rate, color = factor(choice))) +
  scale_color_discrete(labels = c("target", "distractor", "dud")) + labs(color = "choice") +
  scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + xlab("Dud value") +
  scale_y_continuous(breaks = seq(0, 0.9, 0.1), limits = c(0, 0.9)) +
  ylab("High confidence rate") + scale_color_npg() -> p5

# save figures
figure_ppa <- cowplot::plot_grid(p2 + guides(color = F, fill = F) + 
                                   scale_color_npg(guide = guide_legend(nrow = 2)), 
                                 p3 + guides(color = F) +
                                   theme(legend.background = element_rect(fill = NA, colour = NA),
                                         legend.key = element_rect(fill = NA)),
                                 p4 + guides(color = F),
                                 p5 + guides(color = F),
                                 labels = c("a", "b", "c", "d"), label_size = 10, nrow = 2)
figure_ppa 

# sjPlot::save_plot("ppa.jpg", figure_ppa, dpi = 600, height = 9)
# sjPlot::save_plot("ppa_rt.jpg", p1, dpi = 600)
# write.csv(simdat, "simdat_ppa.csv")