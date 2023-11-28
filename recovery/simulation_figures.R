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

Mu_d  <-   list(rep(0.0073, 24),  rnorm(24, 0.0073, 0.0073 / 10),   rnorm(24, 0.0073, 0.0073 / 7))
Sigma <-   list(rep(0.023, 24),   rnorm(24, 0.023, 0.023 / 10),     rnorm(24, 0.023, 0.023 / 7))
Tau_nd <-  list(rep(0.33, 24),    rnorm(24, 0.33, 0.33 / 10),       rnorm(24, 0.33, 0.33 / 7))
Sigma_d <- list(rep(0.00015, 24), rnorm(24, 0.00015, 0.00015 / 10), rnorm(24, 0.00015, 0.00015 / 7))
Beta <-    list(rep(0.25, 24),    rnorm(24, 0.25, 0.25 / 10),       rnorm(24, 0.25, 0.25 / 7))
Theta <-   list(rep(0.17, 24),    rnorm(24, 0.17, 0.17 / 10),       rnorm(24, 0.17, 0.17 / 7))
Weight <-  list(rep(0.75, 24),    rnorm(24, 0.75, 0.75 / 10),       rnorm(24, 0.75, 0.75 / 7))

th <-    1
trial <-  288
sample <- 3501
condition <- c(0, 27, 45, 63, 76, 86)
recdat <- c()

tic()
for (l in 1:3) {
    
    for (k in 1:24) {
        
        mu_d <-    Mu_d[[l]][k]
        sigma <-   Sigma[[l]][k]
        tau_nd <-  Tau_nd[[l]][k]
        sigma_d <- Sigma_d[[l]][k]
        beta <-    Beta[[l]][k]
        theta <-   Theta[[l]][k]
        weight <-  Weight[[l]][k]
        dat <- c()
        
        for (j in condition) {
            
            dr_scaler <- cbind(rnorm(trial, mu_d, sigma_d),
                               rnorm(trial, mu_d, sigma_d), 
                               rnorm(trial, mu_d, sigma_d)) 
            
            a10 <- runif(trial, -beta, beta)
            a20 <- runif(trial, -beta, beta)
            a30 <- runif(trial, -beta, beta)
            
            d <- foreach(i = 1:trial, .combine = "rbind", .packages = c("dplyr", "matrixStats")) %dopar% {
                
                stimulus <- c(100, 90, j) / 100
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
                
                mvn <- matrixStats::rowOrderStats(-accum, which = 2) - matrixStats::rowOrderStats(-accum, which = 1)
                simrt  <- which(mvn > th)[1]
                choice <- which.max(accum[simrt, ])
                third  <- which.max(-accum[simrt, ])
                second <- ifelse(nAlt == 3, setdiff(c(1, 2, 3), c(choice, third)), third)
                e_chosen_post <- ifelse(simrt <= 3000, accum[simrt + 350, choice], NA)
                e_second_post <- ifelse(simrt <= 3000, accum[simrt + 350, second], NA)
                cv <- weight * (mvn[350] - mvn[1]) + (1 - weight) * ((e_chosen_post - e_second_post) - mvn[simrt])
                c(simrt, choice, cv)                                              
                
            }
            
            d <- as.data.frame(d)
            colnames(d) <- c("rt", "choice", "cv")
            d$rt   <- d$rt + tau_nd * 1000
            d <- mutate(d, confidence = ifelse(cv > theta, 1, 0))
            d$condition <- j
            d$subject <- k
            d$iteration <- l
            dat <- rbind(dat, d)
        }
        recdat <- rbind(recdat, dat)
    }
}
toc()

3*24*6*288
table(recdat$iteration)

#'# visualization
df <- read.csv("exp1_data_behavior.csv", header = TRUE)
df <- na.omit(df)

recdat <- na.omit(recdat)
recdat$choice <- ifelse(recdat$choice == 1, "target", 
                        ifelse(recdat$choice == 2, "distractor", "dud"))
recdat$choice <- factor(recdat$choice, levels = c("target", "distractor", "dud"))

# RT histogram
df %>%
    complete(chosenItem, val3) %>%
    dplyr::select(chosenItem, val3, rt) %>%
    rename(choice = chosenItem, condition = val3) %>%
    mutate(rt = rt * 1000) -> emp1
emp1$choice <- factor(emp1$choice, levels = c("target", "distractor", "dud"))

recdat %>%
    complete(choice, condition) %>%
    ggplot() +
    geom_histogram(emp1, mapping = (aes(x = rt, fill = choice, color = choice))) +
    geom_freqpoly(aes(linetype = factor(iteration), x = rt, y = after_stat(count / nrow(recdat) * 3 * nrow(df)))) +
    scale_color_discrete(labels = c("target", "distractor", "dud")) +
    guides(fill = F, color = F) + xlab("RT (s)") + ylab("Count") +
    scale_x_continuous(breaks = c(0, 1000, 2000, 3000), labels = c("0", "1", "2", "3")) +
    facet_wrap(. ~ factor(condition) + factor(choice), nrow = 3) +
    theme(axis.text = element_text(size = 6.5), axis.title = element_text(size = 7)) +
    scale_fill_npg() + scale_color_npg() -> p1
p1

recdat %>%
    complete(choice, condition) %>%
    ggplot() +
    geom_freqpoly(aes(color = factor(iteration), linetype = factor(iteration),
                      x = rt, y = after_stat(count / nrow(recdat) * 3 * nrow(df)))) +
    scale_color_discrete(labels = c("0", "mean / 10", "mean / 7")) +
    xlab("RT (s)") + ylab("Count") + guides(linetype = FALSE) + labs(color = "SD") +
    scale_x_continuous(breaks = c(0, 1000, 2000, 3000), labels = c("0", "1", "2", "3")) +
    facet_wrap(. ~ factor(condition) + factor(choice), nrow = 3) +
    theme(axis.text = element_text(size = 6.5), axis.title = element_text(size = 7)) +
    theme(legend.position = "bottom", 
          legend.text = element_text(size = 7),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.key = element_rect(fill = NA)) -> p2
p2


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

recdat %>%
    group_by(choice, condition, iteration) %>%
    mutate(n1 = n()) %>%
    ungroup(choice, condition, iteration) %>%
    group_by(condition, iteration) %>%
    mutate(n2 = n()) %>%
    dplyr::select(n1, n2, choice, condition, iteration) %>%
    distinct() %>%
    mutate(p = n1/n2) %>%
    ungroup(choice, condition, iteration) %>%
    complete(choice, condition, iteration) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    ggplot() + geom_line(aes(x = condition, y = p, color = factor(choice), linetype = factor(iteration))) + ylab("Choice proportion") +
    scale_color_discrete(labels = c("target", "distractor", "dud")) + labs(color = "Choice", linetype = "SD") +
    theme(legend.text = element_text(size = 5),
          legend.position = c(.5, .63), 
          legend.direction = "horizontal",
          legend.background = element_rect(fill = NA, colour = NA),
          legend.key = element_rect(fill = NA)) +
    scale_linetype_discrete(labels = c("0", "mean / 10", "mean / 7")) + xlab("Dud value") +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(0, 0.8) + scale_color_npg() -> p3
p3


# target vs. distractor choice
df %>%
    dplyr::select(chosenItem, val3, subj) %>%
    rename(choice = chosenItem, condition = val3) %>%
    group_by(condition, subj) %>%
    summarise(n1 = sum(choice == "target"), n2 = sum(choice == "distractor"), n3 = sum(choice == "dud")) %>%
    mutate(r = n1/n2) %>%
    summarize(mean_r = mean(r)) -> emp3

recdat %>%
    group_by(condition, iteration) %>%
    summarise(n1 = sum(choice == "target"), n2 = sum(choice == "distractor"), n3 = sum(choice == "dud")) %>%
    mutate(r = n1/n2) %>% 
    ggplot() + geom_line(aes(x = condition, y = r, linetype = factor(iteration))) + 
    geom_point(emp3, size = 0.8, mapping = aes(x = condition, y = mean_r)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(3, 5) + xlab("Dud value") +
    theme(legend.position = c(.25, .8), legend.background = element_rect(fill = NA, colour = NA)) + 
    guides(linetype = guide_legend(title = "SD")) + scale_color_discrete(labels = c("0", "μ / 10", "μ / 7")) +
    ylab("Pc(target) / Pc(distractor)") -> p4
p4

recdat %>%
    group_by(condition, iteration) %>%
    summarise(n1 = sum(choice == "target"), n2 = sum(choice == "distractor"), n3 = sum(choice == "dud")) %>%
    mutate(r = n1/n2) %>% 
    ggplot() + geom_line(aes(x = condition, y = r, linetype = factor(iteration))) + 
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylim(3, 5) + xlab("Dud value") +
    theme(legend.position = c(.25, .8), legend.background = element_rect(fill = NA, colour = NA)) + 
    guides(linetype = guide_legend(title = "SD")) + scale_linetype_discrete(labels = c("0", "μ / 10", "μ / 7")) +
    ylab("Pc(target) / Pc(distractor)") +
    theme(legend.text = element_text(size = 7),
          legend.position = c(.5, .7), 
          legend.background = element_rect(fill = NA, colour = NA),
          legend.key = element_rect(fill = NA)) -> p5
p5


# mean RT
df %>%
    dplyr::select(rt, chosenItem, val3) %>%
    rename(choice = chosenItem, condition = val3) %>%
    group_by(choice, condition) %>%
    summarise(mean_rt = mean(rt) * 1000) -> emp4
emp4$choice <- factor(emp4$choice, levels = c("target", "distractor", "dud"))

recdat %>%
    group_by(choice, condition, iteration) %>%
    summarise(mean_rt = mean(rt)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_rt / 1000, color = factor(choice), linetype = factor(iteration))) +
  # geom_point(emp4, size = 0.8, mapping = aes(x = condition, y = mean_rt / 1000, color = factor(choice))) +
    scale_color_discrete(labels = c("target", "distractor", "dud")) + labs(color = "choice") +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + xlab("Dud value") +
    guides(linetype = FALSE, color = FALSE) +
    ylab("Mean RT (s)") + ylim(0.6, 1.4) + scale_color_npg() -> p6
p6

# confidence
recdat %>%
    group_by(choice, condition, iteration) %>%
    summarise(high_conf = mean(confidence)) %>%
    ggplot() + geom_line(aes(x = condition, y = high_conf, color = factor(choice), linetype = factor(iteration))) +
    scale_color_discrete(labels = c("target", "distractor", "dud")) + labs(color = "choice") +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + xlab("Dud value") +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75), limits = c(0, 0.75)) +
    guides(linetype = FALSE, color = FALSE) +
    ylab("High confidence rate") + scale_color_npg() -> p7
p7


# save figures
legend <- get_legend(p3 + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "bottom"))
g <- cowplot::plot_grid(p3 + guides(linetype = F, color = F), 
                        p5 + guides(linetype = F), 
                        p6 + guides(color = F),
                        p7 + guides(color = F),
                        labels = c("a", "b", "c", "d"), label_size = 10, nrow = 2)
g <- cowplot::plot_grid(g, legend, ncol = 1, rel_heights = c(1, .1))
g

save_plot("simulated_rt.jpg", p2, dpi = 600)
save_plot("simulated_behavior.jpg", g, dpi = 600)