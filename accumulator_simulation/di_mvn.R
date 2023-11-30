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

mu_d <-    0.007316
sigma <-   0.02336
tau_nd <-  0.3343
sigma_d <- 0.00015
beta <- 0.25
w <-  0
th <- 1

trial <- 20000
sample <- 3501
condition <- c(92, 96, 100)
simdat <- c()

tic()
for (j in condition) {
    
    dr_scaler <- cbind(rnorm(trial, mu_d, sigma_d),
                       rnorm(trial, mu_d, sigma_d)) 
    a10 <- runif(trial, -beta, beta)
    a20 <- runif(trial, -beta, beta) 
    
    d <- foreach(i = 1:trial, .combine = "rbind", .packages = c("magrittr", "dplyr")) %dopar% {
        
        stimulus <- c(j, 90) / 100
        # stimulus <- stimulus / (1 + w * sum(stimulus))
        stimulus <- stimulus * dr_scaler[i, ] # drift rate variability
        
        accum <- matrix(0, nrow = sample + 1, ncol = 2)
        accum[, 1] <- c(a10[i], rep(stimulus[1], sample))
        accum[, 2] <- c(a20[i], rep(stimulus[2], sample))
        accum[-1, 1] <- accum[-1, 1] + rnorm(sample, 0, sigma)
        accum[-1, 2] <- accum[-1, 2] + rnorm(sample, 0, sigma)
        accum <- matrixStats::colCumsums(accum)
        
        mvn <- matrixStats::rowOrderStats(-accum, which = 2) - matrixStats::rowOrderStats(-accum, which = 1)
        
        simrt  <- which(mvn > th)[1]
        choice <- which.max(accum[simrt, ])
        second <- second <- ifelse(choice == 1, 2, 1)
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
simdat <- na.omit(simdat)
simdat <- mutate(simdat, early_conf = mvn_350 - mvn_0)
simdat <- mutate(simdat, post_conf  = e_chosen_post - e_second_post - mvn_choice)
simdat$choice <- ifelse(simdat$choice == 1, "target", "distractor")
simdat$choice <- factor(simdat$choice, levels = c("target", "distractor"))

# RT histogram
simdat %>%
    complete(choice, condition) %>%
    mutate(rt = rt / 1000) %>%
    ggplot() +
    geom_freqpoly(aes(x = rt)) +
    scale_color_discrete(labels = c("target", "distractor")) +
    labs(fill = "choice") + xlab("RT (s)") + ylab("Count") +
    facet_wrap(. ~ factor(condition) + factor(choice), nrow = 5) +
    scale_color_npg() -> p1
p1

# choice proportion
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
    theme(legend.position = c(.65, .5), legend.background = element_rect(fill = NA, colour = NA)) + xlab("Target value") +
    scale_color_discrete(labels = c("target", "distractor")) + labs(color = NULL) +
    scale_x_continuous(breaks = c(92, 96, 100)) + scale_color_npg() -> p2

# mean RT
simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_rt = mean(rt) / 1000) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_rt, color = factor(choice))) +
    scale_color_discrete(labels = c("target", "distractor")) + labs(color = "choice") +
    xlab("Target value") + ylab("Mean RT (s)") + ylim(0.9, 1.2) + 
    scale_x_continuous(breaks = c(92, 96, 100)) + scale_color_npg() -> p3

# mean confidence
simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_ec = mean(early_conf)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_ec, color = factor(choice))) +
    scale_color_discrete(labels = c("target", "distractor")) + labs(color = "choice") +
    xlab("Target value") + ylab("Early-stage MVN") + 
    scale_x_continuous(breaks = c(92, 96, 100)) + 
    scale_y_continuous(breaks = c(0.35, 0.4, 0.45), limits = c(0.35, 0.45)) +
    scale_color_npg() -> p4

simdat %>%
    group_by(choice, condition) %>%
    summarise(mean_pc = mean(post_conf)) %>%
    ggplot() + geom_line(aes(x = condition, y = mean_pc, color = factor(choice))) +
    scale_color_discrete(labels = c("target", "distractor")) + labs(color = "choice") +
    xlab("Target value") + ylab("Post-decisional MVN") + 
    scale_x_continuous(breaks = c(92, 96, 100)) + scale_color_npg() -> p5


figure_di_mvn <- cowplot::plot_grid(p2 + scale_color_npg() +
                                        theme(legend.text = element_text(size = 6),
                                              legend.position = c(.67, .5),
                                              legend.background = element_rect(fill = NA, color = NA),
                                              legend.key = element_rect(fill = NA, color = NA)),
                                    p3 + guides(color = F),
                                    p4 + guides(color = F),
                                    p5 + guides(color = F),
                                    labels = c("a", "b", "c", "d"), label_size = 10)
figure_di_mvn
sjPlot::save_plot("di_mvn.jpg", figure_di_mvn, dpi = 600)
sjPlot::save_plot("di_rt_mvn.jpg", p1, dpi = 600)
# write.csv(simdat, "simdat_di_mvn.csv")