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
    
    d <- foreach(i = 1:trial, .combine = "rbind", .packages = c("dplyr", "matrixStats")) %dopar% {
        
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
        
        mvn <- matrixStats::rowOrderStats(-accum, which = 2) - matrixStats::rowOrderStats(-accum, which = 1)
        simrt  <- which(mvn > th)[1]
        choice <- which.max(accum[simrt, ])
        third  <- which.max(-accum[simrt, ])
        second <- ifelse(nAlt == 3, setdiff(c(1, 2, 3), c(choice, third)), third)
        # pre <- ifelse(simrt <= 350, simrt - 1, 350)
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
df <- read.csv("exp2_data_3afc_processed.csv", header = TRUE)
df <- na.omit(df)
df$rt <- as.numeric(df$rt)
df$chosenItem <- factor(df$chosenItem, levels = c("target", "distractor", "dud"))
df <- mutate(df, confidence = ifelse(conf > 2, 1, 0))

simdat <- na.omit(simdat)
simdat <- mutate(simdat, early_conf = mvn_350 - mvn_0)
simdat <- mutate(simdat, post_conf  = e_chosen_post - e_second_post - mvn_choice)
simdat <- mutate(simdat, cv = 0.4 * early_conf + 0.6 * post_conf)
simdat <- mutate(simdat, confidence = ifelse(cv > theta, "high confidence", "low confidence"))
simdat$choice <- ifelse(simdat$choice == 1, "target", 
                        ifelse(simdat$choice == 2, "distractor", "dud"))
simdat$choice <- factor(simdat$choice, levels = c("target", "distractor", "dud"))
simdat <- rename(simdat, chosenItem = choice, val3 = condition)


# RT histogram
df %>%
    tidyr::complete(chosenItem, val3) %>%
    dplyr::select(chosenItem, val3, rt) %>%
    mutate(rt = rt * 1000) -> emp1
emp1$chosenItem <- factor(emp1$chosenItem, levels = c("target", "distractor", "dud"))

simdat %>%
    complete(chosenItem, val3) %>%
    ggplot() +
    geom_histogram(emp1, mapping = (aes(x = rt, fill = chosenItem, color = chosenItem))) +
    geom_freqpoly(aes(x = rt, y = after_stat(count / nrow(simdat) * nrow(df)))) +
    scale_color_discrete(labels = c("target", "distractor", "dud")) +
    guides(fill = F, color = F) + xlab("RT (s)") + ylab("Count") +
    scale_x_continuous(breaks = c(0, 1000, 2000, 3000), labels = c("0", "1", "2", "3")) +
    facet_wrap(. ~ factor(val3) + factor(chosenItem), nrow = 3) +
    theme(axis.text = element_text(size = 6.5), axis.title = element_text(size = 7)) +
    scale_fill_npg() + scale_color_npg() -> p1
p1

# choice proportion
simdat %>%
    dplyr::select(chosenItem, val3) %>%
    group_by(chosenItem, val3) %>%
    mutate(n1 = n()) %>%
    ungroup(chosenItem, val3) %>%
    group_by(val3) %>%
    mutate(n2 = n()) %>%
    dplyr::select(n1, n2, chosenItem, val3) %>%
    distinct() %>%
    mutate(p = n1/n2) %>%
    ungroup(chosenItem, val3) %>%
    complete(chosenItem, val3) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_all(~replace(., is.na(.), 0)) -> m2
m2$chosenItem <- factor(m2$chosenItem, levels = c("target", "distractor", "dud"))

df %>%
    group_by(chosenItem, val3, subj) %>%
    mutate(n1 = n()) %>%
    ungroup(chosenItem, val3, subj) %>%
    group_by(val3, subj) %>%
    mutate(n2 = n()) %>%
    dplyr::select(n1, n2, val3, chosenItem, subj) %>%
    distinct() %>%
    mutate(p = n1/n2) %>%
    ungroup(chosenItem, val3, subj) %>%
    complete(chosenItem, val3, subj) %>%
    mutate_all(~ replace_na(., 0)) %>%
    ggplot() + 
    geom_line(m2, size = 0.4, mapping = aes(x = val3, y = p, color = chosenItem)) +
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = p, color = chosenItem)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = p, color = chosenItem)) +
    labs(color = "choice") + xlab("Dud value") +
    theme(legend.position = c(.5, .6), legend.text = element_text(size = 10), 
          legend.background = element_rect(fill = NA, colour = NA), legend.key = element_rect(fill = NA)) + 
    guides(color = guide_legend(title = NULL)) + scale_color_npg() +
    coord_cartesian(ylim = c(0, 0.82)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylab("Choice proportion") -> p2
p2

# target vs. distractor choice
simdat %>%
    dplyr::select(chosenItem, val3) %>%
    group_by(val3) %>%
    summarise(n1 = sum(chosenItem == "target"), n2 = sum(chosenItem == "distractor"), n3 = sum(chosenItem == "dud")) %>%
    mutate(r = n1/n2) -> m3

df %>%
    group_by(val3, subj) %>%
    summarise(n1 = sum(chosenItem == "target"), n2 = sum(chosenItem == "distractor"), n3 = sum(chosenItem == "dud")) %>%
    mutate(r = n1/n2) %>%
    ggplot() + 
    geom_point(aes(x = val3, y = r)) +
    geom_line(m3, size = 0.4, mapping = aes(x = val3, y = r)) +
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = r)) + xlab("Dud value") +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = r)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylab("Pc(target) / Pc(distractor)") -> p3
p3

subset(df, df$subj != "subj21") %>%
    group_by(val3, subj) %>%
    summarise(n1 = sum(chosenItem == "target"), n2 = sum(chosenItem == "distractor"), n3 = sum(chosenItem == "dud")) %>%
    mutate(r = n1/n2) %>%
    ggplot() + 
    geom_line(m3, size = 0.4, mapping = aes(x = val3, y = r)) +
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = r)) + xlab("Dud value") +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = r)) + 
    coord_cartesian(ylim = c(3.5, 5.5)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylab("Pc(target) / Pc(distractor)") -> p3
p3

# mean RT
simdat %>%
    dplyr::select(rt, chosenItem, val3) %>%
    group_by(chosenItem, val3) %>%
    summarise(mean_rt = mean(rt) / 1000) -> m4
m4$chosenItem <- factor(m4$chosenItem, levels = c("target", "distractor", "dud"))

df %>% group_by(val3, chosenItem, subj) %>%
    summarise(mean_rt = mean(rt)) %>%
    ggplot() + 
    geom_line(m4, size = 0.4, mapping = aes(x = val3, y = mean_rt, color = chosenItem)) +
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = mean_rt, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = mean_rt, color = chosenItem)) +
    guides(color = F) + xlab("Dud value") + ylab("Mean RT (s)") + scale_color_npg() +
    scale_y_continuous(breaks = seq(0.7, 2.1, 0.3)) + coord_cartesian(ylim = c(0.7, 2.1)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_wrap(chosenItem ~ .) -> p4
p4

# high confidence rate
simdat %>%
    dplyr::select(confidence, chosenItem, val3) %>%
    mutate(conf = ifelse(confidence == "high confidence", 1, 0)) %>%
    group_by(chosenItem, val3) %>%
    summarise(high_conf_rate = mean(conf)) -> m5
m5$chosenItem <- factor(m5$chosenItem, levels = c("target", "distractor", "dud"))

df %>% group_by(val3, chosenItem, subj) %>%
    summarise(mean_confidence = mean(conf)) %>%
    ggplot() + 
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = mean_confidence, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = mean_confidence, color = chosenItem)) +
    xlab("Dud value") + ylab("Mean confidence") + guides(color = F) + scale_color_npg() +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_wrap(chosenItem ~ .) -> p5
p5

df %>% group_by(val3, chosenItem, subj) %>%
    summarise(high_conf = mean(confidence)) %>%
    ggplot() + 
    geom_line(m5, size = 0.4, mapping = aes(x = val3, y = high_conf_rate, color = chosenItem)) +
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = high_conf, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = high_conf, color = chosenItem)) +
    xlab("Dud value") + ylab("High confidence rate") + guides(color = F) + scale_color_npg() +
    scale_y_continuous(breaks = seq(0, 0.8, 0.2)) + coord_cartesian(ylim = c(0, 0.8)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_wrap(chosenItem ~ .) -> p5
p5

# rt for low confidence responses
simdat %>% group_by(chosenItem, val3, confidence) %>%
    summarise(mean_rt = mean(rt) / 1000) %>% 
    filter(confidence == "low confidence") -> m6
m6$chosenItem <- factor(m6$chosenItem, levels = c("target", "distractor", "dud"))

df <- mutate(df, confidence = ifelse(confidence == 1, "high confidence", "low confidence"))
df %>% group_by(val3, chosenItem, subj, confidence) %>%
    summarise(mean_rt = mean(rt)) %>%
    filter(confidence == "low confidence") %>%
    ggplot() + 
    geom_line(m6, size = 0.4, mapping = aes(x = val3, y = mean_rt, color = chosenItem)) +
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = mean_rt, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = mean_rt, color = chosenItem)) +
    xlab("Dud value") + ylab("Mean RT (s)") + guides(color = F) + scale_color_npg() +
    scale_y_continuous(breaks = seq(0.7, 2.1, 0.3)) + coord_cartesian(ylim = c(0.7, 2.1)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_nested(. ~ confidence + chosenItem) -> p6
p6

# rt for high confidence responses
simdat %>% group_by(chosenItem, val3, confidence) %>%
    summarise(mean_rt = mean(rt) / 1000) %>% 
    filter(confidence == "high confidence") -> m7
m7$chosenItem <- factor(m7$chosenItem, levels = c("target", "distractor", "dud"))

df %>% group_by(val3, chosenItem, subj, confidence) %>%
    summarise(mean_rt = mean(rt)) %>%
    filter(confidence == "high confidence") %>%
    ggplot() + 
    geom_line(m7, size = 0.4, mapping = aes(x = val3, y = mean_rt, color = chosenItem)) +
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = mean_rt, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = mean_rt, color = chosenItem)) +
    xlab("Dud value") + ylab("Mean RT (s)") + guides(color = F) + scale_color_npg() +
    scale_y_continuous(breaks = seq(0.7, 2.1, 0.3)) + coord_cartesian(ylim = c(0.7, 2.1)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_nested(. ~ confidence + chosenItem) -> p7
p7

# early confidence
simdat %>%
    group_by(chosenItem, val3) %>%
    summarise(mean_ec = mean(early_conf)) %>%
    ggplot() + geom_line(aes(x = val3, y = mean_ec, color = factor(chosenItem))) +
    scale_color_discrete(labels = c("target", "distractor", "dud")) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + 
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) + coord_cartesian(ylim = c(-0.45, 0.4)) +
    theme(legend.position = c(.37, .46), legend.text = element_text(size = 10), 
          legend.background = element_rect(fill = NA, colour = NA), legend.key = element_rect(fill = NA)) + 
    xlab("Dud value") + ylab("Early-stage MVN") + 
    guides(color = guide_legend(title = NULL, nrow = 3, keyheight = unit(0.5, "lines"))) +
    scale_color_npg() -> p8

# post confidence
simdat %>%
    group_by(chosenItem, val3) %>%
    summarise(mean_pc = mean(post_conf)) %>%
    ggplot() + geom_line(aes(x = val3, y = mean_pc, color = factor(chosenItem))) +
    scale_color_discrete(labels = c("target", "distractor", "dud")) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + 
    xlab("Dud value") + ylab("Post-decisional MVN") + scale_color_npg() -> p9

# save figures
figure_2 <- cowplot::plot_grid(p2 + theme(legend.text = element_text(size = 6), 
                                          legend.background = element_rect(fill = NA, colour = NA),
                                          legend.key = element_rect(fill = NA),
                                          plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                               p3 + theme(plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                               p4 + theme(axis.text.x = element_text(angle = 90),
                                          plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                               p5 + theme(axis.text.x = element_text(angle = 90),
                                          plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")), 
                               p6 + theme(axis.text.x = element_text(angle = 90),
                                          plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                               p7 + theme(axis.text.x = element_text(angle = 90),
                                          plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                               labels = c("a", "b", "c", "d", "e", "f"), label_size = 10, nrow = 3)

figure_5 <- cowplot::plot_grid(p8 + theme(legend.text = element_text(size = 6)),
                               p9 + guides(color = F),
                               labels = c("a", "b"), label_size = 10, nrow = 1)
figure_2
figure_5

sjPlot::save_plot("figure_2.jpg", figure_2, dpi = 600, height = 13, width = 14)
sjPlot::save_plot("figure_5.jpg", figure_5, dpi = 600, height = 4.5, width = 9)
sjPlot::save_plot("figure_A2.jpg", p1, dpi = 600)
write.csv(simdat, "simdat_mvn.csv")

# Loss calculation
rt_q <- c(0, quantile(df$rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9)), 10) # quantiles of grand data
cond <- c(0, 27, 45, 63, 76, 86) # response frequency table
tab_low_list <- c()
tab_high_list <- c()
for (i in 1:6) {
    dat_low <-  subset(df, df$conf < 3 & df$val3 == cond[i])
    dat_high <- subset(df, df$conf > 2 & df$val3 == cond[i])
    cuts_low <-  factor(cut(dat_low$rt,  rt_q, labels = FALSE), levels = 1:(length(rt_q) - 1))
    cuts_high <- factor(cut(dat_high$rt, rt_q, labels = FALSE), levels = 1:(length(rt_q) - 1))
    tab_low <- matrix(c(table(cuts_low[dat_low$chosenItem == "target"]),
                        table(cuts_low[dat_low$chosenItem == "distractor"]),
                        table(cuts_low[dat_low$chosenItem == "dud"])), nrow = 3, byrow = TRUE) + 1 # avoiding zero cells
    tab_high <- matrix(c(table(cuts_high[dat_high$chosenItem == "target"]),
                         table(cuts_high[dat_high$chosenItem == "distractor"]),
                         table(cuts_high[dat_high$chosenItem == "dud"])), nrow = 3, byrow = TRUE) + 1 # avoiding zero cells
    tab_low_list[[i]] <-  tab_low
    tab_high_list[[i]] <- tab_high
}

chisq_stats <- c()
qll_stats <- c()
for (j in 1:6) {
    k <- trial / 6624
    sub <- subset(simdat, simdat$val3 == cond[j])
    sub <- sub %>% dplyr::select(rt, chosenItem, cv)
    sub$rt <- sub$rt / 1000
    sim_low <-  sub[sub[, 3] <= theta, ]
    sim_high <- sub[sub[, 3] >  theta, ]
    simcuts_low <-  factor(cut(sim_low[, 1],  rt_q, labels = FALSE), levels = 1:(length(rt_q) - 1))
    simcuts_high <- factor(cut(sim_high[, 1], rt_q, labels = FALSE), levels = 1:(length(rt_q) - 1))
    tab_sim_low <- matrix(c(table(simcuts_low[sub[, 2] == "target"]), 
                            table(simcuts_low[sub[, 2] == "distractor"]), 
                            table(simcuts_low[sub[, 2] == "dud"])), nrow = 3, byrow = TRUE) / k + 1 # zero cell correction
    tab_sim_high <- matrix(c(table(simcuts_high[sub[, 2] == "target"]), 
                             table(simcuts_high[sub[, 2] == "distractor"]), 
                             table(simcuts_high[sub[, 2] == "dud"])), nrow = 3, byrow = TRUE) / k + 1 # zero cell correction
    chisq_cond <- (tab_sim_low - tab_low_list[[j]]) ^ 2 /  tab_sim_low + 
        (tab_sim_high - tab_high_list[[j]]) ^ 2 /  tab_sim_high
    chisq_cond <- chisq_cond[is.finite(chisq_cond)]
    chisq_stats <- c(chisq_stats, sum(chisq_cond))
    n <- sum(tab_sim_low) + sum(tab_sim_high)
    p_sim_low  <- tab_sim_low / n
    p_sim_high <- tab_sim_high / n
    qll_stats <- c(qll_stats, sum(log(p_sim_low) * tab_low_list[[j]]) + sum(log(p_sim_high) * tab_high_list[[j]]))
}

chisq_stats
sum(chisq_stats)
qll_stats
sum(qll_stats)
-2 * sum(qll_stats) + 2 * 7