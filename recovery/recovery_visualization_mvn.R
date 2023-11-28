#+ message = F
library(data.table)
library(tidyverse)
library(sjPlot)
library(cowplot)

#'# 
files <- dir("10", ".*csv$", full.names = TRUE) 
dat <- c()
win <- c()
for (i in 1:length(files)) {
    
    d <- fread(files[i], header = TRUE)
    colnames(d) <- c("iterarion", "mu_d", "sigma", "tau_nd", "theta", "chisq", "conv", "sigma_d", "beta")
    d$id <- i
    dat <- rbind(dat, d)
    win <- rbind(win, d[which.min(d$chisq), ])
    
}

true_value <- data.frame(value = c(0.76, 0.23, 0.025, 0.36, 0.4, 0.18),
                         name = c("mu_d", "sigma", "sigma_d", "tau_nd", "beta", "theta"))
win <- select(win, mu_d, sigma, sigma_d, tau_nd, beta, theta)

g10 <- pivot_longer(win, cols = c("mu_d", "sigma", "sigma_d", "tau_nd", "beta", "theta")) %>%
    ggplot() + geom_jitter(aes(x = name, y = value, color = name), width = 0.2, height = 0) +
    geom_point(true_value, mapping = aes(x = name, y = value), shape = 4, size = 3) +
    scale_x_discrete(limits = c("sigma", "mu_d", "sigma_d", "beta", "tau_nd", "theta"),
                     labels = c(expression(bold(paste(σ))),
                                expression(bold(paste(μ[d]))),
                                expression(bold(paste(σ[d]))),
                                expression(bold(paste(β))),
                                expression(bold(paste(τ[nd]))),
                                expression(bold(paste(θ))))) + 
    ggtitle(expression(paste(italic("SD"), " = mean/10"))) +
    xlab("") + ylab("Estimated value") + ylim(0, 0.8) + guides(color = FALSE) +
    theme(axis.text.x  = element_text(size = 7))

#'# 
files <- dir("7", ".*csv$", full.names = TRUE) 
dat <- c()
win <- c()
for (i in 1:length(files)) {
    
    d <- fread(files[i], header = TRUE)
    colnames(d) <- c("iterarion", "mu_d", "sigma", "tau_nd", "theta", "chisq", "conv", "sigma_d", "beta")
    d$id <- i
    dat <- rbind(dat, d)
    win <- rbind(win, d[which.min(d$chisq), ])
    
}

true_value <- data.frame(value = c(0.76, 0.23, 0.025, 0.36, 0.4, 0.18),
                         name = c("mu_d", "sigma", "sigma_d", "tau_nd", "beta", "theta"))
win <- select(win, mu_d, sigma, sigma_d, tau_nd, beta, theta)

g7 <- pivot_longer(win, cols = c("mu_d", "sigma", "sigma_d", "tau_nd", "beta", "theta")) %>%
    ggplot() + geom_jitter(aes(x = name, y = value, color = name), width = 0.2, height = 0) +
    geom_point(true_value, mapping = aes(x = name, y = value), shape = 4, size = 3) +
    scale_x_discrete(limits = c("sigma", "mu_d", "sigma_d", "beta", "tau_nd", "theta"),
                     labels = c(expression(bold(paste(σ))),
                                expression(bold(paste(μ[d]))),
                                expression(bold(paste(σ[d]))),
                                expression(bold(paste(β))),
                                expression(bold(paste(τ[nd]))),
                                expression(bold(paste(θ))))) + 
    ggtitle(expression(paste(italic("SD"), " = mean/7"))) +
    xlab("") + ylab("Estimated value") + ylim(0, 0.8) + guides(color = FALSE) +
    theme(axis.text.x  = element_text(size = 7))


#'#
g <- cowplot::plot_grid(g10, g7,
                        labels = c("a", "b"), label_size = 10, nrow = 1)
g
save_plot("recovery_mvn.jpg", g, dpi = 600)