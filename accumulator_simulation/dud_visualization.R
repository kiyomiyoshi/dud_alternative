library(doParallel)
library(magrittr)
library(dplyr)
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

trial <- 50000
sample <- 3501
condition <- c(0, 27)
simdat <- c()

tic()
for (j in condition) {
    
    dr_scaler <- cbind(rnorm(trial, mu_d, sigma_d),
                       rnorm(trial, mu_d, sigma_d), 
                       rnorm(trial, mu_d, sigma_d)) 
    
    a10 <- runif(trial, -beta, beta)
    a20 <- runif(trial, -beta, beta)
    a30 <- runif(trial, -beta, beta)
    
    d <- foreach(i = 1:trial, .combine = "rbind", .packages = c("magrittr", "dplyr")) %dopar% {
        
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
        choice <- which.max(accum[simrt, 1:nAlt])
        e1_0 <- max(accum[1, ])
        e2_0 <- sort(accum[1, ], decreasing = T)[2]
        e3_0 <- min(accum[1, ])
        e1_350 <- max(accum[350, ])
        e2_350 <- sort(accum[350, ], decreasing = T)[2]
        e3_350 <- min(accum[350, ])
        as.numeric(c(simrt, choice, mvn[1], e1_0, e2_0, e3_0, mvn[350], e1_350, e2_350, e3_350))
        
    }
    
    d <- as.data.frame(d)
    colnames(d) <- c("rt", "choice", "mvn_0", "e1_0", "e2_0", "e3_0", "mvn_350", "e1_350", "e2_350", "e3_350")
    d$rt   <- d$rt + tau_nd * 1000
    d$condition <- j
    simdat <- rbind(simdat, filter(d, rt <= 3000))
    
}
toc()

# model visualization
g1 <- ggplot(simdat) + geom_density(aes(x = e1_0,  color = factor(condition))) + xlim(-0.22, 0.22) + ylim(0, 7) +
    xlab(expression(bold(paste(E^0, " (first-place)"))))  + 
    labs(color = "Dud value") +
    theme(legend.position = c(0.5, 1.3),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.key = element_rect(fill = NA),
          legend.key.size = unit(0.25, "cm"),
          legend.title = element_text(face = "plain"))
g2 <- ggplot(simdat) + geom_density(aes(x = e2_0,  color = factor(condition))) + xlim(-0.22, 0.22) + ylim(0, 7) +
    xlab(expression(bold(paste(E^0, " (second-place)")))) + theme(legend.position = "none")
g3 <- ggplot(subset(simdat, simdat$condition == 27)) + geom_density(aes(x = e3_0), color = "#00BFC4") + 
    xlab(expression(bold(paste(E^0, " (third-place)"))))  + xlim(-0.22, 0.22) + ylim(0, 7)
g4 <- ggplot(simdat) + geom_density(aes(x = mvn_0,   color = factor(condition))) + ylim(0, 7) +
    xlab(expression(bold(paste(MVN^0)))) + theme(legend.position = "none")

g5 <- ggplot(simdat) + geom_density(aes(x = e1_350,  color = factor(condition))) + xlim(0, 3) + ylim(0, 1.35) +
    xlab(expression(bold(paste(E^350, " (first-place)"))))  + theme(legend.position = "none")
g6 <- ggplot(simdat) + geom_density(aes(x = e2_350,  color = factor(condition))) + xlim(0, 3) + ylim(0, 1.35) +
    xlab(expression(bold(paste(E^350, " (second-place)")))) + theme(legend.position = "none")
g7 <- ggplot(subset(simdat, simdat$condition == 27)) + geom_density(aes(x = e3_350), color = "#00BFC4") + 
    xlab(expression(bold(paste(E^350, " (third-place)"))))  + xlim(0, 3) + ylim(0, 1.35)
g8 <- ggplot(simdat) + geom_density(aes(x = mvn_350,   color = factor(condition))) + xlim(0, 1.5) +
    xlab(expression(bold(paste(MVN^350)))) + theme(legend.position = "none")

g9 <- ggplot(simdat) + geom_density(aes(x = mvn_350 - mvn_0,  color = factor(condition))) + xlim(-0.5, 1.5) +
    xlab("Early confidence variable") + theme(legend.position = "none")

# save figures
dud_effect <- cowplot::plot_grid(g1, g2, g3, g4, g5, g6, g7, g8, g9,
                                 labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i"),
                                 label_size = 10, scale = 1.1)
dud_effect
cowplot::save_plot("figure_s4.jpg", dud_effect, dpi = 600)
# sjPlot::save_plot("figure_s4.jpg", dud_effect, dpi = 600)