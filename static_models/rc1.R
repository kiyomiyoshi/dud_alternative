#+ message = F
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(evd)


#'# data loading
dat <- fread("exp1_data_behavior.csv", sep = ",", header = T)
dat <- na.omit(dat)


#'# functions
fit_rc1 <- function(d) {
    guess <- c(0.5, 2.22 * 10e-3)
    fit <- suppressWarnings(optim(par = guess, fn = deviance, gr = NULL, method = "L-BFGS-B", d = d,
                                  lower = c(0.001, 2.22 * 10e-16), upper = c(20, 10),
                                  control = list("maxit" = 100000,
                                                 "parscale" = c(0.004, 0.0002))))
    est <- data.frame(sH = fit$par[1], w = fit$par[2], dev = fit$value)
    return(est)
}

deviance <- function(x, d) {
    val3 <- d$val3
    
    sH <- x[1]
    w <- x[2]
    
    Vcu <- as.data.frame(cbind(100, 90, val3))/100
    colnames(Vcu) <- c("V1", "V2", "V3")
    Vcu <- mutate(Vcu, sumVcu = V1 + V2 + V3)
    Vcu <- mutate(Vcu, normalizer = sH + w * sumVcu)
    M <- Vcu[, 1:3]/Vcu$normalizer
    
    M <- mutate(M, P = ifelse(V3 == 0, exp(V1) / (exp(V1) + exp(V2)), 
                              exp(V1) / (exp(V1) + exp(V2) + exp(V3))))
    
    M <- mutate(M, P = ifelse(P < .01, .01, P)) #adjustments to avoid punishing models too much for very unlikely predictions
    M <- mutate(M, P = ifelse(P > .99, .99, P)) #adjustments to avoid punishing models too much for very unlikely predictions
    M$corr <- d$corr
    
    dev <- -2 * sum(M$corr * log(M$P) + (1 - M$corr) * log(1 - M$P))
    return(dev)
}


#'# individual fittings
fits <- c()

for (sub in unique(dat$subj)) {
    d <- subset(dat, dat$subj == sub)
    f <- tryCatch(fit_rc1(d), error = function(e){f = cbind(NA, NA, NA)}) 
    colnames(f) <- c("sH", "w", "dev")
    fits <- try(rbind(fits, cbind(f, sub)))
    fits <- na.omit(fits)
    fits$sH <- as.numeric(fits$sH)
    fits$w <- as.numeric(fits$w)
    fits$dev <- as.numeric(fits$dev)
}


#'# visualization 
prediction <- c()

for (subj in unique(fits$sub)) {
    f <- subset(fits, fits$sub == subj)
    
    n1 <- exp(1 / (f[1] + f[2] * 1.90)) + exp(0.9 / (f[1] + f[2] * 1.90))
    n2 <- exp(1 / (f[1] + f[2] * 2.17)) + exp(0.9 / (f[1] + f[2] * 2.17)) + exp(0.27 / (f[1] + f[2] * 2.17))
    n3 <- exp(1 / (f[1] + f[2] * 2.35)) + exp(0.9 / (f[1] + f[2] * 2.35)) + exp(0.45 / (f[1] + f[2] * 2.35))
    n4 <- exp(1 / (f[1] + f[2] * 2.53)) + exp(0.9 / (f[1] + f[2] * 2.53)) + exp(0.63 / (f[1] + f[2] * 2.53))
    n5 <- exp(1 / (f[1] + f[2] * 2.66)) + exp(0.9 / (f[1] + f[2] * 2.66)) + exp(0.76 / (f[1] + f[2] * 2.66))
    n6 <- exp(1 / (f[1] + f[2] * 2.76)) + exp(0.9 / (f[1] + f[2] * 2.76)) + exp(0.86 / (f[1] + f[2] * 2.76))           
    
    # softmax transformation                                                                      
    pred <- c(
        exp(1 / (f[1] + f[2] * 1.90)) / n1,
        exp(1 / (f[1] + f[2] * 2.17)) / n2,
        exp(1 / (f[1] + f[2] * 2.35)) / n3,
        exp(1 / (f[1] + f[2] * 2.53)) / n4,
        exp(1 / (f[1] + f[2] * 2.66)) / n5,
        exp(1 / (f[1] + f[2] * 2.76)) / n6)
    
    pred <- as.numeric(pred)
    pred <- as.data.frame(pred)
    pred$id <- subj
    pred$condition <- c(0, 27, 45, 63, 76, 86)
    prediction <- rbind(prediction, pred)
    
}


#'# accuracy
dat %>%
    group_by(subj, val3) %>%
    summarise(Accuracy = mean(corr)) -> acc

acc <- subset(acc, acc$subj %in% unique(fits$sub))
acc$pred <- prediction$pred

rc1 <- ggplot(acc) + geom_point(aes(x = val3, y = Accuracy, color = subj)) +
    geom_line(mapping = aes(x = val3, y = pred, color = subj)) + ylim(0.49, 0.9) + xlab("Condition") + 
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) +
    ggtitle("rc1")  + guides(color = guide_legend(title = NULL))
rc1

fits_rc1 <- fits
fits_rc1 <- mutate(fits_rc1, aic = dev + 4)