library(doParallel)
library(tictoc)
library(matrixStats)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)


# import data
df <- read.csv("exp1_data_behavior.csv", header = TRUE)
df <- na.omit(df)
df$chosenItem <- ifelse(df$chosenItem == "target", 1, 
                        ifelse(df$chosenItem == "distractor", 2, 3))


# model fitting function
fit_ppa <- function(dat_behav, sigma_d, beta){
    init_par <- c(0.15, 0.55, 0.3, 0.32)
    fit <- suppressWarnings(optim(par = init_par, fn = chisq, gr = NULL, 
                                  method = "L-BFGS-B", dat1 = dat_behav, sigma_d = sigma_d, beta = beta,
                                  lower = c(0.11, 0.5, 0.2, 0.26), upper = c(0.19, 0.6, 0.4, 0.38),
                                  control = list("maxit" = 3000, "parscale" = c(1, 0.4, 1, 0.4))))
    est <- data.frame(mu_d = fit$par[1], sigma = fit$par[2], tau_nd = fit$par[3], theta = fit$par[4], chisq = fit$value, conv = fit$convergence)
    return(est)
}


# objective function
chisq <- function(par, dat1, sigma_d, beta) {
    
    mu_d    <- par[1]
    sigma <-   par[2]
    tau_nd  <- par[3]
    theta <-   par[4]
    # w     <- 0
    
    th <- -0.1
    sample <- 3351
    
    rt_q <- c(0, quantile(dat1$rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9)), 10) # quantiles of grand data
    
    cond <- c(0, 27, 45, 63, 76, 86) # response frequency table
    tab_low_list <- c()
    tab_high_list <- c()
    for (i in 1:6) {
        dat_low <-  subset(dat1, dat1$conf < 3 & dat1$val3 == cond[i])
        dat_high <- subset(dat1, dat1$conf > 2 & dat1$val3 == cond[i])
        cuts_low <-  factor(cut(dat_low$rt,  rt_q, labels = FALSE), levels = 1:(length(rt_q) - 1))
        cuts_high <- factor(cut(dat_high$rt, rt_q, labels = FALSE), levels = 1:(length(rt_q) - 1))
        tab_low <- matrix(c(table(cuts_low[dat_low$chosenItem == 1]),
                            table(cuts_low[dat_low$chosenItem == 2]),
                            table(cuts_low[dat_low$chosenItem == 3])), nrow = 3, byrow = TRUE) + 1 # avoiding zero cells
        tab_high <- matrix(c(table(cuts_high[dat_high$chosenItem == 1]),
                             table(cuts_high[dat_high$chosenItem == 2]),
                             table(cuts_high[dat_high$chosenItem == 3])), nrow = 3, byrow = TRUE) + 1 # avoiding zero cells
        tab_low_list[[i]] <-  tab_low
        tab_high_list[[i]] <- tab_high
    }
    
    chisq_stats <- c()
    
    for (j in 1:6) {
        
        dr_scaler <- cbind(rnorm(2880 * 20, mu_d, sigma_d),
                           rnorm(2880 * 20, mu_d, sigma_d),
                           rnorm(2880 * 20, mu_d, sigma_d))
        
        a10 <- runif(2880 * 20, -beta, beta)
        a20 <- runif(2880 * 20, -beta, beta)
        a30 <- runif(2880 * 20, -beta, beta)
        
        d <- foreach(i = 1:(2880 * 20), .combine = "rbind", .packages = c("matrixStats")) %dopar% {
            stimulus <- c(100, 90, cond[j]) / 1000
            # stimulus <- stimulus / (1 + w * sum(stimulus))
            stimulus <- stimulus * dr_scaler[i, ] # drift rate variability
            nAlt <- sum(stimulus != 0)
            
            accum <- matrix(0, nrow = sample + 1, ncol = 3)
            accum[, 1] <- c(a10[i], rep(stimulus[1], sample))
            accum[, 2] <- c(a20[i], rep(stimulus[2], sample))
            accum[, 3] <- c(a30[i], rep(stimulus[3], sample))
            accum[-1, ] <- accum[-1, ] + matrix(rnorm(sample * 3, 0, sigma / 10), ncol = 3)
            accum <- accum[, 1:nAlt]
            accum <- matrixStats::colCumsums(accum)
            
            a <- matrixStats::rowMaxs(accum) - log(rowSums(exp(accum)))
            simrt  <- which(ppa > th)[1]
            choice <- which.max(accum[simrt, ])
            ppa_post <- ifelse(simrt <= 3000, accum[simrt + 350, choice] - log(sum(exp(accum[simrt + 350, ]))), NA) # choice-conditional post ppa
            # cv <- 0.5 * (ppa[350] - ppa[1]) + 0.5 * (ppa_post - ppa[simrt])
            cv <- 0.6 * (qnorm(exp(ppa[350])) - qnorm(exp(ppa[1]))) + 0.4 * (qnorm(exp(ppa_post)) - qnorm(exp(ppa[simrt])))
            c(simrt, choice, cv)
        }
        
        d <- na.omit(d)
        d[, 1] <- d[, 1] / 1000 + tau_nd
        d <- d[d[, 1] <= 3, ]
        sim_low <-  d[d[, 3] <= theta, ]
        sim_high <- d[d[, 3] >  theta, ]
        simcuts_low <-  factor(cut(sim_low[, 1],  rt_q, labels = FALSE), levels = 1:(length(rt_q) - 1))
        simcuts_high <- factor(cut(sim_high[, 1], rt_q, labels = FALSE), levels = 1:(length(rt_q) - 1))
        tab_sim_low <- matrix(c(table(simcuts_low[d[, 2] == 1]), 
                                table(simcuts_low[d[, 2] == 2]), 
                                table(simcuts_low[d[, 2] == 3])), nrow = 3, byrow = TRUE) / 20 + 1 # zero cell correction
        tab_sim_high <- matrix(c(table(simcuts_high[d[, 2] == 1]), 
                                 table(simcuts_high[d[, 2] == 2]), 
                                 table(simcuts_high[d[, 2] == 3])), nrow = 3, byrow = TRUE) / 20 + 1 # zero cell correction
        chisq_cond <- (tab_sim_low - tab_low_list[[j]]) ^ 2 /  tab_sim_low + 
            (tab_sim_high - tab_high_list[[j]]) ^ 2 /  tab_sim_high
        chisq_cond <- chisq_cond[is.finite(chisq_cond)]
        chisq_stats <- c(chisq_stats, sum(chisq_cond))
    }
    
    return(sum(chisq_stats))
    
}

# model fit
tic()

result_ppa <- c()
sigma_d <- c(0.001, 0.004, 0.008)
beta <- c(0.1, 0.25, 0.4)
for (k in sigma_d) {
    for (l in beta) {
        fit <- fit_ppa(dat_behav = df, sigma_d = k, beta = l)
        result_ppa <- rbind(result_ppa, cbind(fit, k, l))
        write.csv(result_ppa, "result_ppa.csv")
    }
}

toc()