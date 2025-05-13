#+ message = F
library(tidyverse)
library(cowplot)
library(ggsci)
library(psych)

theme_set(theme_publication()) 


#'# model fitting
#+ message = F
#' logit: multinomial logit without relative coding  
#' rc1:   relative coding under 0 < w < 10  
#' rc2:   relative coding under -5 < w < 5  

source("static_relative_coding_model/logit.R")
source("static_relative_coding_model/rc1.R")
source("static_relative_coding_model/rc2.R")

fits_logit$model <- "logit"
fits_rc1$model   <- "rc1"
fits_rc2$model   <- "rc2"


#'# fitting curves
#+ message = F
g1 <- cowplot::plot_grid(logit + theme(legend.position = "none"),
                rc1 + theme(legend.position = "none"),
                rc2 + theme(legend.position = "none"), 
                get_legend(logit + guides(color = guide_legend(ncol = 2))),
                labels = c(NA, NA, NA, NA), align = "vh")
g1


#'# model comparison
#+ message = F
fittings <- as.data.frame(rbind(fits_logit, fits_rc1, fits_rc2))
fittings %>%
    group_by(model) %>%
    summarise(summed_aic = sum(aic)) %>%
    ggplot() + geom_point(aes(x = model, y = summed_aic), size = 2) + 
    scale_x_discrete(limits = c("logit", "rc1", "rc2")) + ylab("Summed AIC") + ylim(18100, 18150) -> g2


#'# parameter distribution
#+ message = F
# rc2
ttest <- t.test(fits_rc2$w, mu = 0)
ttest
cd <- abs(ttest$statistic) * sqrt(1 / 10) # cohen's d
cd
cohen.d.ci(cd, n1 = 10, alpha = .05)

counts <- rbind(fits_rc1, fits_rc2)
ggplot(counts) + geom_histogram(aes(x = w, fill = model), alpha = 0.4, position = "identity") + 
    labs(fill = "") + theme(legend.background = element_rect(fill = NA, colour = NA)) +
    xlim(-0.1, 0.1) + ylim(0, 10) + scale_y_continuous(breaks = seq(0, 10, 2)) + xlab("ω") +
    theme(legend.position = c(0.82, 0.85), legend.key.size = unit(0.2, 'cm')) + ylab("Count") -> g3

cowplot::plot_grid(logit + theme(legend.position = "none"),
                   rc1   + theme(legend.position = "none"),
                   rc2   + theme(legend.position = "none"), 
                   g2, g3,
                   labels = c("a", "b", "c", "d", "e"), align = "vh")


#'# save figures
fittings <- as.data.frame(rbind(fits_logit, fits_rc2))
fittings %>%
    group_by(model) %>%
    summarise(summed_aic = sum(aic)) %>%
    ggplot() + geom_point(aes(x = model, y = summed_aic), size = 2) + 
    scale_x_discrete(limits = c("rc2", "logit"),
                     labels = c("relative coding", "standard logit")) + 
    scale_y_continuous(breaks = c(18110, 18115, 18120), limits = c(18110, 18120)) +
    xlab("") + ylab("Summed AIC") -> g4
g4

ggplot(fits_rc2) + geom_histogram(aes(x = w), position = "identity") + 
    labs(fill = "") + theme(legend.background = element_rect(fill = NA, colour = NA)) +
    xlim(-0.1, 0.1) + ylim(0, 10) + scale_y_continuous(breaks = seq(0, 10, 2)) + xlab("ω") +
    theme(legend.position = c(0.82, 0.85), legend.key.size = unit(0.2, 'cm')) + ylab("Count") -> g5
g5

relative_coding <- cowplot::plot_grid(rc2   + theme(legend.position = "none") + ggtitle("relative coding") + 
                                          xlab("Dud value") + scale_color_npg(),
                                      logit + theme(legend.position = "none") + ggtitle("standard logit") + 
                                          xlab("Dud value") + scale_color_npg(),
                                      g4 + scale_color_npg() + theme(axis.text.x = element_text(face = "bold", size = 6.5, angle = 0)), 
                                      g5 + scale_color_npg(),
                                      labels = c("a", "b", "c", "d"))
relative_coding
save_plot("relative_coding.jpg", relative_coding, dpi = 600)