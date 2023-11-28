library(magrittr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggsci)

theme_set(theme_publication()) 

dr_scaler <- 0.0051
sigma <-     0.022
tau_nd <-    0.23
beta <-      0.25
# w <-    -0.09
th <- 1

sample <- 3501
simdat <- c()

stimulus <- c(100, 90, 45) / 100
# stimulus <- stimulus / (1 + w * sum(stimulus))
stimulus <- stimulus * dr_scaler
nAlt <- sum(stimulus != 0)

accum <- data.frame(a1 = c(runif(1, -beta, beta), rep(stimulus[1], sample)), 
                    a2 = c(runif(1, -beta, beta), rep(stimulus[2], sample)), 
                    a3 = c(runif(1, -beta, beta), rep(stimulus[3], sample)))
accum[-1, 1] <- accum[-1, 1] + rnorm(sample, 0, sigma)
accum[-1, 2] <- accum[-1, 2] + rnorm(sample, 0, sigma)
accum[-1, 3] <- accum[-1, 3] + rnorm(sample, 0, sigma)
accum <- accum[, 1:nAlt]
accum <- cumsum(accum)
accum <- mutate(accum, mvn = apply(accum, 1, max) - apply(accum, 1, function(x){return(max(x[-which.max(x)]))}))
accum <- mutate(accum, max = apply(accum, 1, max),
                inhib = apply(accum, 1, function(x){return(log(sum(c(exp(c(x))))))}),
                mvaa = max - inhib)

accum <- mutate(accum, t = row_number())
accum <- dplyr::select(accum, a1, a2, a3, mvn, mvaa, t)
colnames(accum) <- c("E1", "E2", "E3", "MVN", "MVAA", "t")
accum <- pivot_longer(accum, names_to = "accumulator", values_to = "value", cols = c("E1", "E2", "E3"))
choice = accum$t[which(accum$MVN > 1)[1]]

g1 <- ggplot(accum) + 
    geom_line(aes(x = t, y = value, color = accumulator), size = 0.3) + ylab("Evidence") +
    annotate("segment", x = 0,            xend = 0,            y = 0, yend = 20, linetype = "dashed", size = 0.3) +
    annotate("segment", x = 350,          xend = 350,          y = 0, yend = 20, linetype = "dashed", size = 0.3) +
    annotate("segment", x = choice,       xend = choice,       y = 0, yend = 20, linetype = "dashed", size = 0.3) +
    annotate("segment", x = choice + 350, xend = choice + 350, y = 0, yend = 20, linetype = "dashed", size = 0.3) +
    annotate("rect", xmin = 0,      xmax = 350,          ymin = 0, ymax = 20, alpha = 0.2, fill = "grey") +
    annotate("rect", xmin = choice, xmax = choice + 350, ymin = 0, ymax = 20, alpha = 0.2, fill = "grey") +
    annotate("text", x = 200,          y = 24, label = "early-stage", size = 1.7) +
    annotate("text", x = 200,          y = 22, label = "confidence", size = 1.7) +
    annotate("text", x = choice + 175, y = 24, label = "post-decisional", size = 1.7) +
    annotate("text", x = choice + 175, y = 22, label = "confidence", size = 1.7) +
    xlim(0, 3000) + ylim(0, 25) + labs(color = NULL) +
    scale_color_npg(labels = c(expression(bold(paste({E[tar]}))),
                               expression(bold(paste({E[dis]}))),
                               expression(bold(paste({E[dud]}))))) +
    theme(legend.text = element_text(size = 5),
          legend.position = c(.5, 1.1),
          legend.background = element_rect(fill = NA, color = NA),
          legend.key = element_rect(fill = NA, color = NA))
g1
g2 <- ggplot(accum) + geom_line(aes(x = t, y = MVN), size = 0.8, color = "#3C5488FF") + 
    xlim(0, 3000) + ylim(0, 2.5) + geom_hline(yintercept = 1, linetype = "dashed")

g3 <- ggplot(accum) + geom_line(aes(x = t, y = MVAA), size = 0.8, color = "#F39B7FFF") + 
    xlim(0, 3000) + ylim(-1.5, 0) + geom_hline(yintercept = -0.14, linetype = "dashed") 

fig <- cowplot::plot_grid(g1 + scale_color_npg(labels = c(expression(bold(paste({E[tar]}))),
                                                          expression(bold(paste({E[dis]}))),
                                                          expression(bold(paste({E[dud]}))))) +
                              theme(legend.text = element_text(size = 6),
                                    legend.position = c(.35, .85),
                                    legend.background = element_rect(fill = NA, color = NA),
                                    legend.key = element_rect(fill = NA, color = NA)),
                          g2,
                          labels = c("a", "b"), label_size = 10)
fig
save_plot("accumulator.jpg", fig, dpi = 600, height = 4.498)
save_plot("figure_2.jpg", g1, dpi = 600, width = 6, height = 4.5)
# write.csv(accum, "model_visualization.csv")