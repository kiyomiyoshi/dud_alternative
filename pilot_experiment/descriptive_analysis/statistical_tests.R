#+ message = F
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(lme4)
library(lmerTest)
library(car)
library(multcomp)
library(sjPlot)
library(ggsignif)
library(cowplot)
library(ggsci)
library(ggh4x)

theme_set(theme_publication()) 


#'# data loading
d <- fread("exp1_data_behavior.csv", sep = ",", header = T)
d$rt <- as.numeric(d$rt)
d <- na.omit(d)
d$chosenItem <- factor(d$chosenItem, levels = c("target", "distractor", "dud"))
d$fVal3 <- as.factor(d$val3)
d <- mutate(d, confidence = ifelse(conf > 2, 1, 0))


#'# data visualization
d %>%
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
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = p, color = chosenItem)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = p, color = chosenItem)) +
    labs(color = "choice") + xlab("Dud value") +
    theme(legend.position = c(.5, .6), legend.text=element_text(size = 10), 
          legend.background = element_rect(fill = NA, colour = NA), legend.key = element_rect(fill = NA)) + 
    guides(color = guide_legend(title = NULL)) + scale_color_npg() +
    coord_cartesian(ylim = c(0, 0.82)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylab("Choice proportion") -> p1
p1

d %>%
    group_by(val3, subj) %>%
    summarise(n1 = sum(chosenItem == "target"), n2 = sum(chosenItem == "distractor"), n3 = sum(chosenItem == "dud")) %>%
    mutate(r = n1/n2) %>%
    ggplot() + 
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = r)) + xlab("Dud value") +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = r)) + 
    coord_cartesian(ylim = c(3, 5)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + ylab("Pc(target) / Pc(distractor)") -> p2
p2

d %>% group_by(val3, chosenItem, subj) %>%
    summarise(mean_rt = mean(rt)) %>%
    ggplot() + 
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = mean_rt, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = mean_rt, color = chosenItem)) +
    guides(color = F) + xlab("Dud value") + ylab("Mean RT (s)") + scale_color_npg() +
    scale_y_continuous(breaks=seq(0.4, 1.8, 0.3)) + coord_cartesian(ylim = c(0.4, 1.8)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_wrap(chosenItem ~ .) -> p3
p3

d %>% group_by(val3, chosenItem, subj) %>%
    summarise(high_conf = mean(confidence)) %>%
    ggplot() + 
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = high_conf, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = high_conf, color = chosenItem)) +
    xlab("Dud value") + ylab("High confidence rate") + guides(color = F) + scale_color_npg() +
    scale_y_continuous(breaks = seq(0, 0.8, 0.2)) + coord_cartesian(ylim = c(0, 0.8)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_wrap(chosenItem ~ .) -> p4
p4

d <- mutate(d, confidence = ifelse(confidence == 1, "high confidence", "low confidence"))
d %>% group_by(val3, chosenItem, subj, confidence) %>%
    summarise(mean_rt = mean(rt)) %>%
    filter(confidence == "low confidence") %>%
    ggplot() + 
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = mean_rt, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = mean_rt, color = chosenItem)) +
    xlab("Dud value") + ylab("Mean RT (s)") + guides(color = F) + scale_color_npg() +
    scale_y_continuous(breaks=seq(0.4, 1.8, 0.3)) + coord_cartesian(ylim = c(0.4, 1.8)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_nested(. ~ confidence + chosenItem) -> p5
p5

d %>% group_by(val3, chosenItem, subj, confidence) %>%
    summarise(mean_rt = mean(rt)) %>%
    filter(confidence == "high confidence") %>%
    ggplot() + 
    stat_summary(fun = "mean", geom = "point", aes(x = val3, y = mean_rt, color = chosenItem)) + 
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(x = val3, y = mean_rt, color = chosenItem)) +
    xlab("Dud value") + ylab("Mean RT (s)") + guides(color = F) + scale_color_npg() +
    scale_y_continuous(breaks=seq(0.4, 1.8, 0.3)) + coord_cartesian(ylim = c(0.4, 1.8)) +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) + facet_nested(. ~ confidence + chosenItem) -> p6
p6


#'# statistical testing
# factorial analysis on accuracy
f1 <- glmer(corr ~ fVal3 + (1 + fVal3|subj), data = d, family = binomial)
summary(f1)
Anova(f1)
# summary(glht(f1, mcp(fVal3 = "Tukey")))
ap <- p.adjust(as.data.frame(summary(f1)[10])[2:6, 4], method = "holm") # adjusted p values for correct choices
pv <- data.frame(x = c(0, 0, 0, 0, 0),
                 y = c(0.84, 0.85, 0.86, 0.87, 0.88),
                 xend = c(27, 45, 63, 76, 86),
                 yend = c(0.84, 0.85, 0.86, 0.87, 0.88),
                 values = (as.character(round(ap, 3))))
g1 <- plot_model(f1, type = "pred", terms = "fVal3", pred.type = "fe") + 
    xlab("Dud value") + ylab("Accuracy") + ggtitle("") + ylim(0.68, 0.9) +
    geom_signif(data = pv, aes(x = x, xend = xend, y = y, yend = yend, annotation = values, textsize = 2), 
                stat = "identity", tip_length = 1) + guides(color = F)
g1

# quadratic trend of accuracy
f2 <- glmer(corr ~ I(val3) + I(val3^2) + (1 + I(val3) + I(val3^2)|subj), data = d, family = binomial)
summary(f2)
Anova(f2)
g2 <- plot_model(f2, type = "pred", terms = "val3[all]", pred.type = "fe") +
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) +
    xlab("Dud value") + ylab("Accuracy") + ggtitle("")
g2

# factorial analysis on reaction time
f3 <- lmer(rt ~ fVal3 * chosenItem + (1 + fVal3 + chosenItem|subj), data = filter(d, chosenItem != "dud"))
summary(f3)
Anova(f3)
dl3 <- difflsmeans(f3) # uncorrected linear contrasts 
dl3[17:21,] # p values for correct choices
ap1 <- p.adjust(dl3[17:21,][, 7] , method = "holm") # adjusted p values for correct choices
dl3[68:72,] # p values for incorrect choices
ap2 <- p.adjust(dl3[68:72,][, 7] , method = "holm") # adjusted p values for incorrect choices
pv1 <- data.frame(x = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  y = c(0.75, 0.7, 0.65, 0.6, 0.55, 1.45, 1.5, 1.55, 1.6, 1.65),
                  xend = c(27, 45, 63, 76, 86, 27, 45, 63, 76, 86),
                  yend = c(0.75, 0.7, 0.65, 0.6, 0.55, 1.45, 1.5, 1.55, 1.6, 1.65),
                  group_col = c(rep("target", 5), rep("distractor", 5)),
                  values = (as.character(c(round(ap1, 3), round(ap2, 3)))))
g3 <- plot_model(f3, type = "pred", terms = c("fVal3", "chosenItem"), pred.type = "fe") +
    xlab("Dud value") +  ylab("Response time") + ggtitle("") + ylim(0.55, 1.7) + 
    scale_x_continuous(breaks = c(0, 27, 45, 63, 76, 86)) +
    geom_signif(data = pv1, aes(x = x, xend = xend, y = y, yend = yend, annotation = values, color = group_col, textsize = 2), 
                stat = "identity", tip_length = 1)
g3

# factorial analysis on confidence
f4 <- lmer(conf ~ fVal3 * chosenItem + (1 + fVal3 + chosenItem|subj), data = filter(d, chosenItem != "dud"))
summary(f4)
Anova(f4)
difflsmeans(f4) # uncorrected linear contrasts
dl4 <- difflsmeans(f4) # uncorrected linear contrasts 
dl4[17:21,] # p values for correct choices
ap1 <- p.adjust(dl4[17:21,][, 7] , method = "holm") # adjusted p values for correct choices
dl4[68:72,] # p values for incorrect choices
ap2 <- p.adjust(dl4[68:72,][, 7] , method = "holm") # adjusted p values for incorrect choices
pv2 <- data.frame(x = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  y = c(3.5, 3.6, 3.7, 3.8, 3.9, 1.9, 1.8, 1.7, 1.6, 1.5),
                  xend = c(27, 45, 63, 76, 86, 27, 45, 63, 76, 86),
                  yend = c(3.5, 3.6, 3.7, 3.8, 3.9, 1.9, 1.8, 1.7, 1.6, 1.5),
                  group_col = c(rep("target", 5), rep("distractor", 5)),
                  values = (as.character(c(round(ap1, 3), round(ap2, 3)))))
g4 <- plot_model(f4, type = "pred", terms = c("fVal3", "chosenItem"), pred.type = "fe") +
    xlab("Dud value") + ylab("Confidence") + ggtitle("") + ylim(1.5, 4) + 
    geom_signif(data = pv2, aes(x = x, xend = xend, y = y, yend = yend, annotation = values, color = group_col, textsize = 2),  
                stat = "identity")
g4


#'# save figures
figure_S1 <- cowplot::plot_grid(p1 + theme(legend.text = element_text(size = 6), 
                                           legend.background = element_rect(fill = NA, colour = NA),
                                           legend.key = element_rect(fill = NA),
                                           plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                                p2 + theme(plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                                p3 + theme(axis.text.x = element_text(angle = 90),
                                           plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                                p4 + theme(axis.text.x = element_text(angle = 90),
                                           plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")), 
                                p5 + theme(axis.text.x = element_text(angle = 90),
                                           plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                                p6 + theme(axis.text.x = element_text(angle = 90),
                                           plot.margin = unit(c(0.35, 0.35, 0.35, 0.35), "cm")),
                                labels = c("a", "b", "c", "d", "e", "f"), label_size = 10, nrow = 3)
figure_S1
sjPlot::save_plot("figure_S1.jpg", figure_S1, dpi = 600, height = 13, width = 14)