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
library(emmeans)

theme_set(theme_publication()) 


#'# data loading
d <- fread("exp2_data_3afc_processed.csv", sep = ",", header = T)
d$rt <- as.numeric(d$rt)
d <- na.omit(d) # no trials omitted here
d$chosenItem <- factor(d$chosenItem, levels = c("target", "distractor", "dud"))
d$fVal3 <- as.factor(d$val3)
d <- mutate(d, confidence = ifelse(conf > 2, 1, 0))
# d <- subset(d, d$subj != "subj21")


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
f4 <- glmer(confidence ~ fVal3 * chosenItem + (1 + fVal3 + chosenItem|subj), 
            data = filter(d, chosenItem != "dud"), family = binomial)
summary(f4)
Anova(f4)
em1 <- emmeans(f4, ~ fVal3 * chosenItem)
contrast(em1, method = "pairwise")
em2 <- emmeans(f4, ~ fVal3)
contrast(em2, method = "pairwise")