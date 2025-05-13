library(tidyverse)
library(lme4)
library(lmerTest)
library(mlogit)

read.csv("exp1_data_eyetracking.csv", header = TRUE) %>% 
  filter(event == "fixation" & fixItem != "other") %>% 
  select(subj, trial, condition, dur, fixItem, countFix_aoi, durFirstFix_fromTrialBegin, fromTrialBegin) %>% 
  mutate(eventEnd_fromTrialBegin = fromTrialBegin + dur) -> df

fix_stats <- data.frame()
prop_fix <- data.frame()
dur_startFix <- list()
fixation_interval <- list()


for (s in unique(df$subj)) {
  d <- filter(df, subj == s)
  i <- which(unique(df$subj) == s)
  
  group_by(d, condition, countFix_aoi) %>% 
    summarise(N = n(), n_target = sum(fixItem == "target"), 
              n_distractor = sum(fixItem == "distractor"), 
              n_dud = sum(fixItem == "dud"), 
              p_target = sum(fixItem == "target") / N, 
              p_distractor = sum(fixItem == "distractor") / N, 
              p_dud = sum(fixItem == "dud") / N) %>%
    as.data.frame() -> prop_fix_indiv
  
  prop_fix <- rbind(prop_fix, mutate(prop_fix_indiv, subj = s))
  
  group_by(d, condition, countFix_aoi, fixItem) %>% 
    summarise(mean_duration = mean(dur)) %>% 
    as.data.frame() %>% 
    pivot_wider(names_from = fixItem, values_from = mean_duration) %>% 
    as.data.frame() %>% 
    rename(dur_distractor = distractor, dur_target = target, dur_dud = dud) %>% 
    select(dur_target, dur_distractor, dur_dud) -> dur_fix_indiv
  dur_fix_indiv[is.na(dur_fix_indiv)] <- 0
  
  bind_cols(prop_fix_indiv, dur_fix_indiv) %>% 
    mutate(subj = s) %>% 
    relocate(subj, .before = condition) -> fix_indiv
  
  fix_stats <- bind_rows(fix_stats, fix_indiv)
  
  group_by(d, trial) %>% 
    summarise(dur_startFix = unique(durFirstFix_fromTrialBegin)) -> dur_startFix[[i]]
  
  fixation_interval_indiv <- list()
  for (n in unique(d$trial)) {
    d_trial <- filter(d, trial == n)
    if (nrow(d_trial != 0)) {
      fixation_interval_indiv[[n]] <- d_trial$fromTrialBegin[2:nrow(d_trial)] - d_trial$eventEnd_fromTrialBegin[1:(nrow(d_trial) - 1)]
    }
  }
  fixation_interval[[i]] <- fixation_interval_indiv
}
write.csv(fix_stats, "exp1_fixation_statistics.csv", row.names = FALSE)

pivot_longer(prop_fix, cols = c(p_target, p_distractor, p_dud), names_to = "fixItem", values_to = "P") %>% 
  filter(countFix_aoi <= 5) -> d_plot
ggplot(d_plot, aes(x = condition, y = P, color = fixItem)) + geom_line(linewidth = 1.2) + 
  facet_wrap(countFix_aoi ~ subj, ncol = 10) + 
  geom_text(aes(group = countFix_aoi, label = N), vjust = -2, data = filter(d_plot, fixItem == "p_target"))