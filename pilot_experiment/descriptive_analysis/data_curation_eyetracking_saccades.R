library(rhdf5)
library(tidyverse)
library(saccades)
library(purrr)

files <- list.files(path = "rawdata_eyetracking", pattern = ".hdf5")
df_behav <- read.csv("exp1_data_behavior.csv", header = TRUE)
df_indiv <- data.frame()
nTrial <- 576
id_list <- c("sub01", "sub02", "sub03", "sub04", "sub05", "sub06", "sub07", "sub08", "sub09", "sub10")
id_list_itr <- 1
N <- 0

fixateItem <- function(x, y){
  dist <- 5
  rad <- 5
  
  pos_up <- c(0, dist * 2 - 1.5)
  pos_left <- c(-sqrt(3) * dist, -dist - 1.5)
  pos_right <- c(sqrt(3) * dist, -dist - 1.5)
  
  dist_up <- sqrt((pos_up[1] - x) ^ 2 + (pos_up[2] - y) ^ 2)
  dist_left <- sqrt((pos_left[1] - x) ^ 2 + (pos_left[2] - y) ^ 2)
  dist_right <- sqrt((pos_right[1] - x) ^ 2 + (pos_right[2] - y) ^ 2)
  
  if (dist_up < rad){
    item <- "up"
  } else if (dist_left < rad){
    item <- "left"
  } else if (dist_right < rad){
    item <- "right"
  } else {
    item <- "noFix"
  }
  
  return (item)
}

for (i in 1:length(files)){
  setwd("rawdata_eyetracking")
  h5read(files[i], "/data_collection/events/experiment/MessageEvent") %>% 
    filter(category == "choice") %>% 
    select(time, category, text) -> event
  if (nrow(event) %% 144 != 0){
    event <- event[1:(floor(nrow(event) / 144) * 144), ]
  }
  event_start <- filter(event, text == "choice_start")
  event_end <- filter(event, text == "choice_end")
  
  h5read(files[i], "/data_collection/events/eyetracker/BinocularEyeSampleEvent") %>% 
    filter(time >= event$time[1]) %>% 
    mutate(gaze_x = (left_gaze_x + right_gaze_x) / 2, gaze_y = (left_gaze_y + right_gaze_y) / 2) %>% 
    select(time, gaze_x, gaze_y) -> gaze_indiv
  
  gaze_indiv <- gaze_indiv[!is.nan(gaze_indiv$gaze_x) & !is.nan(gaze_indiv$gaze_y), ]
  df <- data.frame()
  
  for (j in 1:nrow(event_start)){ #loop over trials
    gaze <- subset(gaze_indiv, gaze_indiv$time >= event_start$time[j] & gaze_indiv$time <= event_end$time[j])
    if (nrow(gaze) == 0){
      next
    }
    gaze$trial <- j + N
    colnames(gaze) <- c("time", "x", "y", "trial")
    tryCatch(
      {fixation <- detect.fixations(gaze)}
      , error = function(e){fixation <- NULL}
    ) 
    df <- bind_rows(df, fixation)
  }
  
  N <- N + nrow(event_start)
  
  df_indiv <- bind_rows(df_indiv, df)
  if (df_indiv$trial[nrow(df_indiv)] >= nTrial * 3){
    df_indiv <- df_indiv[df_indiv$trial <= nTrial * 3, ]
    df_indiv$participant <- id_list[id_list_itr]
    setwd("..")
    setwd("exp1_data_eyetracking_individual_saccades")
    df_behav_indiv <- filter(df_behav, participant == id_list[id_list_itr])
    for (j in 1:nrow(df_indiv)){
      df_indiv$Condition[j] <- df_behav_indiv$Condition[df_indiv$trial[j]]
      df_indiv$UVal[j] <- df_behav_indiv$UVal[df_indiv$trial[j]]
      df_indiv$LVal[j] <- df_behav_indiv$LVal[df_indiv$trial[j]]
      df_indiv$RVal[j] <- df_behav_indiv$RVal[df_indiv$trial[j]]
      df_indiv$ChosenITM[j] <- df_behav_indiv$ChosenITM[df_indiv$trial[j]]
      df_indiv$CorrectITM[j] <- df_behav_indiv$CorrectITM[df_indiv$trial[j]]
      df_indiv$Conf[j] <- df_behav_indiv$Conf[df_indiv$trial[j]]
    }
    df_indiv$item <- unlist(map2(df_indiv$x, df_indiv$y, fixateItem))
    write.csv(df_indiv, paste(id_list[id_list_itr], "eyetracking.csv", sep = "_"), row.names = FALSE)
    id_list_itr <- id_list_itr + 1
    df_indiv <- data.frame()
    N <- 0
    setwd("..")
  } else {
    setwd("..")
  }
}