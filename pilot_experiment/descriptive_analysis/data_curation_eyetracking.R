library(rhdf5)
library(tidyverse)

files <- list.files(path = "rawdata_eyetracking", pattern = ".hdf5")

fixateItem <- function(x, y){
  dist <- 5
  pos_up <- c(0, dist * 2 - 1.5)
  pos_left <- c(-sqrt(3) * dist, -dist - 1.5)
  pos_right <- c(sqrt(3) + dist, -dist - 1.5)
  
  dist_up <- sqrt((pos_up[1] - x) ^ 2 + (pos_up[2] - y) ^ 2)
  dist_left <- sqrt((pos_left[1] - x) ^ 2 + (pos_left[2] - y) ^ 2)
  dist_right <- sqrt((pos_right[1] - x) ^ 2 + (pos_right[2] - y) ^ 2)
  
  if (dist_up < dist){
    item <- "up"
  } else if (dist_left < dist){
    item <- "left"
  } else if (dist_right < dist){
    item <- "right"
  } else {
    item <- "noFix"
  }
  
  return (item)
}

df_indiv <- data.frame()
nTrial <- 576
id_list <- c("sub01", "sub02", "sub03", "sub04", "sub05", "sub06", "sub07", "sub08", "sub09", "sub10")
id_list_itr <- 1
N <- 0

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
    select(time, left_gaze_x, left_gaze_y, left_pupil_measure1, right_gaze_x, right_gaze_y, right_pupil_measure1) %>% 
    mutate(gaze_x = (left_gaze_x + right_gaze_x) / 2, gaze_y = (left_gaze_y + right_gaze_y) / 2) -> gaze_indiv
  
  gaze_indiv <- gaze_indiv[!is.nan(gaze_indiv$gaze_x) & !is.nan(gaze_indiv$gaze_y), ]
  df <- data.frame()
  
  for (j in 1:nrow(event_start)){ #loop over trials
    gaze <- subset(gaze_indiv, gaze_indiv$time >= event_start$time[j] & gaze_indiv$time <= event_end$time[j]) #試行ごとのループ parallel OK j
    gaze$item <- mapply(fixateItem, gaze$gaze_x, gaze$gaze_y)
    time_index <- c(1, cumsum(rle(gaze$item)$length))
    item_trial <- rle(gaze$item)$values
    duration_trial <- vector()
    for (k in 1:length(time_index) - 1){
      duration_trial[k] <- gaze$time[time_index[k + 1]] - gaze$time[k]
    }
    if (nrow(gaze) != 0){
      d <- data.frame(duration = duration_trial, item = item_trial, trial = j + N)
      df <- rbind(df, d)
      
    }
  }
  
  N <- N + nrow(event_start)
  
  df_indiv <- rbind(df_indiv, df)
  if (df_indiv$trial[nrow(df_indiv)] >= nTrial * 3){
    df_indiv <- df_indiv[df_indiv$trial <= nTrial * 3, ]
    df_indiv$participant <- id_list[id_list_itr]
    setwd("..")
    setwd("exp1_data_eyetracking_individual")
    write.csv(df_indiv, paste(id_list[id_list_itr], "eyetracking.csv", sep = "_"), row.names = FALSE)
    id_list_itr <- id_list_itr + 1
    df_indiv <- data.frame()
    N <- 0
    setwd("..")
  } else {
    setwd("..")
  }
}
files_indiv <- list.files(path = "exp1_data_eyetracking_individual", pattern = ".csv")
df_all <- data.frame()
setwd("exp1_data_eyetracking_individual")
for (i in 1:length(files_indiv)){
  df_all <- rbind(df_all, read.csv(files_indiv[i], header = TRUE))
}
setwd("..")
write.csv(df_all, "exp1_data_eyetracking.csv", row.names = FALSE)