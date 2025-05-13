library(tidyverse)
library(data.table)

################################################################################
# files <- list.files(recursive = F, pattern = ".csv")
# d1 <- fread(files[2], header = T)
# d1 <- d1[, 5:38]
# d1 <- d1[, -9]
# d2 <- fread(files[4], header = T)
# d2 <- d2[, 1:33]
# d <- rbind(d1, d2)
# d <- na.omit(d)
# d <- d[1:576, ]
# write.csv(d, files[2])
# 
# d1 <- fread(files[3], header = T)
# d1 <- d1[, 5:38]
# d1 <- d1[, -9]
# d1 <- d1[-nrow(d1), ]
# d2 <- fread(files[5], header = T)
# d2 <- d2[, 1:33]
# d <- rbind(d1, d2)
# d <- na.omit(d)
# write.csv(d, files[3])
# 
# d1 <- fread(files[23], header = T)
# d1 <- d1[, 1:33]
# d1 <- d1[-nrow(d1), ]
# d2 <- fread(files[24], header = T)
# d2 <- d2[, 1:33]
# d <- rbind(d1, d2)
# 
# write.csv(d, files[23])
################################################################################
files <- list.files(recursive = T, pattern = ".csv")
df <- data.frame()
for (i in 1:length(files)) {
  d <- fread(files[i], header = T)
  d <- select(d, participant, Condition, Nalt, Val1, Val2, Val3, UVal, LVal, RVal, ChoiceRT, ChosenITM, CorrectITM, Conf, ConfRT)
  d <- na.omit(d)
  d <- as.data.frame(d)
  df <- rbind(df, d)
}

write.csv(df, "exp1_data_behavior.csv", row.names = FALSE)
