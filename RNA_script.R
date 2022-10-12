if(!require(tidyverse))
  install.packages("tidyverse")
if(!require(ggplot2))
  install.packages("ggplot2")

library(readxl)
library(tidyverse)
library(ggplot2)

rawData <- read_excel("Metformin_RNA.xlsx", col_names = T)

# Even after simplifying the column names, they are hard to understand
# So instead of L,R,N, etc., I will call them 1,2,3 respectively
colnames(rawData) <- 
  c("Gene", "AET-1", "AET-2", "AET-3", "MET-1", "MET-2", "MET-3"
    , "SED-1", "SED-2", "SED-3")

# Since in this experiment N=3, I would like to remove any rows with 0s
rawData[rawData == 0] <- NA
validData <- rawData %>% na.omit()

# Replace all hyphens (-) with underscores (_) in the Gene column
# to avoid ANOVA bugging itself out
validData$Gene <- gsub("-", "_", validData$Gene)

# Change the dataframe from wide form to long form
tidyNormData <- validData %>% 
  pivot_longer(c(2:dim(validData)[2]), names_to = "keys", values_to = "intensity") %>% 
  extract(keys, into = "sample", regex = "([[:alpha:]]{0,3})", remove = F) %>% 
  extract(keys, into = "group", regex = "([[:digit:]]{1,})", remove = T)

#See if data is normally distributed when mutated with log2
tidyNormData <- tidyNormData %>%
  mutate(intensity = log2(intensity))

#Plot histogram for all intensity values.
tidyNormData %>%
  ggplot(aes(x = intensity)) + geom_histogram() + 
  theme_bw()

#Looks pretty normal, just slightly right-skewed. 
#But just to check: Make it a little more normal
tidyRNANorm <- tidyNormData %>%
  group_by(sample) %>%
  mutate(normIntensity = scale(intensity, center = T, scale = T))

#Widen data again from normalized state
statsRNANorm <- tidyRNANorm[,-c(4)] %>% 
  unite("sm", sample:group, sep = "", remove = T) %>% 
  # sep = 3 means the first 3 letters. This value does not matter, can be anything
  separate(sm, into = "group", sep = 3, remove = F)

genRef <- statsRNANorm[,c(1)] %>% distinct()

# This is the most important step for ANOVA - grouping
groupedRNAStats <- statsRNANorm %>% 
  pivot_wider(names_from = "Gene", values_from = "normIntensity")

# Convert into dataframe
groupedRNAStats <- as.data.frame(groupedRNAStats)
row.names(groupedRNAStats) <- groupedRNAStats$sm

groupedRNAStats <- groupedRNAStats[,-1] %>%
  mutate(group = factor(group, levels = c("AET", "MET", "SED")))

#Start anova testing
formula <- 
  lapply(colnames(groupedRNAStats)[2:ncol(groupedRNAStats)],
         function(x) as.formula(paste0(x, " ~group")))
