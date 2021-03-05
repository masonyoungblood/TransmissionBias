library(fastcluster)
library(dynamicTreeCut)
library(WGCNA)
library(lme4)
library(rstanarm)

load("data.RData")
load("clustering.RData")
load("hybrid_cut.RData")

#add syllable types to data file
data$cluster <- as.numeric(hybrid_cut)

#subset data to only include individuals with at least 8 recorded songs
num_songs_per_individual <- sapply(1:length(unique(data$individual)), function(x){length(unique(data$song[which(data$individual == unique(data$individual)[x])]))})
data <- data[which(data$individual %in% unique(data$individual)[which(num_songs_per_individual >= 8)]), ]

#calculate concavity
inds <- 1:nrow(data)
inds <- inds[-c(25001)] #remove 25,001, which is too short to analyze - it receives the correct value of 0 (visual inspection) from the code below
concavity <- rep(0, nrow(data))
for(i in inds){
  temp <- pspline::sm.spline(data$meanfreq[[i]], spar = 5)
  temp <- diff(c(temp$ysmth)) #extract first derivative (slopes)
  concav <- 0
  for(j in 2:length(temp)){
    if((temp[[j]] < 0 & temp[[j-1]] > 0) | (temp[[j]] > 0 & temp[[j-1]] < 0)){
      concavity[i] <- concavity[i] + 1
    }
  }
  concavity[i] <- concavity[i]/length(temp)
}

#calculate excursion
excursion <- rep(0, nrow(data))
for(i in inds){
  temp <- data$meanfreq[[i]]
  excursion[i] <- sum(sapply(2:length(temp), function(x){abs(temp[x]-temp[x-1])}))/length(temp)
}

#add to data table
data$concavity <- concavity
data$excursion <- excursion

#calculate other measures from mean frequency trace
data$average <- sapply(1:nrow(data), function(x){mean(data$meanfreq[[x]])})
data$min <- sapply(1:nrow(data), function(x){min(data$meanfreq[[x]])})
data$max <- sapply(1:nrow(data), function(x){max(data$meanfreq[[x]])})
data$bandwidth <- sapply(1:nrow(data), function(x){max(data$meanfreq[[x]])-min(data$meanfreq[[x]])}) 
data$duration <- sapply(1:nrow(data), function(x){length(data$meanfreq[[x]])}) 

#store as separate data tables for each year
data_1975 <- data[which(substring(data$individual, 1, 4) == "1975"),]
data_2012 <- data[which(substring(data$individual, 1, 4) == "2012"),]
data_2019 <- data[which(substring(data$individual, 1, 4) == "2019"),]

#calculate observed summary statistics
obs_sum_stats <- c(length(which(as.numeric(table(data_1975$cluster)) == 1))/sum(table(data_1975$cluster)), #proportion of syllable types that only appear once (among all syllables, independent of type)
                   max(table(data_1975$cluster))/sum(table(data_1975$cluster)), #proportion of the most common syllable type (among all syllables, independent of type)
                   length(unique(data_1975$cluster)), #number of syllable types
                   vegan::diversity(as.numeric(table(data_1975$cluster)), "simpson"), #simpson's diversity index
                   vegan::diversity(as.numeric(table(data_1975$cluster)), "shannon"), #shannon's diversity index
                   vegan::diversity(as.numeric(table(data_1975$cluster)), "shannon")/log(vegan::specnumber(as.numeric(table(data_1975$cluster)))+0.00001), #pielou's evenness index
                   igraph::power.law.fit(data_1975$cluster)$alpha, #exponent of the fitted power-law function to the progeny distribution
                   length(which(as.numeric(table(data_2012$cluster)) == 1))/sum(table(data_2012$cluster)), #proportion of syllable types that only appear once (among all syllables, independent of type)
                   max(table(data_2012$cluster))/sum(table(data_2012$cluster)), #proportion of the most common syllable type (among all syllables, independent of type)
                   length(unique(data_2012$cluster)), #number of syllable types
                   vegan::diversity(as.numeric(table(data_2012$cluster)), "simpson"), #simpson's diversity index
                   vegan::diversity(as.numeric(table(data_2012$cluster)), "shannon"), #shannon's diversity index
                   vegan::diversity(as.numeric(table(data_2012$cluster)), "shannon")/log(vegan::specnumber(as.numeric(table(data_2012$cluster)))+0.00001), #pielou's evenness index
                   igraph::power.law.fit(data_2012$cluster)$alpha, #exponent of the fitted power-law function to the progeny distribution
                   length(which(as.numeric(table(data_2019$cluster)) == 1))/sum(table(data_2019$cluster)), #proportion of syllable types that only appear once (among all syllables, independent of type)
                   max(table(data_2019$cluster))/sum(table(data_2019$cluster)), #proportion of the most common syllable type (among all syllables, independent of type)
                   length(unique(data_2019$cluster)), #number of syllable types
                   vegan::diversity(as.numeric(table(data_2019$cluster)), "simpson"), #simpson's diversity index
                   vegan::diversity(as.numeric(table(data_2019$cluster)), "shannon"), #shannon's diversity index
                   vegan::diversity(as.numeric(table(data_2019$cluster)), "shannon")/log(vegan::specnumber(as.numeric(table(data_2019$cluster)))+0.00001), #pielou's evenness index
                   igraph::power.law.fit(data_2019$cluster)$alpha) #exponent of the fitted power-law function to the progeny distribution

#save observed summary statistics
save(obs_sum_stats, file = "obs_sum_stats.RData")

#create data frame of syllables that persisted from 1975 to 2019
syl_types_1975_2019 <- data.frame(type = unique(data_1975$cluster),
                                  concavity = rep(NA, length(unique(data_1975$cluster))),
                                  excursion = rep(NA, length(unique(data_1975$cluster))),
                                  average = rep(NA, length(unique(data_1975$cluster))),
                                  min = rep(NA, length(unique(data_1975$cluster))),
                                  max = rep(NA, length(unique(data_1975$cluster))),
                                  bandwidth = rep(NA, length(unique(data_1975$cluster))),
                                  duration = rep(NA, length(unique(data_1975$cluster))))

#get average measures for each syllable type
for(i in 1:nrow(syl_types_1975_2019)){
  syl_types_1975_2019$concavity[i] <- mean(data$concavity[which(data$cluster == syl_types_1975_2019$type[i])])
  syl_types_1975_2019$excursion[i] <- mean(data$excursion[which(data$cluster == syl_types_1975_2019$type[i])])
  syl_types_1975_2019$average[i] <- mean(data$average[which(data$cluster == syl_types_1975_2019$type[i])])
  syl_types_1975_2019$min[i] <- mean(data$min[which(data$cluster == syl_types_1975_2019$type[i])])
  syl_types_1975_2019$max[i] <- mean(data$max[which(data$cluster == syl_types_1975_2019$type[i])])
  syl_types_1975_2019$bandwidth[i] <- mean(data$bandwidth[which(data$cluster == syl_types_1975_2019$type[i])])
  syl_types_1975_2019$duration[i] <- mean(data$duration[which(data$cluster == syl_types_1975_2019$type[i])])
}

#add whether or not the syllable also appears in 2019
syl_types_1975_2019$persisted <- sapply(1:nrow(syl_types_1975_2019), function(x){syl_types_1975_2019$type[x] %in% data_2019$cluster})

#VIF testing for multicollinearity
usdm::vif(syl_types_1975_2019[,c(2, 3, 4, 5, 6, 7, 8)])
usdm::vif(syl_types_1975_2019[,c(2, 3, 5, 8)])

#bayesian logistic regression
glm_results <- stan_glm(persisted ~ scale(concavity) + scale(excursion) + scale(min) + scale(duration), 
                        data = syl_types_1975_2019, family = binomial, 
                        prior = student_t(location = 0, scale = 2.5),
                        prior_intercept = student_t(location = 0, scale = 2.5))
