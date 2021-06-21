library(foreach)
library(TransmissionBias)
library(Rcpp)
library(abcrf)

#load the christmas bird count population trends
load("cbc_trend.RData")

#add final year since trends only go through 2018
cbc_trend <- rbind(cbc_trend, c(2019, 1))
rownames(cbc_trend) <- NULL

#set up for parallelization
cl <- parallel::makeForkCluster(50)
doParallel::registerDoParallel(cl)

#run the following four times, each time adjusting the values of i so that you end up with 200,000 simulations
foreach(i = 1:50) %dopar% {
  set.seed(i)
  
  #set number of iterations (within the four manual loops of 50)
  n_iter <- 1000
  
  #set priors
  priors <- data.frame(ini_pop = round(runif(n_iter, 2000, 10000)), #uniform: 2,000-10,000
                       ini_syls = round(KScorrect::rlunif(n_iter, 596, 800)), #log uniform: 596-800
                       innov = KScorrect::rlunif(n_iter, 0.001, 0.3), #log uniform: 0.001-0.3
                       dem = round(KScorrect::rlunif(n_iter, 2, 10)), #log uniform: 2-10
                       p_att = runif(n_iter, 0.01, 1), #uniform: 0.01-1
                       a = runif(n_iter, 0.25, 3), #uniform: 0.25-3
                       v = KScorrect::rlunif(n_iter, 0.01, 6)) #log uniform: 0.01-6
  
  #run simulations
  simulations <- ABM(priors = priors, pop_trends = cbc_trend$trend, obs_years = c(6, 43, 50), obs_n = c(36, 75, 29), rep_m = 61, rep_sd = 18, burn = 100, n_iter = n_iter, n_cores = 1)
  
  #save output
  save(priors, file = paste0("simulation_output/priors_", i, ".RData"))
  save(simulations, file = paste0("simulation_output/simulations_", i, ".RData"))
}

#end parallelization
parallel::stopCluster(cl)

#load observed summary statistics
load("obs_sum_stats.RData")
obs_sum_stats <- as.data.frame(t(obs_sum_stats))

#load simulation files (from above)
files <- c(1:200)
for(i in files){
  load(paste0("simulation_output/simulations_", i, ".RData"))
  load(paste0("simulation_output/priors_", i, ".RData"))
  if(i == 1){
    all_simulations <- simulations[[2]]
    all_priors <- priors
  }
  if(i > 1){
    all_simulations <- rbind(all_simulations, simulations[[2]])
    all_priors <- rbind(all_priors, priors)
  }
  rm(list = c("simulations", "priors"))
}

#create predictions list
predictions <- list()

#iterate through the following (i <- 1... i <- i + 1) for each parameter - manual looping is required to save density plots
i <- 1

set.seed(i)

#create data frame for abcrf
abcrf_data <- data.frame(param = all_priors[,i], sum_stats = all_simulations)
colnames(obs_sum_stats) <- colnames(abcrf_data)[-1]

#run abcrf
reg_abcrf <- regAbcrf(formula = param ~ ., data = abcrf_data, ntree = 1000, mtry = 7, min.node.size = 5, sampsize = nrow(all_simulations), max.depth = 5, paral = TRUE, ncores = 11)

#create density plot
densityPlot(object = reg_abcrf, obs = obs_sum_stats, training = abcrf_data, paral = TRUE, ncores = 11, ylab = "Density", xlab = "Parameter Value", main = "Title")
density_plot <- recordPlot()

#store everything in predictions file
predictions[[i]] <- list(OOB_MSE = reg_abcrf$model.rf$prediction.error, OOB_NMAE = reg_abcrf$model.rf$NMAE, prediction = predict(object = reg_abcrf, obs = obs_sum_stats, training = abcrf_data, post.err.med = TRUE, paral = TRUE, ncores = 8), var_importance = sort(reg_abcrf$model.rf$variable.importance, decreasing = TRUE), plot = NULL)
predictions[[i]]$plot <- density_plot

#remove objects
rm(list = c("abcrf_data", "reg_abcrf", "density_plot"))

i <- i + 1

#add names to predictions file
names(predictions) <- colnames(all_priors)

#save predictions file
save(predictions, file = "predictions.RData")
