#' @title Burn-in only
#' @description The burn-in component of the agent-based model of cultural transmission.
#'
#' @param priors A dataframe of priors with the follow columns: ini_pop, ini_syls, innov, dem, p_att, a, v. Each row corresponds to an iteration of the model.
#' @param surv Average overwinter survival.
#' @param rep_m Mean repertoire size across observed dataset.
#' @param rep_sd Standard deviation of repertoire size across observed dataset.
#' @param males_only Whether only males are included in simulation (TRUE/FALSE).
#' @param burn Length of burn-in phase.
#' @param n_iter Number of iterations.
#' @param n_cores Number of cores for parallelization.
#'
#' @return Returns the new repertoire of the learner.
#' @export
burn_only <- function(priors, surv = 0.5, rep_m, rep_sd, males_only = TRUE, burn = 100, n_iter = 100, n_cores = 4){
  #store all parameter values
  priors <- data.table::as.data.table(priors)
  
  #set max range of geographic index (arbitrary; higher numbers just give higher resolution)
  geo <- 100
  
  #generate empty burn-in matrix
  burn_stats <- matrix(data = NA, nrow = n_iter, ncol = burn)
  
  i <- 1
  
  #run simulations
  for(i in 1:n_iter){
    #draw parameter values from priors
    ini_pop <- priors$ini_pop[i] #initial population size (N_b)
    ini_syls <- priors$ini_syls[i] #number of syllables (N_s)
    innov <- priors$innov[i] #innovation rate (mu)
    dem <- priors$dem[i] #number of demonstrators (D)
    p_att <- priors$p_att[i] #proportion of syllables that are attractive (p_att)
    a <- priors$a[i] #level of conformity bias (a)
    v <- priors$v[i] #variation in demonstrator attractiveness (v)
    
    #generate vector of pop sizes
    if(males_only){
      pop_size <- round(ini_pop/2)
    }
    if(!males_only){
      pop_size <- ini_pop
    }
    
    #generate vector of syllable attractiveness
    syl_counter <- rep.int(1, ini_syls)
    syl_counter[sample.int(ini_syls, ini_syls*(1-p_att))] <- 0.05
    
    #generate vector of mortality
    mort <- pop_size*(1-surv)
    
    #initialize agents
    agents <- data.table::data.table(t_x = exp(Rfast::Rnorm(pop_size, 0, v)), loc = sample.int(geo, pop_size, replace = TRUE), reps = replicate(pop_size, sample.int(ini_syls, truncnorm::rtruncnorm(1, 1, ini_syls, rep_m, rep_sd)))) #truncated normal dist with lower boundary of 1 and upper boundary of total number of syls to prevent negative or zero rep sizes
    
    #simulate mortality
    agents <- agents[-sample.int(nrow(agents), round(pop_size/2)),]
    
    #generate empty burn diversity list
    burn_div <- integer(burn)
    
    #burn-in phase
    for(j in 1:burn){
      #calculate simpson diversity
      burn_div[j] <- vegan::diversity(table(unlist(agents$reps)), "simpson")
      
      #store number of new agents
      num_new_agents <- pop_size-nrow(agents)
      
      #learn
      locs <- sample.int(geo, num_new_agents, replace = TRUE)
      rep_sizes <- round(truncnorm::rtruncnorm(num_new_agents, 1, Inf, rep_m, rep_sd)) #THIS IS THE CULPRIT
      new_reps <- parallel::mclapply(1:num_new_agents, function(x){learn(rep_sizes[x], agents$reps, agents$loc, locs[x], agents$t_x, syl_counter, dem, a, innov, p_att)}, mc.cores = n_cores)
      
      #add new agents to population
      agents <- data.table::rbindlist(list(agents, data.table::data.table(t_x = exp(Rfast::Rnorm(num_new_agents, 0, v)), loc = locs, reps = new_reps)))
      
      #simulate mortality
      agents <- agents[-sample.int(nrow(agents), round(pop_size/2)),]
    }
    
    #attach burn diversity list to matrix
    burn_stats[i,] <- burn_div
    
    gc()
    
    #print iteration number
    sink("status.txt")
    print(paste0("Iteration ", i, ": ", (i/n_iter)*100, "%"))
    sink()
  }
  
  return(burn_stats)
}
