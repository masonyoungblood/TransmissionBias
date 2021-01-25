#' @title Agent-based model
#' @description An agent-based model of cultural transmission that incorporates content bias, frequency bias, and demonstrator bias, as well as changes in population size.
#'
#' @param priors A dataframe of priors with the follow columns: ini_pop, ini_syls, innov, dem, p_att, a, v. Each row corresponds to an iteration of the model.
#' @param surv Average overwinter survival.
#' @param pop_trends A vector of populations trends.
#' @param obs_years Years with observed data.
#' @param obs_n Sample size in years with observed data.
#' @param rep_m Mean repertoire size across observed dataset.
#' @param rep_sd Standard deviation of repertoire size across observed dataset.
#' @param males_only Whether only males are included in simulation (TRUE/FALSE).
#' @param burn Length of burn-in phase.
#' @param n_iter Number of iterations.
#' @param n_cores Number of cores used in the learning step. For large population and repertoire sizes it may be faster to run the ABM function serially (n_cores = 1) in a parallelized for loop.
#'
#' @return Returns a list of two matrices: (1) the Simpson's diversity from each burn-in year (column) and iteration (row); (2) the summary statistics corresponding to each observed year (subsequent columns) and iteration (row). The returned summary statistics are the proportion of syllables that only appear once, the proportion of the most common syllable type, the number of syllable types, Simpson's diversity index, Shannon's diversity index, Pielou's evenness index, and the exponent of the fitted power-law function to the progeny distribution.
#' @export
ABM <- function(priors, surv = 0.5, pop_trends, obs_years, obs_n, rep_m, rep_sd, males_only = TRUE, burn = 100, n_iter = 100, n_cores = 4){
  #store all parameter values
  priors <- data.table::as.data.table(priors)
  
  #set max range of geographic index (arbitrary; higher numbers just give higher resolution)
  geo <- 100
  
  #generate empty summary statistic matrix
  sum_stats <- matrix(data = NA, nrow = n_iter, ncol = length(obs_years)*7) #length is number of years times seven because of the seven summary statistics
  
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
      pop_sizes <- round((ini_pop*pop_trends)/2)
    }
    if(!males_only){
      pop_sizes <- round(ini_pop*pop_trends)
    }
    
    #set number of years
    n_years <- length(pop_sizes)
    
    #generate vector of syllable attractiveness
    syl_counter <- rep.int(1, ini_syls)
    syl_counter[sample.int(ini_syls, ini_syls*(1-p_att))] <- 0.05
    
    #generate vector of mortality
    mort <- pop_sizes*(1-surv)
    
    #identify years where population crashes over 50% and replace mortality
    crash_years <- which(pop_sizes[2:n_years]/pop_sizes[1:n_years-1] <= 0.5)
    mort[crash_years] <- pop_sizes[crash_years]-pop_sizes[crash_years+1]*surv
    
    #initialize agents
    agents <- data.table::data.table(t_x = exp(Rfast::Rnorm(pop_sizes[1], 0, v)), loc = sample.int(geo, pop_sizes[1], replace = TRUE), reps = replicate(pop_sizes[1], sample.int(ini_syls, truncnorm::rtruncnorm(1, 1, ini_syls, rep_m, rep_sd)))) #truncated normal dist with lower boundary of 1 and upper boundary of total number of syls to prevent negative or zero rep sizes
    
    #simulate mortality
    agents <- agents[-sample.int(nrow(agents), round(pop_sizes[1]/2)),]

    #generate empty burn diversity list
    burn_div <- integer(burn)

    #burn-in phase
    for(j in 1:burn){
      #calculate simpson diversity
      burn_div[j] <- vegan::diversity(table(unlist(agents$reps)), "simpson")

      #store number of new agents
      num_new_agents <- pop_sizes[1]-nrow(agents)

      #learn
      locs <- sample.int(geo, num_new_agents, replace = TRUE)
      rep_sizes <- round(truncnorm::rtruncnorm(num_new_agents, 1, Inf, rep_m, rep_sd)) #THIS IS THE CULPRIT
      new_reps <- parallel::mclapply(1:num_new_agents, function(x){learn(rep_sizes[x], agents$reps, agents$loc, locs[x], agents$t_x, syl_counter, dem, a, innov, p_att)}, mc.cores = n_cores)

      #add new agents to population
      agents <- data.table::rbindlist(list(agents, data.table::data.table(t_x = exp(Rfast::Rnorm(num_new_agents, 0, v)), loc = locs, reps = new_reps)))

      #simulate mortality
      agents <- agents[-sample.int(nrow(agents), round(pop_sizes[1]/2)),]
    }

    #attach burn diversity list to matrix
    burn_stats[i,] <- burn_div

    #create empty sum_stats_temp vector for this set of parameter values
    sum_stats_temp <- c()

    #run iterations
    for(j in 1:n_years){
      #calculate sum_stats
      if(j %in% obs_years){
        raw_data <- unlist(agents$reps[sample.int(length(agents$reps), obs_n[which(obs_years == j)], replace = FALSE)]) #from random subset of obs_n individuals
        raw_data_table <- table(raw_data)
        sum_stats_temp <- c(sum_stats_temp,
                            length(which(as.numeric(raw_data_table) == 1))/sum(raw_data_table), #proportion of syllable types that only appear once (among all syllables, independent of type)
                            max(raw_data_table)/sum(raw_data_table), #proportion of the most common syllable type (among all syllables, independent of type)
                            length(unique(raw_data)), #number of syllable types
                            vegan::diversity(as.numeric(raw_data_table), "simpson"), #simpson's diversity index
                            vegan::diversity(as.numeric(raw_data_table), "shannon"), #shannon's diversity index
                            vegan::diversity(as.numeric(raw_data_table), "shannon")/log(vegan::specnumber(as.numeric(raw_data_table))+0.00001), #pielou's evenness index
                            igraph::power.law.fit(raw_data)$alpha) #exponent of the fitted power-law function to the progeny distribution
      }

      #store number of new agents
      num_new_agents <- pop_sizes[j]-nrow(agents)

      #learn
      locs <- sample.int(geo, num_new_agents, replace = TRUE)
      rep_sizes <- round(truncnorm::rtruncnorm(num_new_agents, 1, Inf, rep_m, rep_sd))
      new_reps <- parallel::mclapply(1:num_new_agents, function(x){learn(rep_sizes[x], agents$reps, agents$loc, locs[x], agents$t_x, syl_counter, dem, a, innov, p_att)}, mc.cores = n_cores)

      #add new agents to population
      agents <- data.table::rbindlist(list(agents, data.table::data.table(t_x = exp(Rfast::Rnorm(num_new_agents, 0, v)), loc = locs, reps = new_reps)))

      #simulate mortality
      agents <- agents[-sample.int(pop_sizes[j], mort[j]),]
    }

    #add summary statistics to matrix
    sum_stats[i,] <- sum_stats_temp

    gc()

    #print iteration number
    sink("status.txt")
    print(paste0("Iteration ", i, ": ", (i/n_iter)*100, "%"))
    sink()
  }

  return(list(burn_stats, sum_stats))
}
