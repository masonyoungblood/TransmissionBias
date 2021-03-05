# TransmissionBias
*Agent-Based Model of Biased Cultural Transmission*

This package includes an agent-based model (ABM) of cultural transmission that simulates the three main forms of transmission bias: content bias, frequency bias, and demonstrator bias. While it was designed to simulate the cultural transmission of birdsong, it is theoretically generalizable to any cultural system where individuals learn a repertoire (or a set) of cultural traits. The ABM is based on [Lachlan et al. (2018)](https://www.nature.com/articles/s41467-018-04728-1), but has been expanded to allow for dynamic population size. The back-end of learning step of the ABM was written in C++ to minimize computation time.

To install the package run the following in R:

```
install.packages("devtools")
devtools::install_github("masonyoungblood/TransmissionBias”)
```

A basic example of the core ABM function is below. In this example, the priors allow for a relatively small population size (100-200), and the simulation is run for three years with no burn-in phase. You can access the detailed documentation of the ABM function by entering ?ABM in R:

```
priors <- data.frame(ini_pop = round(runif(n_iter, 100, 200)), #uniform: 100-200
  ini_syls = round(KScorrect::rlunif(n_iter, 5, 20)), #log uniform: 5-20
  innov = KScorrect::rlunif(n_iter, 0.001, 0.1), #log uniform: 0.001-0.1
  dem = round(KScorrect::rlunif(n_iter, 2, 10)), #log uniform: 2-10
  p_att = runif(n_iter, 0.01, 1), #uniform: 0.01-1
  a = runif(n_iter, 0.25, 3), #uniform: 0.25-3
  v = KScorrect::rlunif(n_iter, 0.01, 3)) #log uniform: 0.01-3
  
simulations <- ABM(priors = priors, pop_trends = pop_trends,
  obs_years = c(1, 2, 3), obs_n = c(20, 20, 20),
  rep_m = 10, rep_sd = 2, burn = 0, n_iter = 1, n_cores = 1)
```

For more details about the methods in this package, check out the corresponding manuscript:

Youngblood, M., & Lahti, D. (2020). Content bias in the cultural evolution of house finch song. *In preparation*.
