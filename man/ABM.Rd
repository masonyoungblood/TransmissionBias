% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABM.R
\name{ABM}
\alias{ABM}
\title{Agent-based model}
\usage{
ABM(
  priors,
  surv = 0.5,
  pop_trends,
  obs_years,
  obs_n,
  rep_m,
  rep_sd,
  males_only = TRUE,
  burn = 100,
  n_iter = 100,
  n_cores = 4
)
}
\arguments{
\item{priors}{A dataframe of priors with the follow columns: ini_pop, ini_syls, innov, dem, p_att, a, v. Each row corresponds to an iteration of the model.}

\item{surv}{Average overwinter survival.}

\item{pop_trends}{A vector of populations trends.}

\item{obs_years}{Years with observed data.}

\item{obs_n}{Sample size in years with observed data.}

\item{rep_m}{Mean repertoire size across observed dataset.}

\item{rep_sd}{Standard deviation of repertoire size across observed dataset.}

\item{males_only}{Whether only males are included in simulation (TRUE/FALSE).}

\item{burn}{Length of burn-in phase.}

\item{n_iter}{Number of iterations.}

\item{n_cores}{Number of cores used in the learning step. For large population and repertoire sizes it may be faster to run the ABM function serially (n_cores = 1) in a parallelized for loop.}
}
\value{
Returns a list of two matrices: (1) the Simpson's diversity from each burn-in year (column) and iteration (row); (2) the summary statistics corresponding to each observed year (subsequent columns) and iteration (row). The returned summary statistics are the proportion of syllables that only appear once, the proportion of the most common syllable type, the number of syllable types, Simpson's diversity index, Shannon's diversity index, Pielou's evenness index, and the exponent of the fitted power-law function to the progeny distribution.
}
\description{
An agent-based model of cultural transmission that incorporates content bias, frequency bias, and demonstrator bias, as well as changes in population size.
}
