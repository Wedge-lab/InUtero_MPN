require(rstan)
require(bayesplot)

# Compile model
mpn_model <- rstan::stan_model(file = "MPN.stan", model_name = "MPN poisson process")

# List of model parameters + input data
dat <- list(
    n = 2,  # number of post-MRCA mutation estimates
    m = c(506, 697),  # post-MRCA mutation estimates
    M = 8,  # pre-MRCA mutation estimates (use array to pass a single value)
    min_permitted_lower_bound = 0,  # lower bound on uniform prior on 'T' in years
    max_permitted_upper_bound = 2044/52,  # upper bound on uniform prior on 'T' in years
    expon_rate = 0.0579  # parameter of exponential prior on mean mutation rate (Derived from source data for Extended Data Fig. 8 in Gerstung et al. The evolutionary history of 2,658 cancers. Nature 578, 122–128 (2020))
)

# Get MCMC samples from the model
trace <- rstan::sampling(mpn_model, data = dat, iter = 50000, chains = 4)

# List parameter estimates
summary(trace, prob = c(0.025, 0.975))$summary

# Plot 
bayesplot::mcmc_areas(as.array(trace), pars = "t_mrca", prob = .95) + ggplot2::labs(title = "Time to MPN MRCA (years)"); ggsave("timeToMRCA.pdf")
bayesplot::mcmc_areas(as.array(trace), pars = "rate", prob = .95) + ggplot2::labs(title = "Mutation rate (total substitutions per year)"); ggsave("mutationRate.pdf")
