require(rstan)
require(bayesplot)
library(bridgesampling)

my_model <- rstan::stan_model(file = "poissonProcess_MPN_CPG.stan", model_name = "MPN poisson process")

# List of model parameters + input data
dat <- list(
    n = 2,  # number of post-MRCA mutation estimates
    m = c(91, 90),  # post-MRCA mutation estimates
    M = 7,  # pre-MRCA mutation estimates (use array to pass a single value)
    min_permitted_lower_bound = 0,
    max_permitted_lower_bound = 37, #age of youngest twin at presentation
    expon_rate = 0.09,  # parameter of exponential prior on mean mutation rate 
    TWIN_A_AGE = 40, # age of TC1 at sampling in years
    TWIN_B_AGE = 38.48 # age of TC2 at sampling in years
)

# Get MCMC samples from the model
trace <- rstan::sampling(my_model, data = dat, iter = 100000, chains = 5, warmup=20000,seed=483892929)

# List parameter estimates
summary(trace, prob = c(0.025, 0.975))$summary

BS <- bridge_sampler(trace)
print(BS)
bridgesampling::error_measures(BS)

bayesplot::mcmc_areas(as.array(trace), pars = "t_mrca", prob = .95) + ggplot2::labs(title = "Time to MRCA of MPNs") + ggplot2::xlim(0, 5)
bayesplot::mcmc_areas(as.array(trace), pars = "t[1]", prob = .95) + ggplot2::labs(title = "t_TwinB")
bayesplot::mcmc_areas(as.array(trace), pars = "t[2]", prob = .95) + ggplot2::labs(title = "t_TwinA")
bayesplot::mcmc_areas(as.array(trace), pars = "rate", prob = .95, prob_outer=1) + ggplot2::labs(title = "Mutation rate") + ggplot2::xlim(2, 3)
#ggsave("mutationRate.pdf") + ggplot2::expand_limits(x=c(2, 5))
