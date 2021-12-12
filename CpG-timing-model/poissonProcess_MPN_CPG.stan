data {
    int n; // number of data points for m (samples)
    int m[n]; // number of mutations per sample to MRCA of mutations in that sample
    int M; // number of mutations up until MRCA of MPN cells in both samples
    
    // prior parameters:
    real<lower=0> min_permitted_lower_bound;
    real<lower=0> max_permitted_lower_bound;
    real<lower=0> expon_rate;
    real TWIN_A_AGE;
    real TWIN_B_AGE;
}

parameters {
    real<lower=0> t_mrca;
    vector<lower=0, upper=40>[n] i_raw; // used to construct post-MRCA interval
    real<lower=0, upper=3> rate; // rate parameter of the Poisson process (mutations per genome per unit time [year])
}
transformed parameters {
	vector<lower=0>[n] t = i_raw - t_mrca; // i is the post-MRCA interval; then have no sampling statement for p in the model block, then you are implicitly assigning a uniform [0,1] prior on i; https://mc-stan.org/docs/2_18/stan-users-guide/some-differences-in-the-statistical-models-that-are-allowed.html
}
model {
    t_mrca ~ uniform(min_permitted_lower_bound, max_permitted_lower_bound);

    // prior on mutation rate is Exponential ( [0, inf) )
    // NB: Expectation(mutation rate) = 1/expon_rate
    rate ~ exponential(expon_rate);

    for (j in 1:n) {
        m[j] ~ poisson(rate * t[j]);
    }
    // likelihood of 'M' is Poisson, rate scaled by t_mrca
    M ~ poisson(rate * t_mrca);
}
