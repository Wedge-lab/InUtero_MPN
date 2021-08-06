data {
    int n; // number of data points for m (samples)
    int m[n]; // number of mutations per sample in short time period (t)
    int M; // number of mutations since origin of MPN, until MRCA of all samples (long time period, T)
    
    // prior parameters:
   
    real<lower=0> expon_rate;
}
transformed data {
    real TWIN_A_AGE = 1992 / 52; // age of TC1 at sampling in years
    real TWIN_B_AGE = 2044 / 52; // age of TC2 at sampling in years
}
parameters {
    real<lower=0> t_mrca;
         // upper=TWIN_B_AGE> t_mrca; // time to MRCA of both samples 
    
    vector<lower=0,
           upper=TWIN_B_AGE>[n] t; // length of time (years) before present to each sample's MRCA
    
    real<lower=1.5, upper=77> rate; // rate parameter of the Poisson process (mutations per unit time [year])
   
}

model {
   t_mrca ~ exponential(0.029) T[, 0.769]; // 95% of prior weight on 40 weeks of gestation

    for (j in 1:n) {
    	t[j] ~ exponential(0.029) T[0.769, 35];  // 95% of prior on period from birth to age of diagnosis
    }

    // prior on mutation rate is Exponential ( [0, inf) )
    // NB: Expectation(mutation rate) = 1/expon_rate
    rate ~ exponential(expon_rate);
 
    m[1] ~ poisson(rate * t[1]);
    m[2] ~ poisson(rate * t[2]);
    // likelihood of 'M' is Poisson, rate scaled by 'T'
    M ~ poisson(rate * t_mrca);
}
