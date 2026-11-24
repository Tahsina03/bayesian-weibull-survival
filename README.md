# Bayesian Weibull Survival Model for Lung Cancer Data

A Bayesian Weibull AFT model for right-censored survival data. I wrote the Metropolis sampler myself in R, without Stan or JAGS, and compared the results with maximum likelihood.

## Question

In patients with advanced lung cancer, how are age, sex, and ECOG performance status related to survival time?

## Data

`lung` from the R `survival` package: 227 patients with complete data, 164 deaths.

## Method

- Weibull AFT model on the log-time scale: $\log T_i = x_i^\top \beta + \sigma W_i$
- Covariates: age (per 10 years, centred at 60), sex, ECOG status (0, 1, 2+)
- Weakly informative priors: $\beta_j \sim N(0, 10^2)$, $\log\sigma \sim N(0, 2^2)$
- Random-walk Metropolis, 4 chains, 20,000 iterations each (5,000 burn-in), acceptance rate about 0.29
- Convergence: all R-hat below 1.01, effective sample size about 2,900 to 3,100

![Trace plots](output/figures/trace_plots.png)

## Results

| Covariate | Time ratio (95% CrI) |
|---|---|
| Age, per 10 years | 0.93 (0.81, 1.06) |
| Female | 1.51 (1.18, 1.96) |
| ECOG 1 vs 0 | 0.74 (0.54, 0.99) |
| ECOG 2+ vs 0 | 0.50 (0.36, 0.71) |

- Women survive about 1.5 times longer than men with the same age and ECOG status.
- Patients with ECOG 2+ survive about half as long as those with ECOG 0.
- Age has a weak effect once sex and ECOG are in the model.
- The shape parameter is about 1.33, so the risk of death goes up over time.
- The posterior estimates are very close to the MLE.

Posterior median survival for a 60-year-old ranges from about 6 months (man, ECOG 2+) to about 19 months (woman, ECOG 0).

![Posterior survival curves](output/figures/posterior_survival.png)

## How to run

Needs R 4.1+ with `survival`, `coda`, and `ggplot2`.

```bash
Rscript R/analysis.R
```

Results are saved in `output/`. The seed is fixed.

## References

- Ibrahim, J. G., Chen, M.-H., & Sinha, D. (2001). *Bayesian Survival Analysis*. Springer.
- Roberts, G. O., Gelman, A., & Gilks, W. R. (1997). Weak convergence and optimal scaling of random walk Metropolis algorithms. *Annals of Applied Probability*, 7, 110–120.

## Author

Tahsina Fariha Ohi, MS Statistics, Ball State University
