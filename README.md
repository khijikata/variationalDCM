# variationalDCM

**variationalDCM** is an R package that performs recently-developed variational Bayesian inference for diagnostic classification models (DCMs), which are a family of statistical models for collecting, analyzing, and reporting diagnostic information in Education and Psychology.

To install this package, use the following commands.

``` r
if(!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("khijikata/variationalDCM")
```

## Models

The following groups of models are currently supported: 
- DINA model 
- DINO model 
- Multiple-choice-DINA model 
- Saturated DCM 
- Hidden Markov DCM

## Acknowledgements

This package was developed as part of the project supported by JST, PRESTO Grant Number JPMJPR21C3, Japan and JSPS KAKENHI Grant Number 21H00936.

## References

-   Hijikata, K., Oka, M., Yamaguchi, K., & Okada, K. (2023). variationalDCM: An R package for variational Bayesian inference in diagnostic classification models. PsyArXive. <https://psyarxiv.com/f2sqd>

-   Yamaguchi, K., & Okada, K. (2020). Variational Bayes Inference for the DINA Model. *Journal of Educational and Behavioral Statistics*. <https://doi.org/10.3102/1076998620911934>

-   Yamaguchi, K., & Okada, K. (2020). Variational Bayes Inference Algorithm for the Saturated Diagnostic Classification Model. *Psychometrika*. <https://doi.org/10.1007/s11336-020-09739-w>

-   Yamaguchi, K. (2020). Variational Bayesian inference for the multiple-choice DINA model. *Behaviormetrika*. <https://doi.org/10.1007/s41237-020-00104-w>

-   Yamaguchi, K., & Martinez, A. J. (2023). Variational Bayes inference for hidden Markov diagnostic classification models. *British Journal of Mathematical and Statistical Psychology*. <https://doi.org/10.1111/bmsp.12308>
