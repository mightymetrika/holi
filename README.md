
<!-- README.md is generated from README.Rmd. Please edit that file -->

# holi

<!-- badges: start -->
<!-- badges: end -->

The goal of holi is to provide web applications for higher order
likelihood inference.

## Installation

You can install the released version of ‘holi’ from
[CRAN](https://CRAN.R-project.org):

``` r
install.packages("holi")
```

You can install the development version of holi from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mightymetrika/holi")
```

## Example

This is a basic example which shows you how to compare the p-value from
stats::glm() and the r\* p-value from holi::rstar_glm() when analyzing
‘mtcars’. The holi::rstar_glm() function relies on
likelihoodAsy::rstar().

``` r
library(holi)

# Fit model
rs_linear <- rstar_glm(mpg ~ wt + hp, .data = mtcars, .model = "linear")
#> get mle ....     get mle under the null.... 
#> start Monte Carlo computation 
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%

# See results from stats::glm()
rs_linear$fit_glm |> summary()
#> 
#> Call:
#> stats::glm(formula = .formula, family = stats::gaussian, data = .data)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 37.22727    1.59879  23.285  < 2e-16 ***
#> wt          -3.87783    0.63273  -6.129 1.12e-06 ***
#> hp          -0.03177    0.00903  -3.519  0.00145 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 6.725785)
#> 
#>     Null deviance: 1126.05  on 31  degrees of freedom
#> Residual deviance:  195.05  on 29  degrees of freedom
#> AIC: 156.65
#> 
#> Number of Fisher Scoring iterations: 2

# See r* results
rs_linear$fit_glm |> summary()
#> 
#> Call:
#> stats::glm(formula = .formula, family = stats::gaussian, data = .data)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 37.22727    1.59879  23.285  < 2e-16 ***
#> wt          -3.87783    0.63273  -6.129 1.12e-06 ***
#> hp          -0.03177    0.00903  -3.519  0.00145 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for gaussian family taken to be 6.725785)
#> 
#>     Null deviance: 1126.05  on 31  degrees of freedom
#> Residual deviance:  195.05  on 29  degrees of freedom
#> AIC: 156.65
#> 
#> Number of Fisher Scoring iterations: 2
```

In this example, the p-value for r\* (5.556e-07) is smaller than the
p-value for stats::glm() (1.12e-06).

## References

Pierce, D. A., & Bellio, R. (2017). Modern Likelihood-Frequentist
Inference. International Statistical Review / Revue Internationale de
Statistique, 85(3), 519–541. <doi:10.1111/insr.12232>
