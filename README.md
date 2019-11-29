# survELtest  <img src='man/figures/logo_new.png' align="right" height="0" />
This R package contains routines for computing the one-sided/two-sided integrated/maximally selected EL statistics for simultaneous testing, the one-sided/two-sided EL tests for pointwise testing, and an initial test that precedes one-sided testing to exclude the possibility of crossings or alternative orderings.
### Data Frame
- `hepatitis` obtained by digitizing the published Kaplan-Meier curves
- `threearm` obtained by perturbing the time-to-remission data from patients in a three-arm randomized clinical trial for	 			the treatment of major depression
### Function
- `intELtest` the integrated EL test
- `supELtest` the maximally selected EL test
- `ptwiseELtest` the pointwise EL testing
- `nocrossings` the test that excludes the possibility of crossings or alternative order- ings
		
## Installation
``` r
# install package directly 
install.packages("survELtest")
library(survELtest)
```
### Development version
``` r
#install.packages("devtools")
devtools::install_github("news11/survELtest")
library(survELtest)
```

## Usage
- `intELtest(data, group_order = sort(unique(data[, 3])), t1 = 0, t2 = Inf, sided = 2, nboot = 1000, wt = "p.event",
alpha = 0.05, seed = 1011, nlimit = 200)`
- `supELtest(data, group_order = sort(unique(data[, 3])), t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200)`
- `ptwiseELtest(data, group_order = sort(unique(data[, 3])), t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200)`
- `nocrossings(data, group_order = sort(unique(data[, 3])), t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200)`
### Arguments
- `data` a data frame/matrix with 3 columns: column 1 contains the observed survival and censoring times, column 2 the censoring indicator, and column 3 the group- ing variable. This is a compulsory input.
- `group_order` a k-vector containing the values of the grouping variable (in column 3 of the data frame/matrix), with the j-th element being the group hypothesized to have the j-th highest survival rates, j = 1, . . . , k. The default is the vector of sorted grouping variables.
- `t1` the first endpoint of a prespecified time interval, if any, during which the com- parison of the survival functions is restricted to. The default value is 0.
- `t2` the second endpoint of a prespecified time interval, if any, during which the comparison of the survival functions is restricted to. The default value is ∞.
- `sided` 2 if two-sided test, and 1 if one-sided test. The default value is 2.
- `nboot` the number of bootstrap replications in calculating critical values for the tests.
The defualt value is 1000.
- `wt` the name of the weight for the integrated EL statistics: "p.event", "dF", or
"dt". The default is "p.event".
- `alpha` the pre-specified significance level of the tests. The default value is 0.05.
- `seed` the seed of random number generation in R for generating bootstrap samples needed to calculate critical values for the tests. The default value is 1011.
- `nlimit` a number used to calculate nsplit= m/nlimit, the number of parts we split the calculation of the nboot bootstrap replications into. This can make computation faster when the number of time points m is too large. The default value for nlimit is 200.

## More Information
Find the reference manual for more details: https://cran.r-project.org/web/packages/survELtest/survELtest.pdf

