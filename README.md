# survELtest  <img src='man/figures/logo_new.png' align="right" width = "280" height="300" />
Computing the one-sided/two-sided integrated/maximally selected EL statistics for simultaneous testing, the one-sided/two-sided EL tests for pointwise testing, and an initial test that precedes one-sided testing to exclude the possibility of crossings or alternative orderings among the survival functions.
### Data Frame
- `hepatitis` Severe alcoholic hepatitis data
- `threearm` Time to first remission data
- `hazardcross` Simulated survival data with crossing hazard functions from the piece- wise exponential distribution
- `hazardcross_Weibull` Simulated survival data with crossing hazard functions from the Weibull distribution
### Function
- `intELtest` The integrated EL test
- `supELtest` The maximally selected EL test
- `ptwiseELtest` The pointwise EL testing
- `nocrossings` The test that excludes the possibility of crossings or alternative orderings among the survival functions
- `print.intELtest` Print an intELtest object
- `print.supELtest` Print a supELtest object
- `print.ptwiseELtest` Print a ptwiseELtest object
- `print.nocrossings` Print a nocrossings object
- `summary.intELtest` Summary function for intELtest object
- `summary.supELtest` Summary function for supELtest object
- `summary.ptwiseELtest` Summary function for ptwiseELtest object
- `summary.nocrossings` Summary function for nocrossings object
## Installation
``` r
# install package directly 
install.packages("survELtest")
library(survELtest)
```
### Development version
``` r
#install.packages("devtools", dependencies = TRUE)
devtools::install_github("news11/survELtest")
library(survELtest)
```
## Usage
- `intELtest(formula, data = NULL, group_order = NULL, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, wt = "p.event", alpha = 0.05, seed = 1011, nlimit = 200)`
- `supELtest(formula, data = NULL, group_order = NULL, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200)`
- `ptwiseELtest(formula, data = NULL, group_order = NULL, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200)`
- `nocrossings(formula, data = NULL, group_order = NULL, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, alpha = 0.05, seed = 1011, nlimit = 200)`
- `print.intELtest(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)`
- `print.supELtest(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)`
- `print.ptwiseELtest(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)`
- `print.nocrossings(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)`
- `summary.intELtest(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)`
- `summary.supELtest(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)`
- `summary.ptwiseELtest(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)`
- `summary.nocrossings(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)`
### Arguments
- `formula` a formula object with a Surv object as the response on the left of the ~ operator and the grouping variable as the term on the right. The Surv object involves two variables: the observed survival and censoring times, and the censoring indicator, which takes a value of 1 if the observed time is uncensored and 0 otherwise. The grouping variable takes different values for different groups.
- `data` an optional data frame containing the variables in the formula: the observed survival and censoring times, the censoring indicator, and the grouping variable. If not found in data, the variables in the formula should be already defined by the user or in attached R objects. The default is the data frame with three columns of variables taken from the formula: column 1 contains the observed survival and censoring times, column 2 the censoring indicator, and column 3 the grouping variable.
- `group_order` a k-vector containing the values of the grouping variable, with the j-th element being the group hypothesized to have the j-th highest survival rates, j = 1, . . . , k. The default is the vector of sorted grouping variables.
- `t1` the first endpoint of a prespecified time interval, if any, to which the comparison of the survival functions is restricted. The default value is 0.
- `t2` the second endpoint of a prespecified time interval, if any, to which the comparison of the survival functions is restricted. The default value is ∞.
- `sided` 2 if two-sided test, and 1 if one-sided test. The default value is 2.
- `nboot` the number of bootstrap replications in calculating critical values for the tests. The default value is 1000.
- `wt` the name of the weight for the integrated EL statistics: "p.event", "dF", or "dt". The default is "p.event".
- `alpha` the pre-specified significance level of the tests. The default value is 0.05.
- `seed` the seed for the random number generator in R, for generating bootstrap samples needed to calculate the critical values for the tests. The default value is 1011.
- `nlimit` a number used to calculate nsplit= m/nlimit, the number of parts into which the calculation of the nboot bootstrap replications is split. The use of this vari- able can make computation faster when the number of time points m is large. The default value for nlimit is 200.
- `x` the result of a call to the intELtest/supELtest/ptwiseELtest/nocrossings function.
- `object` the result of a call to the intELtest/supELtest/ptwiseELtest/nocrossings function.
- `digits` significant digits to print, the default value is max(3L,getOption("digits")-3L).
- `quiet` a logical indicating whether to reduce the amount of output or not, the default value is FALSE.
- `...` for future method.
## More Information
Find the reference manual for more details: https://cran.r-project.org/web/packages/survELtest/survELtest.pdf

