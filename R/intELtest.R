#' The integrated EL test
#' 
#' @description \code{intELtest} gives a class of integrated EL statistics:
#' \deqn{\sum_{i=1}^{m}w_i\cdot \{-2\log R(t_i)\},}
#' where \eqn{R(t)} is the EL ratio that compares the survival functions at each given time \eqn{t},
#' \eqn{w_i} is the weight at each \eqn{t_i}, and \eqn{ 0<t_1<\ldots<t_m<\infty} are the (ordered)
#' observed uncensored times at which the Kaplan--Meier estimate is positive and less than \eqn{1} for
#' each sample.
#' @name intELtest
#' @param formula a formula object with a \code{Surv} object as the response on the left of the \code{~} operator 
#' and the grouping variable as the term on the right. The \code{Surv} object involves two variables: the 
#' observed survival and censoring times, and the censoring indicator (it takes a value of \eqn{1} if the 
#' observed time is uncensored and \eqn{0} otherwise). The grouping variable takes different values for different groups. 
#' @param data an optional data frame containing the variables in the \code{formula}: the observed survival and censoring 
#' times, the censoring indicator, and the grouping variable. If not found in \code{data}, the variables in the \code{formula} 
#' should be already defined by the user or in attached \code{R} objects. The default is the data frame with 3 columns of 
#' variables taken from the \code{formula}: column 1 contains the observed survival and censoring times, column 2 the 
#' censoring indicator, and column 3 the grouping variable.
#' @param formula a formula object with a \code{Surv} object as the response on the left of the \code{~} operator and 
#' the grouping variable as the term on the right. The \code{Surv} object involves two variables: the observed survival 
#' and censoring times, and the censoring indicator, which takes a value of \eqn{1} if the observed time is uncensored 
#' and \eqn{0} otherwise. The grouping variable takes different values for different groups. 
#' @param data an optional data frame containing the variables in the \code{formula}: the observed survival and censoring times,
#'  the censoring indicator, and the grouping variable. If not found in \code{data}, the variables in the \code{formula} should 
#'  be already defined by the user or in attached \code{R} objects. The default is the data frame with three columns of 
#'  variables taken from the formula: column 1 contains the observed survival and censoring times, column 2 the censoring 
#'  indicator, and column 3 the grouping variable.
#' @param group_order a \eqn{k}-vector containing the values of the grouping variable, with the \eqn{j}-th element being the group
#' hypothesized to have the \eqn{j}-th highest survival rates, \eqn{j=1,\ldots,k}. The default is the vector of sorted grouping variables.
#' @param t1  the first endpoint of a prespecified time interval, if any, to which the comparison of the survival functions is restricted. 
#' The default value is \eqn{0}.
#' @param t2 the second endpoint of a prespecified time interval, if any, to which the comparison of the survival 
#' functions is restricted. The default value is \eqn{\infty}.
#' @param sided \eqn{2} if two-sided test, and \eqn{1} if one-sided test. The default value is \eqn{2}.
#' @param nboot the number of bootstrap replications in calculating critical values for the tests. 
#' The default value is \eqn{1000}.
#' @param wt the name of the weight for the integrated EL statistics:
#' \code{"p.event"}, \code{"dF"}, or \code{"dt"}. The default is \code{"p.event"}.
#' @param alpha the pre-specified significance level of the tests. The default value is \eqn{0.05}.
#' @param seed  the seed for the random number generator in \code{R}, for generating bootstrap samples needed 
#' to calculate the critical values for the tests. The default value is \eqn{1011}.
#' @param nlimit a number used to calculate \code{nsplit=} \eqn{m}\code{/nlimit}, the number of 
#' parts into which the calculation of the \code{nboot} bootstrap replications is split. 
#' The use of this variable can make computation faster when the number of time points \eqn{m} is large. 
#' The default value for \code{nlimit} is 200.
#' @return \code{intELtest} returns a \code{intELtest} object, a list with 15 elements:
#' \itemize{
#'    \item \code{call} the function call
#'    \item \code{teststat} the resulting integrated EL statistics
#'    \item \code{critval} the critical value based on bootstrap
#'    \item \code{pvalue} the p-value of the test
#'    \item \code{formula} the value of the input argument of intELtest
#'    \item \code{data} the value of the input argument of intELtest
#'    \item \code{group_order} the value of the input argument of intELtest
#'    \item \code{t1} the value of the input argument of intELtest
#'    \item \code{t2} the value of the input argument of intELtest
#'    \item \code{sided} the value of the input argument of intELtest
#'    \item \code{nboot} the value of the input argument of intELtest
#'    \item \code{wt} the value of the input argument of intELtest
#'    \item \code{alpha} the value of the input argument of intELtest
#'    \item \code{seed} the value of the input argument of intELtest
#'    \item \code{nlimit} the value of the input argument of intELtest
#' }
#' Methods defined for \code{intELtest} objects are provided for \code{print} and \code{summary}.
#' @details There are three options for the weight \eqn{w_i}:
#' \itemize{
#'     \item (\code{wt = "p.event"}) \cr
#'      This default option is an objective weight,
#'      \deqn{w_i=\frac{d_i}{n},}
#'      which assigns weight proportional to the number of events \eqn{d_i}
#'      at each observed uncensored time \eqn{t_i}. Here \eqn{n} is the total sample size.
#'    \item (\code{wt = "dF"}) \cr
#'      Inspired by the integral-type statistics considered in Barmi and McKeague (2013),
#'      another weigth function is
#'      \deqn{w_i= \hat{F}(t_i)-\hat{F}(t_{i-1}),}
#'      for \eqn{i=1,\ldots,m}, where \eqn{\hat{F}(t)=1-\hat{S}(t)}, \eqn{\hat{S}(t)} is the pooled KM estimator, and \eqn{t_0 \equiv 0}.
#'      This reduces to the objective weight when there is no censoring. The resulting \eqn{I_n} can be seen as an empirical
#'      version of the expected negative two times log EL ratio under \eqn{H_0}.
#'    \item (\code{wt = "dt"}) \cr
#'      Inspired by the integral-type statistics considered in Pepe and Fleming (1989), another weight function is
#'      \deqn{w_i= t_{i+1}-t_i,} 
#'      for \eqn{i=1,\ldots,m}, where \eqn{t_{m+1} \equiv t_{m}}. This gives more weight to the time intervals where 
#'      there are fewer observed uncensored times, but can be affected by extreme observations.
#' }
#' @references
#' \itemize{
#'    \item H. Chang, I.W. McKeague, "Nonparametric testing for multiple survival functions 
#'      with non-inferiority margins," \emph{Annals of Statistics}, Vol. 47, No. 1, pp. 205-232, (2019).
#'    \item M. S. Pepe and T. R. Fleming, "Weighted Kaplan-Meier
#'      Statistics: A Class of Distance Tests for Censored Survival Data," \emph{Biometrics},
#'      Vol. 45, No. 2, pp. 497-507 (1989).
#'     \url{https://www.jstor.org/stable/2531492?seq=1#page_scan_tab_contents}
#'    \item H. E. Barmi and I.W. McKeague, "Empirical likelihood-based tests
#'      for stochastic ordering," \emph{Bernoulli}, Vol. 19, No. 1, pp. 295-307 (2013).
#'      \url{https://projecteuclid.org/euclid.bj/1358531751}
#' }
#' @seealso \code{\link{hepatitis}}, \code{\link{supELtest}}, \code{\link{ptwiseELtest}}, \code{\link{nocrossings}}, \code{\link{print.intELtest}}, \code{\link{summary.intELtest}}
#' @examples
#' library(survELtest)
#' intELtest(survival::Surv(hepatitis$time, hepatitis$censor) ~ hepatitis$group)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## intELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
#' ##     hepatitis$group)
#' ## 
#' ## Two-sided integrated EL test statistic = 1.42, p = 0.007
#' @export

intELtest = function(formula, data = NULL, group_order = NULL, t1 = 0, t2 = Inf, sided = 2, nboot = 1000, wt = "p.event", alpha = 0.05, seed = 1011, nlimit = 200) {
    #check and process inputs
    call = match.call()
    options(warn = -1)
    checkInput_formula_data("intELtest", formula, data)
    group_order = processInput_group_order(formula, group_order)
    formula = processInput_formula(formula, data)
    data = processInput_data(formula, data)
    checkInputs(data, group_order, t1, t2, sided, nboot, wt, alpha, seed, nlimit)
    options(warn = 0)
    
    k = length(unique(data$group))
    if (k == 2){
        at_ts = neg2ELratio(formula, data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
    }else{
        at_ts = teststat(formula, data, group_order, t1, t2, sided, nboot, alpha, details.return = TRUE, seed, nlimit)
    }
    
    if (is.null(at_ts)) return (NULL)
    if (wt == "p.event") {
        critval  = at_ts$int_dbarNtEL_SOcrit
        teststat = at_ts$inttest_dbarNt
        pvalue   = at_ts$p_value_inttest_dbarNt
    } else if (wt == "dt") {
        critval  = at_ts$int_dtEL_SOcrit
        teststat = at_ts$inttest_dt
        pvalue   = at_ts$p_value_inttest_dt
    } else if (wt == "dF") {
        critval  = at_ts$int_dFEL_SOcrit
        teststat = at_ts$inttest_dF
        pvalue   = at_ts$p_value_inttest_dF
    }
    
    result = list(call = call, teststat = teststat, critval = critval, pvalue = pvalue, formula = formula, data = data, group_order = group_order, t1 = t1, t2 = t2, sided = sided, nboot = nboot, wt = wt, alpha = alpha, seed = seed, nlimit = nlimit)
    attr(result, "class") = "intELtest"
    return (result)
}
