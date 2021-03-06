#' Summary function for ptwiseELtest object
#' 
#' @description Returns a list with a data frame containing the observed uncensored time points, and the decisions, 
#' statistics, and critical values of the pointwise EL tests at those time points.
#' @name summary.ptwiseELtest
#' @param object the result of a call to the \code{ptwiseELtest} function
#' @param digits significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}
#' @param quiet a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}
#' @param ... for future method
#' @return \code{summary.ptwiseELtest} returns a list with following components:
#' \itemize{
#'    \item \code{call} the statement used to create the \code{ptwiseELtest} object
#'    \item \code{result_dataframe} a dataframe with \code{time_pts} in the first column, \code{decision} 
#'    in the second column, \code{stat_ptwise} in the third column and \code{critval_ptwise} in the fourth column.
#'    \itemize{
#'          \item \code{time_pts} a vector containing the observed uncensored time points at which the
#'          Kaplan—Meier estimate is positive and less than \eqn{1} for each sample.
#'          \item \code{decision} a vector containing the decisions of the pointwise EL tests at \code{time_pts}.
#'          The decision at each of \code{time_pts} is \eqn{1} for rejection of the null hypothesis that the survival
#'          functions are the same at the specific time point, and \eqn{0} otherwise.
#'          \item \code{stat_ptwise} a vector containing the pointwise EL statistics at \code{time_pts}.
#'          \item \code{critval_ptwise} a vector containing the critical values for pointwise EL testing
#'          at \code{time_pts}.
#'    }
#' }
#' @seealso \code{\link{hepatitis}}, \code{\link{ptwiseELtest}}, \code{\link{print.ptwiseELtest}}
#' @examples
#' library(survELtest)
#' result = ptwiseELtest(survival::Surv(hepatitis$time, hepatitis$censor)~
#'              hepatitis$group, sided = 1)
#' summary(result)
#' 
#' ## OUTPUT:
#' ## Call:
#' ## ptwiseELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~  
#' ##     hepatitis$group, sided = 1)
#' ## 
#' ##    time_pts decision stat_ptwise critval_ptwise
#' ## 1       5.2        0      0.3005          2.951
#' ## 2       9.7        0      0.0000          2.833
#' ## 3      12.9        0      0.1627          2.748
#' ## 4      14.0        0      0.6114          2.583
#' ## 5      14.9        0      2.0010          2.780
#' ## 6      15.7        1      3.7873          2.764
#' ## 7      18.0        1      3.0722          2.652
#' ## 8      18.9        0      1.8878          2.454
#' ## 9      19.2        1      2.5896          2.339
#' ## 10     19.7        0      1.6133          2.601
#' ## 11     20.0        0      2.2393          2.383
#' ## 12     21.7        1      3.6936          2.192
#' ## 13     24.0        1      4.5083          2.300
#' ## 14     24.9        1      5.3743          2.391
#' ## 15     26.0        1      6.2879          2.253
#' ## 16     26.9        1      9.2827          2.117
#' ## 17     27.8        1     10.3581          2.209
#' ## 18     28.0        1      6.9862          2.317
#' ## 19     30.0        1      7.9190          2.346
#' ## 20     31.2        1      6.5074          2.318
#' ## 21     32.1        1      4.9709          2.310
#' ## 22     34.1        1      5.7455          2.360
#' ## 23     36.1        1      6.5627          2.244
#' ## 24     44.9        1      5.4374          2.363
#' ## 25     45.2        1      6.2240          2.416
#' ## 26     47.8        1      7.0519          2.409
#' ## 27     54.1        1      7.9198          2.427
#' ## 28     54.9        1      6.7260          2.310
#' ## 29     58.1        1      7.5667          2.456
#' ## 30     59.8        1      7.2524          2.483
#' ## 31     63.2        1      6.1770          2.511
#' ## 32     70.4        1      5.2110          2.562
#' ## 33     76.1        1      4.3461          2.683
#' ## 34     80.1        1      3.5753          2.744
#' ## 35     81.3        1      2.8926          2.467
#' ## 36     82.1        0      2.2925          2.669
#' ## 37     90.1        1      2.7908          2.543
#' ## 38     92.1        0      2.2120          2.523
#' ## 39     95.0        0      1.7079          2.755
#' ## 40     99.0        0      2.1383          2.762
#' ## 41    108.2        0      2.6206          2.652
#' ## 42    109.9        1      3.1475          2.630
#' ## 43    117.0        0      2.5398          2.646
#' ## 44    148.8        1      3.0555          2.685
#' ## 45    153.1        0      2.4658          2.774
#' @export
## generic function to summary an S3 object of class "supELtest"
summary.ptwiseELtest = function(object, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...){
  if(quiet == FALSE){
    cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    print(object$result_dataframe, digits = digits)
  }
  result = list(call = object$call, result_dataframe = object$result_dataframe)
  invisible(result)
}
