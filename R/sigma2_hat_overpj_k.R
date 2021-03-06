
sigma2_hat_overpj_k = function(t, fit, fitls) {
    n = fit$n
    k = length(fitls)
    iter = 1:k
    out = matrix(0, nrow = length(t), ncol = k)
    for (i in 1:length(t)) {
        d = lapply(fitls, FUN = function (f) {
            f$n.event[f$time <= t[i] & f$n.event != 0]
        })
        r = lapply(fitls, FUN = function (f) {
            f$n.risk[f$time <= t[i] & f$n.event != 0]
        })
        out_func = function(f){
            return (n * sum(division00(d[[f]], r[[f]] * (r[[f]] - d[[f]]))))
        }
        out[i,] = mapply(out_func,iter)
    }
    return(out)
}
