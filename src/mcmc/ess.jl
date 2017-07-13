"""
Lag-`k` autocorrelation of `x` from a variogram. `v` is the variance
of `x`, used when supplied.
"""
function ρ(x, k, v=var(x))
    x1 = @view(x[1:(end-k)])
    x2 = @view(x[(1+k):end])
    V = sum((x1 .- x2).^2)/length(x1)
    1-V/(2*v)
end

"""
Factor `τ` for effective sample size and the last lag which was used
in the estimation.
"""
function ess_factor(x)
    v = var(x)
    invfactor = 1 + 2*ρ(x, 1, v)
    K = 2
    while K < length(x)-2
        increment = ρ(x, K, v) + ρ(x, K+1, v)
        if increment < 0
            break
        else
            K += 2
            invfactor += 2*increment
        end
    end
    length(x)/invfactor, K #1/invfactor, K
end
