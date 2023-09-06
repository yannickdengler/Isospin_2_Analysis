# see equation (10) in arXiv:1607.06654 [hep-lat]
# (Notation in Gattringer/Lang is misleading!)
function _meff_at_t(c::AbstractVector,t,T;sign=+1)
    # non-implicit mass as initial guess
    m0 = 1 # log(abs(c[t]/c[mod1(t+1,T)]))
    t0 = mod1(t-1,T)
    r0 = abs(c[t0]/c[t])
    # correlator at large times (dropped overall factor)
    cor_lt(m,T,t) = exp(-m*t) + sign*exp(-m*(t-T/2))
    # function to fit the effective mass
    g(m,T,t,t0) = cor_lt(m,T,t0)/cor_lt(m,T,t)
    # Use the more simpler algorithms from the Roots.jl package
    # find_zero() has more overhead and fails if the algorithm does not converged
    # Here we just use two simple, derivative free methods. If they do not converge
    # they return NaN. If that is the case then we try a slightly more robust algorithm.
    m = Roots.secant_method(x->g(x,T,t,t0)-r0,m0;maxevals=5000)
    if isnan(m)
       m = Roots.dfree(x->g(x,T,t,t0)-r0,m0)
    end
    return m
end
function implicit_meff(c::AbstractVector;sign=+1)
    T = length(c)
    m = zeros(eltype(c),div(T,2))
    for t in 1:div(T,2)
        m[t] = _meff_at_t(c,t,T;sign=sign)
    end
    return m
end
function implicit_meff(c::AbstractVector,Δc::AbstractVector;sign=+1)
    m1 = implicit_meff(c;sign=sign)
    m2 = implicit_meff(c + Δc;sign=sign)
    m3 = implicit_meff(c - Δc;sign=sign)
    Δm = @. abs(m3-m2)/2
    return m1, Δm
end
function implicit_meff_jackknife(corrs::AbstractMatrix;sign=+1)
    T, N = size(corrs)
    # create arrays for decay constant
    corrs_delete1 = zeros(T,N-1)
    meff = zeros(N,T÷2)
    # set up jack-knife (one deletion)
    for i in 1:N
        for t in 1:T
            for j in 1:N
               (j < i) && (corrs_delete1[t,j]   = corrs[t,j])
               (j > i) && (corrs_delete1[t,j-1] = corrs[t,j])
            end
        end
        # perform averaging for fitting weights
        C = reshape(mean(corrs_delete1,dims=2),T)
        meff[i,:] = implicit_meff(C;sign)
    end
    return apply_jackknife(meff)
end