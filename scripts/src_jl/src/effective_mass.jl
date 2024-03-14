"""
    meff(c::AbstractVector;sign=+1)

Calculates the effective mass of the correlator `c` according to equation (10)
of arXiv:1607.06654. By default, a periodic correlator is assumed. For an 
anti-periodic correlator choose `sign=-1`.

    meff(c::AbstractVector,Δc::AbstractVector;sign=+1)

Gives a crude estimate of the uncertainty of the effective mass, when the 
standard uncertainty `Δc` is provided. For a more reliable result use the 
jackknife method.
"""
function meff(c::AbstractVector;sign=+1)
    T = length(c)
    m = similar(c)
    for t in 1:T
        m[t] = _meff_at_t(c,t,T;sign=sign)
    end
    return m
end
function meff(c::AbstractVector,Δc::AbstractVector;sign=+1)
    m1 = meff(c;sign=sign)
    m2 = meff(c + Δc;sign=sign)
    m3 = meff(c - Δc;sign=sign)
    Δm = @. abs(m3-m2)/2
    return m1, Δm
end
# see equation (10) in arXiv:1607.06654 [hep-lat]
# (Notation in Gattringer/Lang is misleading!)
function _meff_at_t(c::AbstractVector,t,T;sign=+1)
    # non-implicit mass as initial guess
    t0 = mod1(t+1,T)
    r0 = c[t0]/c[t]
    m0 = log(abs(r0))
    # correlator at large times (dropped overall factor)
    cor_lt(m,T,t) = exp(-m*t) + sign*exp(-m*(T-t))
    g(m,T,t,t0) = cor_lt(m,T,t0)/cor_lt(m,T,t) -r0
    # Use the more simpler algorithms from the Roots.jl package
    # find_zero() has more overhead and fails if the algorithm does not converged
    # Here we just use two simple, derivative free methods. If they do not converge
    # they return NaN. If that is the case then we try a slightly more robust algorithm.
    m = Roots.secant_method(x->g(x,T,t,t0),m0;maxevals=5000)
    if isnan(m)
       m = Roots.dfree(x->g(x,T,t,t0),m0)
    end
    return abs(m)
end
"""
    meff_jackknife(c::AbstractArray;sign=+1)

Calculates the effective mass of the correlator `c` according to equation (10)
of arXiv:1607.06654 and provides an estimator of the uncertainty using a 
jackknife analysis. 

By default, a periodic correlator is assumed. For an  anti-periodic correlator 
choose `sign=-1`.

The data `c` is ssumed to be a matrix where the first index corresponds to the
Euclidean time and the second one to the Monte-Carlo samples. 
"""
function meff_jackknife(corrs::AbstractMatrix;sign=+1)
    T, N = size(corrs)
    # create arrays for decay constant
    corrs_delete1 = zeros(T,N-1)
    m_eff = zeros(N,T)
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
        m_eff[i,:] = meff(C;sign)
    end
    return apply_jackknife(m_eff)
end
# This function takes care of arrays with additional inidces. The additional 
# indices are collected as ind by the use of the `axes(corr)` function. This 
# returns a set of iterators for each additional dimension of the array. We then
# loop over every additional index using the `Iterators.product` utility 
function meff_jackknife(corrs::AbstractArray;sign=+1)
    d = ndims(corrs)
    it, iMC, ind... = axes(corrs)
    size_meff_array = (size(corrs,1),size(corrs)[3:end]...)
    m_eff  = zeros(eltype(corrs),size_meff_array)
    Δm_eff = zeros(eltype(corrs),size_meff_array)
    #@showprogress "Jackknife analysis for effective mass:" for i in Iterators.product(ind...)
    for i in Iterators.product(ind...)
        # slurping is needed to correctly insert the tuple i into an index
        c = @view corrs[:,:,i...]
        m_eff[:,i...], Δmeff[:,i...] = meff_jackknife(c;sign)
    end
    return m_eff, Δm_eff
end
function apply_jackknife(obs)
    N  = size(obs)[1]
    O  = vec(nanmean(obs,dims=1))
    ΔO = vec(sqrt(N-1)*nanstd(obs,corrected=false,dims=1))
    return O, ΔO
end