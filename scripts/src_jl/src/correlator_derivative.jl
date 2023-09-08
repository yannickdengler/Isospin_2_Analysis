"""
    correlator_derivative(c::AbstractArray;t_dim=1)

Perform a three-point symmetric numerical derivative of the correlator data in 
the Euclidean time direction. This is intended to remove constant terms in the
correlator.

Ba default, the temporal dimension of the array containing the correlator data 
`c` is assumed to be in the 1st dimension. This can be changed by providing an 
explicit dimension `t_dim` as a keyword argument
"""
function correlator_derivative(c::AbstractArray;t_dim=1)
    c_deriv = similar(c)
    # Get the number of Euclidean timeslices in the specified dimension
    T = size(c)[t_dim]
    # get the overall number of dimensions for contructing a unit vector in the
    # temporal direction.
    n = ndims(c)
    δt(i) = ifelse(i == t_dim,1,0)
    # ntuple applies the function in the first argument δt to every index `1:n`
    # Thus, this corresponds to a delta-step in the temporal direction, i.e. a 
    # unit vector in Euclidean time. 
    unit_t = CartesianIndex(ntuple(δt,n))
    # Use julia's Cartesian indices to index arrays of arbitrary dimension 
    for i in CartesianIndices(c)
        # special case the first and the last index: 
        if i[t_dim] == 1
           c_deriv[i] = c[i+unit_t] - c[i]
        elseif i[t_dim] == T
           c_deriv[i] = c[i] - c[i-unit_t]
        else
           c_deriv[i] = (c[i+unit_t] - c[i-unit_t])/2
        end    
    end
    return c_deriv
end
