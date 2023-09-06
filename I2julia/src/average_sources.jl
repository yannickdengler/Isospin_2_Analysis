"""
        average_sources(hdf_file) 

    Takes an hdf5 file according and extract the correlator information from it.
    Then an average over all stochastic sources is performed. Instead of the 
    filename and array containing the correlator data can also be provided.

    NOTE: This assumes that the correlator consist of a sum of single diagrams.
    In case of a product of multiple diagrams this averaging is not necessarily 
    correct. E.g. for disconnected meson contributions this does not yield the 
    correct Euclidean time dependence. 
"""
function average_sources(hdf_file::AbstractString) 
    file_id = h5open(hdf_file)
    corr = read(file_id,"correlators")
    nsrc = read(file_id,"N_hits")
    close(file_id)
    return average_sources(corr)
end
function average_sources(corr::AbstractArray) 
    # According to our specification this correlator needs to have 5 indices
    # corresponding to (t,tMC,nsrc,nops)
    @assert ndims(corr) == 4
    nsrc = size(corr)[3]
    source_averaged = mean(corr,dims=3)
    # this removes the source-dimension from the array
    source_averaged = dropdims(source_averaged,dims=3)
    return source_averaged, nsrc
end

"""
        average_configurations(corr::AbstractArray)

    Average over all configurations and calculate the estimated variance. This 
    assumes that the average over all stochastic sources has been performed.

        average_configurations(hdf_file::AbstractArray)

    Averages over the correlator data from a specified hdf5 files. Assumes that
    stochastic source average has been performed. 

    It is assumed that the configurations are thermalised and independent, i.e. 
    that no autocorrelation exists.
"""
function average_configurations(corr::AbstractArray) 
    # According to our specification this correlator needs to have 4 indices
    # corresponding to (t,tMC,nmom,nops)
    @assert ndims(corr) == 3
    ncfg = size(corr)[2]
    avg_corr  = mean(corr,dims=2)
    Δavg_corr = std(corr,dims=2)./sqrt(ncfg)
    # this removes the source-dimension from the array
    avg_corr  = dropdims(avg_corr,dims=2)
    Δavg_corr = dropdims(Δavg_corr,dims=2)
    return avg_corr, Δavg_corr, ncfg
end
function average_configurations(hdf_file::AbstractString) 
    file_id = h5open(hdf_file)
    corr = read(file_id,"correlators")
    close(file_id)
    # If the array is 4-dimensional we assume that the average over all sources 
    # has already been performed
    @assert ndims(corr) == 3 || ndims(corr) == 4
    if ndims(corr) == 3
        return average_configurations(corr)
    elseif ndims(corr) == 4
        source_averaged, nsrc = average_sources(corr)
        return average_configurations(source_averaged)
    end
end
