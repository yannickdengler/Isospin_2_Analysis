"""
        write_ensemble_list(hdf5list) 

    Take an iterable `hdf5list` containing the absolute, full filenames of the
    hdf5 files and creates a csv table of all corresponding  ensembles. 
"""
function write_ensemble_list(h5file;outdir="output/tables",filename="ensembles.csv") 
    # open file for writing list of ensembles
    ispath(outdir) || mkpath(outdir)
    file = joinpath(outdir,filename)
    io = open(file,"w")
    write(io,"β,m,L,T,nsrc,ncfg,name,tmin,tmax\n");

    fid = h5open(h5file,"r")

    # loop over all relevant
    for ensemble in keys(fid)
        
        T  = read(fid,joinpath(ensemble,"N_T"))[1]
        L  = read(fid,joinpath(ensemble,"N_L"))[1]
        β  = read(fid,joinpath(ensemble,"beta"))[1]
        m1 = read(fid,joinpath(ensemble,"m_1"))[1]
        m2 = read(fid,joinpath(ensemble,"m_2"))[1]

        corr = read(fid,joinpath(ensemble,"correlators"))
        nsrc = size(corr)[3]
        ncfg = size(corr)[2]

        # we always assume degenerate fermions
        @assert m1 == m2
        write(io,"$β,$m1,$L,$T,$nsrc,$ncfg,$(ensemble)/,1,$(T÷2)\n");
    end
    close(fid)
    close(io)
    sort_ensemble_file!(file)
end
function sort_ensemble_file!(filename)
    data, header = readdlm(filename,';',Any,header=true)
    data = sortslices(data,dims=1)
    writedlm(filename,vcat(header,data),';')
end
