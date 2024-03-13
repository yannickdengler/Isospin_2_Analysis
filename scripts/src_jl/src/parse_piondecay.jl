function parse_spectrum(file,type;)
    T = latticesize(file)[1]
    corr = zeros(T) # preallocate array for parsing of correlator
    dict = Dict{String,Vector{Float64}}()
    dictarray = Dict{String,Vector{Float64}}[]
    conf0 = 0
    src0  = 0
    # keep track of position in file for progress meter
    p = Progress(countlines(file); dt=1, desc="Match $type: Progress:")
    for line in eachline(file)
        next!(p)
        if occursin(type,line)
            # get configuration number
            pos_num = findfirst('#',line)
            end_num = findnext(' ',line,pos_num)
            conf = parse(Int,line[pos_num+1:end_num-1])
            # find number of the source if available
            src = 0
            # find last '=' sign which separates values from Γ structure
            # TODO this does not work for momenta
            pos_eq = findlast('=',line)
            #key_st = findprev(' ',line,pos_eq)
            key_st = last(findfirst(type,line))+1
            key = line[key_st+1:pos_eq-1]
            # parse corrrelator values
            pos_0 = findnext(' ',line,pos_eq)
            for t in 1:T
                pos_1 = findnext(' ',line,pos_0+1)
                corr[t] = Parsers.parse(Float64,line[pos_0:pos_1])
                pos_0 = pos_1
            end
            dict[key] = copy(corr)
            conf0 = conf
            src0  = src
        end
        # If we only have one source at a time and possibly one configuration
        # at a time: the method used to separate distinct measurements fails. 
        # In this case the end of measurement on a given confiuration is 
        # signalled by a line that reads:
        # [MAIN][0]Configuration #N: analysed in [a sec b usec]
        if occursin("analysed",line)
            if !isempty(dict)
                push!(dictarray,dict)
                dict = Dict{String,Vector{Float64}}()
            end
        end
    end
    if !isempty(dict)
        push!(dictarray,dict)
    end
    return _reshape_connected(dictarray)
end
function latticesize(file)
    for line in eachline(file)
        if occursin("Global size is",line)
            pos  = last(findfirst("Global size is",line))+1
            sizestring  = lstrip(line[pos:end])
            latticesize = parse.(Int,split(sizestring,"x"))
            return latticesize
        end
    end
end
function _reshape_connected(dict)
    corrs = Dict{String,Array{Float64}}()
    for k in keys(dict[1])
        corrs[k] = permutedims(reduce(hcat,getindex.(dict,k)))
    end
    return corrs
end