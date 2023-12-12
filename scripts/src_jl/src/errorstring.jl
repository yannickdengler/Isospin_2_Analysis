latticestring(Λ) = "$(Λ[1])×$(Λ[2])³"
decimals(s) = length(last(split(s,".")))
function dig(x)
    int = Int(round(abs(x)))
    dec = abs(x) - int
    d = 0
    if int > 0
        while int > 1
            d   -= 1
            int /= 10
        end
        return d
    else
        # print two significant digits for error in 10-29
        while dec < 3
            d   += 1
            dec *= 10
        end
        return d
    end
end

function errorstring_old(x,Δx)
    (isnan(x) || isnan(Δx)) && return "NaN"
    (isinf(x) || isinf(Δx)) && return "Inf"
    e = round(Δx,sigdigits=2)
    d = dig(e)
    if d > 0
        r = string(round(x;digits=d))
    else
        r = string(Int(round(x)))
    end
    if d > 0
        # make sure trailing 0 are printed if significant
        while length(last(split(r,"."))) < d
            r = r*"0"
        end
    end
    if d > 0
        err = string(round(Int,e*10.0^abs(d)))
        s = last(err) == '0' ? r[1:end-1]*"("*err[1:end-1]*")" : r*"("*err*")"
        return s
    else
        s = r*"($(Int(round(e))))"
        return s
    end
end
function leading_zeros(string)
    n_zeros = 0
    for i in string
        i=='0' && (n_zeros +=1)
        i!='0' && break
    end
    return n_zeros
end
# CASE A: 2 < uncertainty
function errorstring_bigger_two(a,Δa;nsig=2)
    central = Integer(round(a,digits=0))
    uncertainty = Integer(round(Δa,digits=0))
    res = "$central($uncertainty)"
    return res
end
# CASE B: 1 < uncertainty < 2
function errorstring_smaller_two(a,Δa;nsig=2)
    central = round(a,digits=nsig-1)
    uncertainty = round(Δa,digits=nsig-1)
    res = "$central($uncertainty)"
    return res
end
# CASE C: uncertainty < 1
function errorstring_smaller_one(a,Δa;nsig=2,roundthreshold=3*10^(nsig-1))
    # 1) take fractional part 
    # 2) drop the decimal point and the leading zero
    # 3) count the number of leading zeros 
    # 4) pad on the right with zeros
    Δstring = string(Δa)
    fractionalΔstring =  Δstring[3:end]
    n0 = leading_zeros(fractionalΔstring)
    fractionalΔstring = rpad(fractionalΔstring,n0+1+nsig,'0')
    # round the uncertainty correctly:
    # 1) take the number of nsig digits plus one and parse as integer
    # 2) divide by 10 and then round
    tmp = parse(Int,fractionalΔstring[n0+1:n0+1+nsig])
    uncertainty = round(Int,tmp/10)
    # if the uncertainty exceeds the roundthreshold and we have more than one
    # significant digit, then remove one significant digit
    # 1) find closest multiple of 10 by dividing by ten and rounding to the
    #    closest integer.
    # 2) reduce nsig by one for counting the number of digits of the main value
    if nsig > 1 && uncertainty >= roundthreshold
        uncertainty = Integer(round(uncertainty/10))
        nsig = nsig - 1
    end 
    # In total we have n0 + nsig digits to display after the decimal point
    # 1) Start with the central value and split it into integer and fractional part
    # 2) then take (n0 + nsig) decimal points and stitch them together
    # 3) Add uncertainty in brackets
    int, frac = split(string(a),'.') 
    res = int*'.'*frac[1:n0+nsig]*"($uncertainty)"
    return res
end
function errorstring(a,Δa)
    (isnan(a) || isnan(Δa)) && return "NaN"
    (isinf(a) || isinf(Δa)) && return "Inf"
    0.0001 < a < 10.0^6 || return errorstring_old(a,Δa)   
    0.0001 < Δa < 10.0^6 || return errorstring_old(a,Δa)   
    #0.0001 < a < 10.0^6 || error("Bracket notation: range not supported, value = $a")   
    #0.0001 < Δa < 10.0^6 || error("Bracket notation: range not supported, value = $Δa")   
    (Δa >= 2) && (return errorstring_bigger_two(a,Δa;nsig=2))
    (Δa >= 1) && (return errorstring_smaller_two(a,Δa;nsig=2))
    return errorstring_smaller_one(a,Δa;nsig=2)
end

