using Base.Iterators: flatten

abstract type MarineHeatWave end

struct MHTemp{T <: AbstractFloat, Ti <: Integer} <: MarineHeatWave
    temp::VecOrMat{T}
    dates::StepRange{Date, Day}
    lyday::Vector{Ti}
    mask
    clima::VecOrMat
    thresh::VecOrMat
    argfn::Function# = nanargmax
    anomfn::Function # = nanmaximum
end

struct MCTemp{T <: AbstractFloat, Ti <: Integer} <: MarineHeatWave
    temp::VecOrMat{T}
    dates::StepRange{Date, Day}
    lyday::Vector{Ti}
    mask
    clima::VecOrMat
    thresh::VecOrMat
    excfn::Function# = greater/less than
    argfn::Function# = nanargmax
    anomfn::Function # = nanmaximum
end


struct MarineHW{T, N} <: MarineHeatWave
    temp::Array{T, N}
    category::Array{T, N}
    exceed::Array{Union{Missing, Bool}, N}
    annuals::NamedTuple
    means::NamedTuple
    trends::NamedTuple
    pvalues::NamedTuple
    pmetrics::NamedTuple
end


const winwidth = 5
const pctwidth = 31
const mindur = 5
const maxgap = 2

Base.length(fl::Base.Iterators.Flatten{<:AbstractVector{<:UnitRange}}) = sum(length, fl.it)

Base.IteratorSize(::Base.Iterators.Flatten{<:AbstractVector{<:UnitRange}}) = Base.HasLength()

"""
    seamask(sst) -> BitVecOrMat

Returns only sea areas. Optional dimension. For arrays, expected dimension is the 3rd dimension.
"""
seamask(xy, dims=3) = dropdims(count(!isnan, xy, dims=dims) .> 0, dims=dims)


"""
    _subtemp(sst::Array{T, N}, sstindex) -> subsst, mask indices

    This function is used to subset the sst to the indices for the mhwdate and/or climatology date. It's not exported.
"""
function _subtemp(sst::Vector, mhwix) 
   xz = seamask(sst, 1)[1]
   if xz
        msst = sst[mhwix]
       return msst, xz
   else
       @warn "The temperature data is all NaNs. Please check again."
   end
end

function _subtemp(sst::Array, mhwix)
    xz = seamask(sst)
    xyz = CartesianIndices(xz)[xz]
    msst = sst[xyz, mhwix]
    return msst, xz, xyz
end

function subtemp(sst, sstdate, evdate)
    etix, elyd = timeindices(sstdate, evdate)
    evsst, mask = _subtemp(sst, etix)
    return evsst, mask, elyd
end

function findices(av, dr)
    (av[i] for i in flatten(dr))
end

function tresh(input::Matrix, drange, thresh)
    inp = deepcopy(input)
    outq = [quantile!.(eachrow(inp[:,vcat(drange[d]...)]), thresh) for d in eachindex(drange)]
    outq[60] = mean((outq[59], outq[61]))
    return outq
end

function clim(input::Vector, drange)
   return [mean(findices(input, dd)) for dd in drange]
end

function tresh(input::Vector, drange, thresh)
    return [quantile(findices(input, dd), thresh) for dd in drange]
end

function clim(input::Matrix, drange)
    clima = [vec(mean(input[:,vcat(drange[i]...)], dims = 2)) for i in eachindex(drange)]
    clima[60] = mean((clima[59], clima[61]))
    return clima
end



function clthr(input::Vector, drange, thresh)
    clima = clim(input, drange)
    climq = tresh(input, drange, thresh)
    _smoothdata!(clima)
    _smoothdata!(climq)
    return clima, climq
end

function clthr(input::Matrix, drange, thresh)
    clima = reduce(hcat, clim(input, drange))
    climq = reduce(hcat, tresh(input, drange, thresh))
    _smoothdata!.(eachrow(clima))
    _smoothdata!.(eachrow(climq))
    return clima, climq
end


"""
    _smoothdata!(ctarray::Union{Vector, SubArray}, pw) -> ctarray

    This function modifies `ctarray` in place. `ctarray` is the climatology/threshold array. It calculates and returns the moving mean. `pw` is the smoothing window. By default, this value is `31`.
    This is a test: lastindex(B) == 3N || error("The smoothing failed. Please check the `_smoothdata!` function.")
"""
function _smoothdata!(ctarray, pw=pctwidth) 
    B = repeat(ctarray, 3)
    N = length(ctarray)
    ctarray[begin:end] = view(movmean(B, pw), N+1:2N)
    return ctarray
end

_exceed(exfn, tempdata::Vector, threshdata::Vector, lyd) = exfn.(tempdata, threshdata[lyd])

_exceed(exfn, tempdata::Matrix, threshdata::Matrix, lyd) = exfn.(tempdata, threshdata[:, lyd])

function exceed(x::MCTemp)
    mexc = _exceed(x.excfn, x.temp, x.thresh, x.lyday)
    return mexc
end
""" 
S == E || error("Something's not right with the starts and ends in the exceedance vector.")

exceed(x::MCTemp) -> exceedance::VecOrMat, excdiff::VecOrMat, 
Label the starts and ends using the exceedance and their differences.
"""
#=function exceed(x::MCTemp)
    mexc = _exceed(x.excfn, x.temp, x.thresh, x.lyday)
    N = ndims(x.temp)
    mexd = N == 1 ? diff(mexc) : diff(mexc, dims = N)
    if N == 1
        mexd = diff(mexc)
        ss = ifelse(isone(first(mexc)), count(==(1), mexd) + 1, count(==(1), mexd))
        es = ifelse(isone(last(mexc)), count(==(-1), mexd) + 1, count(==(-1), mexd))
    else 
        mexd = diff(mexc, dims = N)
        fvec = first.(eachrow(mexc))
        lvec = last.(eachrow(mexc))
    
        ss = [count(==(1), row) for row in eachrow(mexd)]
        es = [count(==(-1), row) for row in eachrow(mexd)]
        
        ss[isone.(fvec)] .+= 1
        es[isone.(lvec)] .+= 1
    end
    return mexc, mexd, es, ss
end

function endlabel(exceed)
    _, mexd, E, = exceed
    menders = ones(Int, E)
    L = length(mexd) + 1
    menders[end] = L

    for (i, (k, _)) in enumerate(Iterators.filter(p -> isequal(-1, p.second), pairs(mexd)))
        menders[i] = k
    end
   return menders
end

function startlabel(exceed)
    mexc, mexd, _, S = exceed
    mstarts = ones(Int, S)
    
    for (i, (k, _)) in enumerate(Iterators.filter(p -> isone(p.second), pairs(mexd)))
        isone(first(mexc)) ? mstarts[i+1] = Int(k) : mstarts[i] = Int(k)
    end
    return mstarts
end


function startlabel(exceed) 
    mexc, mexd, _, S = exceed
    mstarts = [ones(Int, s) for s in S]
    fvec = first.(eachrow(mexc))
    fsts = mstarts[isone.(fvec)]
    ntfs = mstarts[.!isone.(fvec)] 
    nfsts = eachrow(mexd)[.!isone.(first.(eachrow(kedar[1])))]
    fstss = eachrow(mexd)[isone.(first.(eachrow(kedar[1])))]
    for (v, vrow) in ntfs
         vrow[begin:end] = Int[i for (i, a) in enumerate(nfsts[v]) if isone(a)]
    end
    for (v, vrow) in fstss
         vrow[begin+1:end] = Int[i for (i, a) in enumerate(fstss[v]) if isone(a)]
    end
    return mstarts
end

function endlabel(exceed) 
 _, mexd, E, _ = exceed
    mends = [ones(Int, e) for e in E]
    l = size(mexd, ndims(mexd)) + 1
    lvec = last.(eachrow(mexc))
    for (v, vrow) in mends
        vrow[begin:end] = Int[i for (i, a) in enumerate(eachrow(mexd)[v]) if isequal(-1, a)]
     end
    return mends
end
=#


function __mylabeling(mexc, mexcd)
        S = ifelse(isone(first(mexc)), count(==(1), mexcd) + 1, count(==(1), mexcd))
        E = ifelse(isone(last(mexc)), count(==(-1), mexcd) + 1, count(==(-1), mexcd))
        S == E || error("Something's not right with the starts and ends in the exceedance vector.")
        mstarts, menders = ntuple(_ -> ones(Int, S), 2)
        menders[end] = lastindex(mexc)
        for (i, (k,)) in enumerate(Iterators.filter(p -> isone(p.second), pairs(mexcd)))
                isone(first(mexc)) ? mstarts[i+1] = k : mstarts[i] = k
        end
        for (i, (k,)) in enumerate(Iterators.filter(p -> isequal(-1, p.second), pairs(mexcd)))
                menders[i] = k
        end
        return mstarts, menders
end

function _mylabeling(mexc::BitVector)
    mexcd = diff(mexc)
    mstarts, menders = __mylabeling(mexc, mexcd)
    return mstarts, menders
end

function _mylabeling(mexc::BitMatrix) 
    mexcd = diff(mexc, dims=2)
    mstarts, menders = ([Int[] for _ in axes(mexc, 1)] for _ in 1:2)
    for (m, (colc, cold)) in enumerate(zip(eachrow(mexc), eachrow(mexcd)))
        mstarts[m], menders[m] = __mylabeling(colc, cold)
    end
    mstarts, menders
end

"""
    The indices in the vector/matrix where the events begin and end.

_indices(sts, ends, minduration, maximumgap) -> (stixs, enixs)::Tuple{2, Vector{Integer}}
"""
#=function _indices(mstarts, menders, mindur, maxgap)
    oldurations = menders - mstarts
    mstts, mends = (ix[oldurations .≥ mindur] for ix in (mstarts, menders))
    hna = mstts[2:end] - mends[1:end-1] .> maxgap
    return (; mstts, mends, hna)
end

function startindices(mstts, hna)
    mstartsxs = ones(Int, sum(hna)+1)
    mstartsxs[1] = first(mstts)
    mstartsxs[2:end] = mstts[2:end][hna] .+ 1
    return mstartsxs
end

function endindices(mends, hna)
    mendsxs = ones(Int, sum(hna)+1)
    mendsxs[end] = last(mends)
    mendsxs[begin:end-1] = mends[begin:end-1][hna]
    return mendsxs
end
=#

function _indices(mstarts, menders, mindur, maxgap)
    oldurations = menders - mstarts
    mstts = [ix[od .≥ mindur] for (ix, od) in zip(mstarts, oldurations)]
    mends = [ix[od .≥ mindur] for (ix, od) in zip(menders, oldurations)]
    hna = [(ms[2:end] - me[1:end-1] .> maxgap) for (ms, me) in zip(mstts, mends)]
    mstartsxs = [vcat(ms[2:end][hn] .+ 1) for (ms, hn) in zip(mstts, hna)]
    mendersxs = [vcat(me[begin:end-1][hn] .+ 1) for (me, hn) in zip(mends, hna)]
    for (mst, mse, ms, me) in zip(mstartsxs, mendersxs, mstts, mends)
        pushfirst!(mst, first(ms))
        push!(mse, last(me))
    end
    mstartsxs, mendersxs
end

"""
    This function calculates the anomalies for the temperature array and the clim/threshold.
"""
function _anoms(sst, clim, thsh, lyd, mst, mse) 
    anom = sst[mst:mse] - clim[lyd[mst:mse]]
    thsd = thsh[lyd[mst:mse]] - clim[lyd[mst:mse]]
    anbf = sst[mst-1] - clim[lyd[mst-1]]
    anft = sst[mse+1] - clim[lyd[mse+1]]
    lnxt = lastindex(sst)
    return (; anom, anbf, anft, lnxt, thsd)
end


_categorys(anoms) = min(4, maximum(fld.(anoms.anom, anoms.thsd)))

function ronset(anoms, mst)
    anom, anbf = anoms.anom, anoms.anbf
    nmx, fan, ngx = maximum(anom), first(anom), argmax(anom)
    lnmx, snmx = nmx - mean((fan, anbf)), nmx - fan
    ron = mst > 1 ? /(lnmx, (ngx + 0.5)) : /(snmx, ngx)
    return ron
end

function rdecline(anoms, mse) 
    anft, lnx, nmx, lan, ngx, lnt = (anoms.anft, anoms.lnxt, maximum(anoms.anom), last(anoms.anom), argmax(anoms.anom), length(anoms.anom))
    lnmx, snmx = (nmx - mean((lan, anft)), nmx - lan)
    rdc = mse < lnx ? /(lnmx, (lnt - ngx + 0.5)) : ngx == lnx ? /(snmx, 1.0) : /(snmx, (lnt - ngx))
    return rdc
end


function _eventmetrics(sst, clim, thrs, lyd, evdate, mstarts, mends)
    l = length(mstarts)
    cats, rons, rdcs = ntuple(_ -> Vector{eltype(sst)}(undef, l), 3)
    evanom = Vector{Vector{eltype(sst)}}(undef, l)
    mhwout, catout = (zeros(size(sst)) for _ in 1:2)
    stdate, endate = ntuple(_ -> Vector{Date}(undef, l), 2) 
    for (m, (mst, mse)) in enumerate(zip(mstarts, mends))
        eanoms = _anoms(sst, clim, thrs, lyd, mst, mse)
        evanom[m] = eanoms.anom
        cats[m] = _categorys(eanoms)
        rons[m] = ronset(eanoms, mst)
        rdcs[m] = rdecline(eanoms, mse)
        stdate[m] = evdate[mst]
        endate[m] = evdate[mse]
        mhwout[mst:mse] = eanoms.anom
        catout[mst:mse] .= cats[m]
    end
    return evanom, rons, rdcs, cats, stdate, endate, mhwout, catout
end


"""
    Version for single vector which could be part of a loop if we do the values of one pixel all the way to the end.

    eventmetrics(sst, climthresh, startendsindices) -> Vector{Vector{T}}[events, mhwout, mhwcat]
"""
eventmetrics(msst::MCTemp, mstartsxs, mendsxs) = _eventmetrics(msst.temp, msst.clima, msst.thresh, msst.lyday, msst.dates, mstartsxs, mendsxs) 


function eventmetricsa(msst::MCTemp, mstartsxs, mendsxs)
    outvectors = ntuple(_ -> Vector(undef, length(mstartsxs)), 8)
    
    for (m, (sstv, clim, thrs)) in enumerate(zip(eachrow(msst.temp), eachrow(msst.clima), eachrow(msst.thresh)))
        innervec = _eventmetrics(sstv, clim, thrs, msst.lyday, msst.dates, mstartsxs[m], mendsxs[m])
        outvectors[1][m] = innervec[1] # events
        outvectors[2][m] = innervec[2] # mhwout
        outvectors[3][m] = innervec[3] # catout
        outvectors[4][m] = innervec[4] # cats
        outvectors[5][m] = innervec[5] # rons
        outvectors[6][m] = innervec[6] # rdcs
        outvectors[7][m] = innervec[7] # rdcs
        outvectors[8][m] = innervec[8] # rdcs
    end

    return outvectors
end
function meanmetrics(evanom, rons, rdec, fullyears)
    lfy = length(fullyears)
     metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maximum, :days, :frequency)
    # Per Event Metrics
     meanint = mean.(evanom)
     cumint = sum.(evanom)
     maxint = maximum.(evanom)
     duration = length.(evanom)
     varint = std.(evanom)
    evmets = (; meanint, cumint, maxint, duration, varint)

    # Overall Mean Metrics
     mmeanint = mean(Iterators.flatten(evanom))
     mcumint = mean(cumint)
     mmaxint = maximum(Iterators.flatten(evanom))
     mronset = mean(rons)
     mrdecline = mean(rdec)
     mduration = mean(length.(evanom))
     mdays = sum(duration) / lfy
     mfrequency = length(duration) / lfy
    
    meanmets = (; mmeanint, mcumint, mmaxint,
        mronset, mrdecline, mduration, 
        mdays, mfrequency)

    return meanmets, evmets
end


function annualmetrics(evanom,ronset, rdecline,  stdate, endate, fullyears)
    lfy = length(fullyears)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maximum, :days, :frequency)
    evyears = collect(flatten(map((x,y) -> x:y, stdate, endate))) .|> year
    styr, enyr = year.(stdate), year.(endate)
    evanomf = collect(flatten(evanom))
    annuals = ntuple(_ -> zeros(lfy), length(metrics))
    for (ey, yr) in enumerate(fullyears)
        yearixs = findall(isequal(yr), evyears)
        yrix = [i for (i, (a, b)) in enumerate(zip(styr, enyr)) if a == yr || b == yr]
        if isempty(yearixs)
            continue
        end
        annuals[1][ey] = mean(evanomf[yearixs]) # correct nanmean(anom[yearixs])
        annuals[2][ey] = mean(sum.(evanom[yrix]))  #  nanmean(nansum(anom[yearixs]))
        annuals[3][ey] = mean(ronset[yrix])
        annuals[4][ey] = mean(rdecline[yrix])
        annuals[5][ey] = mean(length.(evanom[yrix])) # duration
        annuals[6][ey] = maximum(evanomf[yearixs]) # maxint
        annuals[7][ey] = length(evanomf[yearixs]) # days
        annuals[8][ey] = length(yrix) # frequency
    end
    return NamedTuple{metrics}(annuals)
end

 function __trendv(metric)
    trend, pval, pmet = ntuple(_ -> ones(1), 3)
    m = size(metric, 1)
    X = 1:m
    y = replace(metric, NaN => missing)
    data = DataFrame(X=X, Y=y)
    lr = lm(@formula(Y ~ X), data)
    trend[1] = coef(lr)[2]
    pval[1] = coeftable(lr).cols[4][2]
    pmet[1] = Float64(prod(confint(lr)[2, :]) ≥ 0)
    return trend[1], pval[1], pmet[1]
 end


function trends(annualmetrics)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maximum, :days, :frequency)
    trds = (:trends, :pvalues, :pmetrics)
    m = length(annualmetrics)
    tss = ntuple(_ -> ones(m), 3)
    for (m, metric) in enumerate(annualmetrics)
        tss[1][m], tss[2][m], tss[3][m] = __trendv(metric)
    end
    tss = (NamedTuple{metrics}(x) for x in tss)
    return (; zip(trds, tss)...)
end

function matevents(msst::MCTemp, mstartsxs, mendsxs)

    evmets, meanmets, annmets, trendmets = ntuple(_ -> Vector(undef, length(mstartsxs)), 4)

    fullyears = unique(year.(msst.dates))

     for i in eachindex(mstartsxs)
        evmets[i] = _eventmetrics(eachrow(msst.temp)[i], eachrow(msst.clima)[i], eachrow(msst.thresh)[i], msst.lyday, msst.dates, mstartsxs[i], mendsxs[i])
        meanmets[i] = meanmetrics(evmets[i][1], evmets[i][2], evmets[i][3], fullyears)
        annmets[i] = annualmetrics(evmets[i][1],evmets[i][2], evmets[i][3], evmets[i][5], evmets[i][6], fullyears)
        trendmets = trends(annmets[i])
    end
    return evmets, meanmets, annmets, trendmets
end


#    

function mhctemp(sst, sstdate::StepRange{Date, Day}, mdate::StepRange{Date, Day}, cdate::StepRange{Date, Day}; threshold=threshold)
    mhsst, mask, mlyd = subtemp(sst, sstdate, mdate)
    clsst, mask, clyd = subtemp(sst, sstdate, cdate)
    dvec = daterange(clyd, winwidth)
    clima, climq = clthr(clsst, dvec, threshold)
    return mhsst, mdate, mlyd, mask, clima, climq
end

function edetect(sst::Vector, sstdate, mdate, cdate; threshold=0.9)
    excdfn = >=
    mhwin = MCTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    mexcd = exceed(mhwin)
    # mstarts = startlabel(mexcd)
    # mends = endlabel(mexcd)
    mstarts, mends = _mylabeling(mexcd)
    msxs, mexs = _indices(mstarts, mends, mindur, maxgap)
    # msxs = startindices(mhxs.mstts, mhxs.hna)
    # mexs = endindices(mhxs.mends, mhxs.hna)
    evanom, rons, rdcs, cats, stdate, endate, mhwtemp, cattemp = eventmetrics(mhwin, msxs, mexs)
    fullyears = unique(year.(mdate))
    mmets, evmets = meanmetrics(evanom, rons, rdcs, fullyears)
    anmets = annualmetrics(evanom, stdate, endate, rons, rdcs, fullyears)
    trmets = trends(anmets)
    return (mhtemp = mhwtemp, mhanoms = evanom, mhexd = mexcd, categories = cattemp, events = evmets, means = mmets, annuals = anmets, trend = trmets.trends, pvalues = trmets.pvalues, pmetrics = trmets.pmetrics, ronsets = rons, rdeclines = rdcs, startds = stdate, endds = endate)

end


function edetect(sst::Array, sstdate, mdate, cdate; threshold=0.9)
    excdfn = >=
    mhwin = MCTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    mexcd = exceed(mhwin)
    mstarts, mends = _mylabeling(mexcd)
    msxs, mexs = _indices(mstarts, mends, mindur, maxgap)
    evmets, meanmets, annmets, trendmets = matevents(mhwin, msxs, mexs)
    return evmets, meanmets, annmets, trendmets


    #(mhtemp = mhwtemp, mhanoms = evanom, mhexd = mexcd, categories = cattemp, events = evmets, means = mmets, annuals = anmets, trend = trmets.trends, pvalues = trmets.pvalues, pmetrics = trmets.pmetrics, ronsets = rons, rdeclines = rdcs, startds = stdate, endds = endate)

end

# function vecarr(sst, msst)
# N = ndims(sst)
# TN = Array{Union{T, Missing}, N}
# TB = Array{Union{Bool, Missing}, N}
# const NT = NamedTuple
# mask
# mhwtemp output: x,y,z
# x, y = size(sst)
# z = length(msst.dates)
# zz = length(unique(year.(msst.dates))) # for the annual metrics
#
# trds = (:means, :trends, :pvalues, :pmetrics)
# metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maximum, :days, :frequency)
# lt, lm = length(trds), length(metrics)

# mhwexd = TB(missing, dims) # for vector sst
# mhwtemp, mhwcat = ntuple(_ -> TN(missing,  ds), 2) # for sst as vector
# annuals = NT{metrics}(ntuple(_ -> TN(missing, zz), lm))
# mets = NT{trds}(ntuple(_ -> NT{metrics}(ntuple(_ -> TN(missing, 1), lm)), lt))

# mhwexd = TB(missing, ds, z)
# mhwtemp, mhwcat = ntuple(_ -> TN(missing, ds, z), 2)
# annuals = NT{metrics}(ntuple(_ -> TN(missing, ds, zz), lm))
# mets = NT{trds}(ntuple(_ -> NT{metrics}(ntuple(_ -> TN(missing, ds), lm)), lt))
