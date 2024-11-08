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

climv(input::Vector, drange) = [mean(findices(input, dd)) for dd in drange]

treshv(input::Vector, drange, thresh) = [quantile(findices(input, dd), thresh) for dd in drange]

function clim(input::Matrix, drange)
    clima = [vec(mean(input[:,vcat(drange[i]...)], dims = 2)) for i in eachindex(drange)]
    clima[60] = mean((clima[59], clima[61]))
    return clima
end



function clthr(input::Vector, drange, thresh)
    clima = climv(input, drange)
    climq = treshv(input, drange, thresh)
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
""" 
S == E || error("Something's not right with the starts and ends in the exceedance vector.")

exceed(x::MCTemp) -> exceedance::VecOrMat, excdiff::VecOrMat, 
Label the starts and ends using the exceedance and their differences.
"""
function exceed(x::MCTemp)
    mexc = _exceed(x.excfn, x.temp, x.thresh, x.lyday)
    N = ndims(x.temp)
    mexd = N == 1 ? diff(mexc) : diff(mexc, dims = N)
    S = ifelse(isone(first(mexc)), count(==(1), mexd) + 1, count(==(1), mexd))
    E = ifelse(isone(last(mexc)), count(==(-1), mexd) + 1, count(==(-1), mexd))
   return mexc, mexd, E, S
end

function endlabel(exceed)
    _, mexd, E, = exceed
    menders = ones(Int, E)
    L = length(mexd) + 1
    menders[end] = L

    for (i, (k,)) in enumerate(Iterators.filter(p -> isequal(-1, p.second), pairs(mexd)))
        menders[i] = k
    end
   return menders
end

function startlabel(exceed)
    mexc, mexd, _, S = exceed
   
    mstarts = ones(Int, S)
    
    for (i, (k,)) in enumerate(Iterators.filter(p -> isone(p.second), pairs(mexd)))
        isone(first(mexc)) ? mstarts[i+1] = Int(k) : mstarts[i] = Int(k)
    end
    return mstarts
end

"""
    The indices in the vector/matrix where the events begin and end.

_indices(sts, ends, minduration, maximumgap) -> (stixs, enixs)::Tuple{2, Vector{Integer}}
"""
function _indices(mstarts, menders, mindur, maxgap)
    oldurations = menders - mstarts
    mstts, mends = (ix[oldurations .≥ mindur] for ix in (mstarts, menders))
    hna = mstts[2:end] - mends[1:end-1] .> maxgap
    return (; mstts, mends, hna)
end

function startindices(mstts, hna)
    mstartsxs =  ones(Int, sum(hna)+1)
    mstartsxs[1] = first(mstts)
    mstartsxs[2:end] = mstts[2:end][hna] .+ 1
    return mstartsxs
end

function endindices(mends, hna)
    mendsxs =  ones(Int, sum(hna)+1)
    mendsxs[end] = last(mends)
    mendsxs[begin:end-1] = mends[begin:end-1][hna]
    return mendsxs
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

# _thsd(thsh, clim, lyd, mst, mse) = thsh[lyd[mst:mse]] - clim[lyd[mst:mse]]

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
    stdate, endate = ntuple(_ -> Vector{Date}(undef, l), 2) 
    for (m, (mst, mse)) in enumerate(zip(mstarts, mends))
        eanoms = _anoms(sst, clim, thrs, lyd, mst, mse)
        evanom[m] = eanoms.anom
        cats[m] = _categorys(eanoms)
        rons[m] = ronset(eanoms, mst)
        rdcs[m] = rdecline(eanoms, mse)
        stdate[m] = evdate[mst]
        endate[m] = evdate[mse]
    end
    return evanom, rons, rdcs, cats, stdate, endate
end


"""
    Version for single vector which could be part of a loop if we do the values of one pixel all the way to the end.

    eventmetrics(sst, climthresh, startendsindices) -> Vector{Vector{T}}[events, mhwout, mhwcat]
"""
eventmetrics(msst::MCTemp, mstartsxs, mendsxs) = _eventmetrics(msst.temp, msst.clima, msst.thresh, msst.lyday, msst.dates, mstartsxs, mendsxs) 

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
        mronset, mrdecline, mduration, mdays)

    return meanmets, evmets
end


function annualmetrics(evanom, stdate, endate, ronset, rdecline, fullyears)
    lfy = length(fullyears)
    eyrs = stdate .|> year 
     metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maximum, :days, :frequency)
    evyears = collect(flatten(map((x,y) -> x:y, stdate, endate))) .|> year
    evanom = collect(flatten(evanom))
    annuals = ntuple(_ -> zeros(lfy), length(metrics))
    for ey in eachindex(fullyears)
        yearixs = findall(isequal(fullyears[ey]), evyears)
        yrix = findall(isequal(fullyears[ey]), eyrs)
        if isempty(yearixs)
            continue
        end
        annuals[1][ey] = mean(evanom[yearixs]) # nanmean(anom[yearixs])
        annuals[2][ey] = mean(sum.(evanom[yearixs]))  # nanmean(nansum(anom[yearixs]))
        annuals[3][ey] = mean(ronset[yrix])
        annuals[4][ey] = mean(rdecline[yrix])
        annuals[5][ey] = mean(length.(evanom[yearixs])) # duration
        annuals[6][ey] = maximum(evanom[yearixs]) # maxint
        annuals[7][ey] = sum(length.(evanom[yearixs])) # days
        annuals[8][ey] = length(evanom[yearixs]) # frequency
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



function mhctemp(sst, sstdate::StepRange{Date, Day}, mdate::StepRange{Date, Day}, cdate::StepRange{Date, Day}; threshold=threshold)
    mhsst, mask, mlyd = subtemp(sst, sstdate, mdate)
    clsst, mask, clyd = subtemp(sst, sstdate, cdate)
    dvec = daterange(clyd, winwidth)
    clima, climq = clthr(clsst, dvec, threshold)
    return mhsst, mdate, mlyd, mask, clima, climq
end

function edetect(sst, sstdate, mdate, cdate; threshold=0.9)
    excdfn = >=
    mhwin = MCTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    mexcd = exceed(mhwin)
    excdarray = mexcd[1]
    mstarts = startlabel(mexcd)
    mends = endlabel(mexcd)
    mhxs = _indices(mstarts, mends, mindur, maxgap)
    msxs = startindices(mhxs.mstts, mhxs.hna)
    mexs = endindices(mhxs.mends, mhxs.hna)
    evanom, rons, rdcs, cats, stdate, endate = eventmetrics(mhwin, msxs, mexs)
    fullyears = unique(year.(mdate))
    mmets, evmets = meanmetrics(evanom, rons, rdcs, fullyears)
    @show stdate, endate
    anmets = annualmetrics(evanom, stdate, endate, rons, rdcs, fullyears)
    trmets = trends(anmets)
    return (mhanoms = evanom, mhexd = excdarray, categories = cats, events = evmets, means = mmets, annuals = anmets, trend = trmets.trends, pvalues = trmets.pvalues, pmetrics = trmets.pmetrics, ronsets = rons, rdeclines = rdcs, startds = stdate, endds = endate)

end
