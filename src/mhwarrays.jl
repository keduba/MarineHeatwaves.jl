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
       msst = view(sst, mhwix)
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

function clim(input::Matrix, drange)
    clima = [vec(mean(input[:,vcat(drange[i]...)], dims = 2)) for i in eachindex(drange)]
    clima[60] = mean((clima[59], clima[61]))
    return clima
end

clim(input::Vector, drange) = [mean(findices(input, dd)) for dd in drange]
tresh(input::Vector, drange, thresh) = [quantile(findices(input, dd), thresh) for dd in drange]

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
    return mstts, mends, hna
end

function startindices(mstts, hna)
    mstartsxs =  ones(length(hna)+1)
    mstartsxs[1] = first(mstts)
    mstartsxs[2:end] = mstts[2:end][hna] .+ 1
    return mstartsxs
end

function endindices(mends, hna)
    mendsxs =  ones(length(hna)+1)
    mendsxs[end] = last(mends)
    mendsxs[begin:end-1] = mends[begin:end-1][hna]
    return mendsxs
end


"""
    This function calculates the anomalies for the temperature array and the clim/threshold.
"""
function anoms(sst, clim, thsh, lyd, mst, mse) 
    anom = sst[mst:mse] - clim[lyd[mst:mse]]
    atcd = thsh[lyd[mst:mse]] - clim[lyd[mst:mse]]
    anbf = sst[mst-1] - clim[lyd[mst-1]]
    anft = sst[mse+1] - clim[lyd[mse+1]]
    lnxt = lastindex(sst)
    return (; anom, atcd, anbf, anft, lnxt)
end

categorys(eanoms) = nanmin(4, nanmaximum(fld.(eanoms.anom, eanoms.atcd)))

function ronset(eanoms, mst)
    anom, anbf = eanoms.anom, eanoms.anbf
    nmx, fan, ngx = anomfn(anom), first(anom), argfn(anom)
    lnmx, snmx = nmx - nanmean((fan, anbf)), nmx - fan
    ron = mst > 1 ? /(lnmx, (ngx + 0.5)) : /(snmx, ngx)
    return ron
end

function rdecline(eanoms, mse) 
    anft, lnx, nmx, lan, ngx, lnt = (eanoms.anft, eanoms.lnxt, anomfn(eanoms.anom), last(eanoms.anom), argfn(eanoms.anom), length(eanoms.anom))
    lnmx, snmx = (nmx - nanmean((lan, anft)), nmx - lan)
    rdc = mse < lnx ? /(lnmx, (lnt - ngx + 0.5)) : ngx == lnx ? /(snmx, 1.0) : /(snmx, (lnt - ngx))
    return rdc
end


function _eventmetrics(sst, clim, thrs, lyd, evdates, mstarts, mends)
    evanom, cats, rons, rdcs = ntuple(_ -> Vector{eltype(sst)}(undef, length(mstarts)), 4)
    
    mhwout, catout = (zeros(size(sst)) for _ in 1:2)

    for (m, (mst, mse)) in enumerate(zip(mstarts, mends))
        eanoms = anoms(sst, clim, thrs, lyd, mst, mse)
        evanom[m] = eanoms.anom
        cats[m] = categorys(eanoms)
        rons[m] = ronset(eanoms, mst)
        rdcs[m] = rdecline(eanoms, mse)
        mhwout[mst:mse] = eanoms.anom
        catout[mst:mse] .= cats[m]
    end
    return mhwout, catout, evanom, rons, rdcs, cats
end

"""
    Version for single vector which could be part of a loop if we do the values of one pixel all the way to the end.

    eventmetrics(sst, climthresh, startendsindices) -> Vector{Vector{T}}[events, mhwout, mhwcat]
"""
eventmetrics(msst::MCTemp, mstartsxs, mendsxs) = _eventmetrics(msst.temp, msst.clima, msst.thresh, msst.lyday, msst.dates, mstartsxs, mendsxs) 

# review from here.

function _meanmetrics(events, metrics, fullyears)
    lfy = length(fullyears)
    # metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maximum, :days, :frequency)
    meanmets = map(nanmean, (events.meananom, events.cumanom, events.ronset, events.rdecline, events.durations))
    meanmax = nanmaximum(events.maxanom)
    meandays = /(nansum(events.durations), lfy)
    meanfreq = /(length(events.durations), lfy)
    mmets =  NamedTuple{metrics}((meanmets..., meanmax, meandays, meanfreq))
    return mmets
end

function annualmetrics(events, metrics, fullyears)
    lfy = length(fullyears)
    evyears = events.alldates .|> year
    annuals = ntuple(_ -> zeros(lfy), length(metrics))
    for ey in eachindex(fullyears)
        yearixs = findall(isequal(fullyears[ey]), evyears)
        annuals[1][ey] = nanmean(events.meananom[yearixs]) # nanmean(anom[yearixs])
        annuals[2][ey] = nanmean(events.cumanom[yearixs])  # nanmean(nansum(anom[yearixs]))
        annuals[3][ey] = nanmean(events.ronset[yearixs])
        annuals[4][ey] = nanmean(events.rdecline[yearixs])
        annuals[5][ey] = nanmean(events.durations[yearixs])
        annuals[6][ey] = nanmaximum(events.maxanom[yearixs])
        annuals[7][ey] = nansum(events.durations[yearixs])
        annuals[8][ey] = length(events.durations[yearixs])
    end
end

function mhctemp(sst, sstdate::StepRange{Date, Day}, mdate::StepRange{Date, Day}, cdate::StepRange{Date, Day}; threshold=threshold)
    mhsst, mask, mlyd = subtemp(sst, sstdate, mdate)
    clsst, mask, clyd = subtemp(sst, sstdate, cdate)
    dvec = daterange(clyd, winwidth)
    clima, climq = clthr(clsst, dvec, threshold)
    return mhsst, mdate, mlyd, mask, clima, climq
end

function MHTemp(sst, sstdate, mdate, cdate; threshold=0.9)
    mhwin = MHTemp(mhctemp(sst, sstdate, mdate, cdate, threshold), argmax, maximum)
    mexcd = exceed(mhwin)
    excdarray = mexcd[1]
    mstarts = startlabel(mexcd)
    mends = endlabel(mexcd)
    mhxs = _indices(mstarts, mends, mindur, maxgap)
    msxs = startindices(mhxs.mstts, mhxs.hna)
    mexs = endindices(mhxs.mends, mhxs.hna)
    mhwout, catout, pxanoms, r2cat = eventmetrics(mhwin, msxs, mexs)
end
