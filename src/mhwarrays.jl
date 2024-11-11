using Base.Iterators: flatten

abstract type MarineHeatWave end

struct MHTemp{T <: AbstractFloat, Ti <: Integer} <: MarineHeatWave
    temp::VecOrMat{T}
    dates::StepRange{Date, Day}
    lyday::Vector{Ti}
    mask
    clima::VecOrMat
    thresh::VecOrMat
     excfn # = ≥
     argfn # = argmax
    anomfn # = maximum
end

struct MCTemp{T <: AbstractFloat, Ti <: Integer} <: MarineHeatWave
    temp::VecOrMat{T}
    dates::StepRange{Date, Day}
    lyday::Vector{Ti}
    mask
    clima::VecOrMat
    thresh::VecOrMat
    excfn # = ≤ 
    argfn # = argmin
    anomfn # = minimum
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
    return msst, xyz
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

function exceed(x::Union{MCTemp, MHTemp})
    return _exceed(x.excfn, x.temp, x.thresh, x.lyday)
end

""" 
exceed(x::MCTemp) -> exceedance::VecOrMat
Label the starts and ends using the exceedance and their differences.
"""
function endlabel(mexd)
   menders = Int[i for (i, m) in enumerate(mexd) if isequal(-1, m)]
   return menders
end

function startlabel(mexd)
    mstarts = Int[i for (i, m) in enumerate(mexd) if isone(m)]
    return mstarts
end

function _mylabeling(mexc::BitVector)
    mexd = diff(mexc)
    mstarts = startlabel(mexd)
    menders = endlabel(mexd)
    mstts, mends, hna = _indices(mstarts, menders, mindur, maxgap)
    mstartxs = startindices(mstts, hna)
    mendsxs = endindices(mends, hna)
    return mstartxs, mendsxs
end

function _mylabeling(mexc::BitMatrix) 
    mexd = diff(mexc, dims=2)
    ss = [count(==(1), row) for row in eachrow(mexd)]
    es = [count(==(-1), row) for row in eachrow(mexd)]
    ss == es || throw("Tengo frio! Something's off with the starts and ends in the `_mylabeling`")
    
    mstarts, menders = ntuple(_ -> [ones(Int, e) for e in es], 2)
    for (m, mxd) in enumerate(eachrow(mexd))
        mstarts[m] = startlabel(mxd)
        menders[m] = endlabel(mxd)
    end
    
    mstartxs, mendsxs = ntuple(_ -> typeof(mstarts)(undef, size(mexc, 1)), 2)
    for s in eachindex(mstarts, menders)
        mstts, mends, hna = _indices(mstarts[s], menders[s], mindur, maxgap)
        mstartxs[s] = startindices(mstts, hna)
        mendsxs[s] = endindices(mends, hna)
    end
    mstartxs, mendsxs
end

"""
    The indices in the vector/matrix where the events begin and end.

_indices(sts, ends, minduration, maximumgap) -> (stixs, enixs)::Tuple{2, Vector{Integer}}
"""
function _indices(mstarts, menders, mindur, maxgap)
    oldurations = menders - mstarts
    mstts, mends = (ix[oldurations .≥ mindur] for ix in (mstarts, menders))
    hna = mstts[2:end] - mends[1:end-1] .> maxgap
    return (mstts, mends, hna)
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

# DOUBLE CHECK THIS FOR THE COLD SPELL
_categorys(anoms) = min(4, maximum(fld.(anoms.anom, anoms.thsd)))

function ronset(anoms, mst, anomfn, argfn)
    anom, anbf = anoms.anom, anoms.anbf
    nmx, fan, ngx = anomfn(anom), first(anom), argfn(anom)
    lnmx, snmx = nmx - mean((fan, anbf)), nmx - fan
    ron = mst > 1 ? /(lnmx, (ngx + 0.5)) : /(snmx, ngx)
    return ron
end

function rdecline(anoms, mse, anomfn, argfn) 
    anft, lnx, nmx, lan, ngx, lnt = (anoms.anft, anoms.lnxt, anomfn(anoms.anom), last(anoms.anom), argfn(anoms.anom), length(anoms.anom))
    lnmx, snmx = (nmx - mean((lan, anft)), nmx - lan)
    rdc = mse < lnx ? /(lnmx, (lnt - ngx + 0.5)) : ngx == lnx ? /(snmx, 1.0) : /(snmx, (lnt - ngx))
    return rdc
end


function _eventmetrics(sst, clim, thrs, lyd, anomfn, argfn, evdate, mstarts, mends)
    l = length(mstarts)
    cats, rons, rdcs = ntuple(_ -> Vector{eltype(sst)}(undef, l), 3)
    evanom = Vector{Vector{eltype(sst)}}(undef, l)
    stdate, endate = ntuple(_ -> Vector{Date}(undef, l), 2) 
    mhwout, catout = (zeros(size(sst)) for _ in 1:2)
    for (m, (mst, mse)) in enumerate(zip(mstarts, mends))
        eanoms = _anoms(sst, clim, thrs, lyd, mst, mse)
        evanom[m] = eanoms.anom
        cats[m] = _categorys(eanoms)
        rons[m] = ronset(eanoms, mst, anomfn, argfn)
        rdcs[m] = rdecline(eanoms, mse, anomfn, argfn)
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
eventmetrics(msst::MCTemp, mstartsxs, mendsxs) = _eventmetrics(msst.temp, msst.clima, msst.thresh, msst.lyday, msst.anomfn, msst.argfn, msst.dates, mstartsxs, mendsxs) 

function meanmetrics(evanom, rons, rdec, fullyears)
    lfy = length(fullyears)
    # Per Event Metrics
     meanint = mean.(evanom)
     cumint = sum.(evanom)
     maxint = maximum.(evanom)
     duration = length.(evanom)
     varint = std.(evanom)
    evmets = (; meanint, cumint, maxint, duration, varint)

    # Overall Mean Metrics
     meanint = mean(Iterators.flatten(evanom))
     cumint = mean(cumint)
     maxint = maximum(Iterators.flatten(evanom))
     ronset = mean(rons)
     rdecline = mean(rdec)
     duration = mean(length.(evanom))
     days = sum(duration) / lfy
     frequency = length(duration) / lfy
    
    meanmets = (; meanint, cumint, maxint,
        ronset, rdecline, duration, 
        days, frequency)

    return meanmets, evmets
end

function annualmetrics(evanom, ronset, rdecline, stdate, endate, fullyears)
    lfy = length(fullyears)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maxint, :days, :frequency)
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
        annuals[1][ey] = mean(evanomf[yearixs])      # meanint
        annuals[2][ey] = mean(sum.(evanom[yrix]))    # cumint
        annuals[3][ey] = mean(ronset[yrix]) 
        annuals[4][ey] = mean(rdecline[yrix])
        annuals[5][ey] = mean(length.(evanom[yrix])) # duration
        annuals[6][ey] = maximum(evanomf[yearixs])   # maxint
        annuals[7][ey] = length(evanomf[yearixs])    # days
        annuals[8][ey] = length(yrix)                # frequency
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
    trds = (:trends, :pvalues, :pmetrics)
    m = length(annualmetrics)
    tss = ntuple(_ -> ones(m), 3)
    for (m, metric) in enumerate(annualmetrics)
        tss[1][m], tss[2][m], tss[3][m] = __trendv(metric)
    end
    tss = (NamedTuple{keys(annualmetrics)}(x) for x in tss)
    return (; zip(trds, tss)...)
end

function matevents(msms::MarineHW, msst::Union{MCTemp, MHTemp}, mexc, mstartsxs, mendsxs)
    evmets = Vector(undef, length(mstartsxs))
    fullyears = unique(year.(msst.dates))
    mask = msst.mask
    msms.exceed[mask, :] = mexc
     for i in eachindex(mstartsxs)
        evmeta = _eventmetrics(eachrow(msst.temp)[i], eachrow(msst.clima)[i], eachrow(msst.thresh)[i], msst.lyday, msst.anomfn, msst.argfn, msst.dates, mstartsxs[i], mendsxs[i])
        msms.temp[mask[i], :] = evmeta[7]
        msms.category[mask[i], :] = evmeta[8]
        meanmets, evmets[i] = meanmetrics(evmeta[1], evmeta[2], evmeta[3], fullyears)
        anmets = annualmetrics(evmeta[1], evmeta[2], evmeta[3], evmeta[5], evmeta[6], fullyears)
        trmets = trends(anmets)
        for m in keys(msms.annuals)
            msms.annuals[m][mask[i], :] = anmets[m]
            msms.means[m][mask[i]] = meanmets[m]
            msms.trends[m][mask[i]] = trmets.trends[m]
            msms.pvalues[m][mask[i]] = trmets.pvalues[m]
            msms.pmetrics[m][mask[i]] = trmets.pmetrics[m]
        end
    end
    return msms, evmets
end

function mhctemp(sst, sstdate::StepRange{Date, Day}, mdate::StepRange{Date, Day}, cdate::StepRange{Date, Day}; threshold=threshold)
    mhsst, mask, mlyd = subtemp(sst, sstdate, mdate)
    clsst, mask, clyd = subtemp(sst, sstdate, cdate)
    dvec = daterange(clyd, winwidth)
    clima, climq = clthr(clsst, dvec, threshold)
    return mhsst, mdate, mlyd, mask, clima, climq
end

function edetect(sst::Vector, sstdate, mdate, cdate; threshold=0.9)
   excdfn = ≥ 
    mhwin = MHTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    msms = vecarr(sst, mdate)
    mexc = exceed(mhwin)
    msxs, mexs = _mylabeling(mexc)
    evanom, rons, rdcs, cats, stdate, endate, mhwout, catout  = eventmetrics(mhwin, msxs, mexs)
    fullyears = unique(year.(mdate))
    meanmets, evmets = meanmetrics(evanom, rons, rdcs, fullyears)
    anmets = annualmetrics(evanom, rons, rdcs, stdate, endate, fullyears)
    trmets = trends(anmets)
    msms.exceed .= mexc
    msms.temp .= mhwout
    msms.category .= catout
    # Annuals
    for m in keys(msms.annuals)
        msms.annuals[m][begin:end] = anmets[m]
        msms.means[m] .= meanmets[m]
        msms.trends[m][1] = trmets.trends[m]
        msms.pvalues[m][1] = trmets.pvalues[m]
        msms.pmetrics[m][1] = trmets.pmetrics[m]
    end

    return msms, evmets
end

function edetect(sst::Array, sstdate, mdate, cdate; threshold=0.9)
    excdfn = >=
    mhwin = MHTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    msms = vecarr(sst, mdate)
    mexc = exceed(mhwin)
    msxs, mexs = _mylabeling(mexc)
    msms, evmets = matevents(msms, mhwin, mexc, msxs, mexs)
    return msms, evmets
end

function vecarr(sst, mhwdate)
    N = ndims(sst)
    N ∈ (1,3) || throw("Oh heat! The dimension of the `sst` should be 1 or 3, we got $(N) instead!")
    T = eltype(sst) 
    TN = Array{Union{T, Missing}, N}
    TM = Matrix{Union{T, Missing}}
    TB = Array{Union{Bool, Missing}, N}
    NT = NamedTuple
   
    trds = (:means, :trends, :pvalues, :pmetrics)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maxint, :days, :frequency)
    lt, lm = length(trds), length(metrics)
    x, y = ifelse(N == 3, size(sst), (1,1))
    z = length(mhwdate)
    zz = length(unique(year.(mhwdate))) 
    ds = N == 1 ? z : (x, y)

    mhwexd = N == 1 ? TB(missing, ds) : TB(missing, ds..., z)
    mhwtemp, mhwcat = N == 1 ? ntuple(_ -> TN(missing,  ds), 2) : ntuple(_ -> TN(missing, ds..., z), 2)
    annuals = N == 1 ? NT{metrics}(ntuple(_ -> TN(missing, zz), lm)) : NT{metrics}(ntuple(_ -> TN(missing, ds..., zz), lm))
    mets = N == 1 ? ntuple(_ -> NT{metrics}(ntuple(_ -> TN(missing, N), lm)), lt) : ntuple(_ -> NT{metrics}(ntuple(_ -> TM(missing, ds...), lm)), lt)
    return MarineHW(mhwtemp, mhwcat, mhwexd, annuals, mets...)
end

# eventtable(eventanoms)  
