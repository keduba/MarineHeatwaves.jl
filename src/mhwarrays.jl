using Base.Iterators: flatten

abstract type MarineHeatWave end

struct MHTemp{T<:AbstractFloat,Ti<:Integer} <: MarineHeatWave
    temp::VecOrMat{T}
    dates::StepRange{Date,Day}
    lyday::Vector{Ti}
    mask
    clima::VecOrMat
    thresh::VecOrMat
    excfn # = ≥
    argfn # = argmax
    anomfn # = maximum
end

struct MCTemp{T<:AbstractFloat,Ti<:Integer} <: MarineHeatWave
    temp::VecOrMat{T}
    dates::StepRange{Date,Day}
    lyday::Vector{Ti}
    mask
    clima::VecOrMat
    thresh::VecOrMat
    excfn # = ≤ 
    argfn # = argmin
    anomfn # = minimum
end


struct MarineHW{T,N} <: MarineHeatWave
    temp::Array{T,N}
    category::Array{T,N}
    exceed::Array{Union{Missing,Bool},N}
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

function matevents(msms::MarineHW, msst::MarineHeatWave, mexc::BitMatrix, mstartsxs, mendsxs)
    mni, cmi, mxi, drt, cats, rons, rdcs, vri, sdt, edt, cls, rws = ntuple(_ -> Vector(undef, length(mstartsxs)), 12)
    fullyears = unique(year.(msst.dates))
    mask = msst.mask
    x, y = size(msms.temp)
    climout, threshout = ntuple(_ -> similar(msms.temp, x, y, 366), 2)
    msms.exceed[mask, :] = mexc
    for i in eachindex(mstartsxs)
        evmeta = _eventmetrics(eachrow(msst.temp)[i], eachrow(msst.clima)[i], eachrow(msst.thresh)[i], msst.lyday, msst.anomfn, msst.argfn, msst.dates, mstartsxs[i], mendsxs[i])
        cats[i] = evmeta[4]
        msms.temp[mask[i], :] = evmeta[7]
        msms.category[mask[i], :] = evmeta[8]
        climout[mask[i], :] = eachrow(msst.clima)[i]
        threshout[mask[i], :] = eachrow(msst.thresh)[i]
        meanmets, evmets = meanmetrics(evmeta[1], evmeta[2], evmeta[3], fullyears)
        rons[i], rdcs[i], sdt[i], edt[i] = evmeta[2], evmeta[3], evmeta[5], evmeta[6]
        rws[i] = repeat([mask[i][1]], length(evmeta[2]))
        cls[i] = repeat([mask[i][2]], length(evmeta[2]))
        mni[i], cmi[i], mxi[i], drt[i], vri[i] = evmets
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

    edf = (MeanInt=reduce(vcat, mni), CumInt=reduce(vcat, cmi), MaxInt=reduce(vcat, mxi), Duration=reduce(vcat, drt), Category=reduce(vcat, cats), ROnset=reduce(vcat, rons), RDecline=reduce(vcat, rdcs), VarInt=reduce(vcat, vri), StartDate=reduce(vcat, sdt), EndDate=reduce(vcat, edt), Xind=reduce(vcat, rws), Yind=reduce(vcat, cls))

    return msms, edf, climout, threshout
end


function testarrays(sst, sstdate)
    errmsgdims = "Dimensions of temperature array and time do not match: "
    size(sst, ndims(sst)) == size(sstdate, ndims(sstdate)) || throw(DimensionMismatch("$(errmsgdims) temperature: $(size(sst, ndims(sst))), time: $(size(sstdate, ndims(sstdate)))."))
end

function testdates(sstdate, mhwdate, climdate)
    errmsglen = "The end should be greater than the start: "
    errmsgform = "The date should be either a string as in 'yyyy-mm-dd' -> '1990-01-02' or a date as in `Date(yyyy,mm,dd)` -> `Date(1990,3,12)`."
    errmsgrange = "TimeError: Period is outside temperature date period: "
    # errmsglen2 = "The length of period is greater than the length of `sstdate`. Check: "

    # Format/type tests
    (first(sstdate) isa String || first(sstdate) isa Date) || errmsgform

    (first(climdate) isa String || first(climdate) isa Date) || errmsgform
    (first(mhwdate) isa String || first(mhwdate) isa Date) || errmsgform

    # Assure dates are date objects
    ymd = "yyyy-mm-dd"
    sstdate = eltype(sstdate) == String ? range(Date(first(sstdate), ymd), Date(last(sstdate), ymd)) : sstdate
    mhwdate = eltype(mhwdate) == String ? range(Date(first(mhwdate), ymd), Date(last(mhwdate), ymd)) : mhwdate
    climdate = eltype(climdate) == String ? range(Date(first(climdate), ymd), Date(last(climdate), ymd)) : climdate

    # Make sure the order of the dates is correct
    first(sstdate) < last(sstdate) || errmsglen
    first(mhwdate) < last(mhwdate) || errmsglen
    first(climdate) < last(climdate) || errmsglen

    # Range tests
    first(mhwdate) ≥ first(sstdate) && last(mhwdate) ≤ last(sstdate) || error(errmsgrange, mhwdate)

    first(climdate) ≥ first(sstdate) && last(climdate) ≤ last(sstdate) || error(errmsgrange, climdate)

    climdate = range(max(climdate[1] - Day(5), first(sstdate)), min(climdate[end] + Day(5), last(sstdate)))
    return sstdate, mhwdate, climdate
end

function mhctemp(sst, sstdate::StepRange{Date,Day}, mdate::StepRange{Date,Day}, cdate::StepRange{Date,Day}; threshold=threshold)
    sstdate, mdate, cdate = testdates(sstdate, mdate, cdate)
    testarrays(sst, sstdate)
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
    evanom, rons, rdcs, cats, stdate, endate, mhwout, catout = eventmetrics(mhwin, msxs, mexs)
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
    edf = (MeanInt=evmets[1], CumInt=evmets[2], MaxInt=evmets[3], Duration=evmets[4], Category=cats, ROnset=rons, RDecline=rdcs, VarInt=evmets[5], StartDate=stdate, EndDate=endate)
    return msms, edf, mhwin.clima, mhwin.thresh
end

function edetect(sst::Array, sstdate, mdate, cdate; threshold=0.9)
    excdfn = >=
    mhwin = MHTemp(mhctemp(sst, sstdate, mdate, cdate; threshold)..., excdfn, argmax, maximum)
    msms = vecarr(sst, mdate)
    mexc = exceed(mhwin)
    msxs, mexs = _mylabeling(mexc)
    msms, edf, cl, th = matevents(msms, mhwin, mexc, msxs, mexs)
    return msms, edf, cl, th
end

function vecarr(sst, mhwdate)
    N = ndims(sst)
    N ∈ (1, 3) || throw("Oh heat! The dimension of the `sst` should be 1 or 3, we got $(N) instead!")
    T = eltype(sst)
    TN = Array{Union{T,Missing},N}
    TM = Matrix{Union{T,Missing}}
    TB = Array{Union{Bool,Missing},N}
    NT = NamedTuple

    trds = (:means, :trends, :pvalues, :pmetrics)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maxint, :days, :frequency)
    lt, lm = length(trds), length(metrics)
    x, y = ifelse(N == 3, size(sst), (1, 1))
    z = length(mhwdate)
    zz = length(unique(year.(mhwdate)))
    ds = N == 1 ? z : (x, y)

    mhwexd = N == 1 ? TB(missing, ds) : TB(missing, ds..., z)
    mhwtemp, mhwcat = N == 1 ? ntuple(_ -> TN(missing, ds), 2) : ntuple(_ -> TN(missing, ds..., z), 2)
    annuals = N == 1 ? NT{metrics}(ntuple(_ -> TN(missing, zz), lm)) : NT{metrics}(ntuple(_ -> TN(missing, ds..., zz), lm))
    mets = N == 1 ? ntuple(_ -> NT{metrics}(ntuple(_ -> TN(missing, N), lm)), lt) : ntuple(_ -> NT{metrics}(ntuple(_ -> TM(missing, ds...), lm)), lt)
    return MarineHW(mhwtemp, mhwcat, mhwexd, annuals, mets...)
end

function MarineHW(sst::Array, sdate, mdate, cdate; threshold=0.9)
    excdfn, argfn, anomfn = ≥, argmax, maximum
    mhwin = MHTemp(mhctemp(sst, sdate, mdate, cdate; threshold)..., excdfn, argfn, anomfn)
    msms = vecarr(sst, mdate)
    msms, edf, cl, th = edetect(mhwin, msms)
    return msms, edf, cl, th
end

function MarineCS(sst::Array, sdate, mdate, cdate; threshold=0.1)
    excdfn, argfn, anomfn = ≤, argmin, minimum
    mhwin = MCTemp(mhctemp(sst, sdate, mdate, cdate; threshold)..., excdfn, argfn, anomfn)
    msms = vecarr(sst, mdate)
    msms, edf, cl, th = edetect(mhwin, msms)
    return msms, edf, cl, th
end

function edetect(mhwin::MarineHeatWave, msms::MarineHW)
    mexc = exceed(mhwin)
    msxs, mexs = _mylabeling(mexc)
    msms, edf, cl, th = matevents(msms, mhwin, mexc, msxs, mexs)
    return msms, edf, cl, th
end

function matevents(msms::MarineHW, msst::MarineHeatWave, mexc::BitVector, mstartsxs, mendsxs)
    evanom, rons, rdcs, cats, stdate, endate, mhwout, catout = eventmetrics(msst, mstartsxs, mendsxs)
    fullyears = unique(year.(msst.dates))
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
    edf = (MeanInt=evmets[1], CumInt=evmets[2], MaxInt=evmets[3], Duration=evmets[4], Category=cats, ROnset=rons, RDecline=rdcs, VarInt=evmets[5], StartDate=stdate, EndDate=endate)
    return msms, edf, msst.clima, msst.thresh
end

evtable(mhw) = DataFrame(mhw[2])
evclim(mhw) = mhw[3]
evthresh(mhw) = mhw[4]

