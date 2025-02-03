# TODO 
# rewrite all the variables to have a consistent name for tracking across the file.
#
#-

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
    mask = CartesianIndices(xz)[xz]
    msst = sst[mask, mhwix]'
    return msst, mask
end

function subtemp(sst, sstdate, evdate)
    mhwix, elyd = timeindices(sstdate, evdate)
    emsst, mask = _subtemp(sst, mhwix)
    return emsst, mask, elyd
end

function findices(av, dr)
    (av[i] for i in flatten(dr))
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

evtable(mhw) = DataFrame(mhw[2])
evclim(mhw) = mhw[3]
evthresh(mhw) = mhw[4]

