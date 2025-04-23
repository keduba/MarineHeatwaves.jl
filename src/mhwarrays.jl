# TODO: 
# rewrite all the variables to have a consistent name for tracking across the file.
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

# Defaults
win_width = 5    # for the climatology
pct_width = 31   # for smoothing data 
min_dur = 5      # minimum duration of marine heatwave event
max_gap = 2      # maximum gap between successive marine heatwave event

Base.length(fl::Base.Iterators.Flatten{<:AbstractVector{<:UnitRange}}) = sum(length, fl.it)

Base.IteratorSize(::Base.Iterators.Flatten{<:AbstractVector{<:UnitRange}}) = Base.HasLength()

"""
    seamask(sst) -> BitVecOrMat

Returns only sea areas. Optional dimension. For arrays, expected dimension is the 3rd dimension.
"""
seamask(sst, dims=3) = dropdims(count(!isnan, sst, dims=dims) .> 0, dims=dims)

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

"""
    subtemp(sst, sstdate, eventdate) -> eventsst, mask, eventleapyearday

`subtemp` uses the `timeindices` and `_subtemp` to return array, landmask and leapyearday.
"""
function subtemp(sst, sstdate, evdate)
    mhwix = timeindices(sstdate, evdate)
    emsst, mask = _subtemp(sst, mhwix)
    elyd = leapyearday.(sstdate[mhwix])
    return emsst, mask, elyd
end

"""
helper function to select subset a vector `av` with indices from a range `dr`.
"""
function findices(av, dr)
    (av[i] for i in flatten(dr))
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

