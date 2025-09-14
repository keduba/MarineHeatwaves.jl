
import Statistics : mean, quantile, std
import Distributions: cdf, FDist, TDist
using Dates
using NCDatasets
using SparseArrays

const T = Float32
const TI = Int16

struct MHWrapper{T} 
    anom::Array{T}
    threshanom::Array{T}
    onset::T
    decline::T
end

struct MCWrapper{T}
    anom::Array{T}
    threshanom::Array{T}
    onset::T
    decline::T
end

struct EventsHW{T}
    means::Array{T, 1}
    minimaxes::Array{T, 1}
    onset::Array{T, 1}
    decline::Array{T, 1}
    duration::Array{T, 1}
    sums::Array{T, 1}
end

struct EventsCS{T}
    means::Array{T, 1}
    minimaxes::Array{T, 1}
    onset::Array{T, 1}
    decline::Array{T, 1}
    duration::Array{T, 1}
    sums::Array{T, 1}
end

struct MCS{T<:AbstractVecOrMat{<:AbstractFloat}}
    temp::T
    clim::T
    thresh::T
    exceeds::Array{Bool}
end

function MHW(temp, clim, thresh)
    exc = _excess(temp, thresh)
    MCS{typeof(temp), typeof(clim), typeof(thresh), typeof(exc)}(temp, clim, thresh, exc)
end

function MCS(temp, clim, thresh)
    exc = _excess(thresh, temp)
    MCS{typeof(temp), typeof(clim), typeof(thresh), typeof(exc)}(temp, clim, thresh, exc)
end

function MCS(temp, clim, thresh, event=:mhw)
    exc = event == :mhw ? _excess(temp, thresh) : _excess(thresh, temp)
    MCS{typeof(temp), typeof(clim), typeof(thresh), typeof(exc)}(temp, clim, thresh, exc)
end

struct MHWCSO{T<:AbstractVecOrMat{<:AbstractFloat}}
    outtemp::T
    outclim::T
    outhresh::T
    outexceeds::Array{Bool}
    outmeans::T
    outannuals::T
    outpvalues::T
    outcoeff::T
    outrsquared::T
end

for op = (:length,  :maximum, :minimum, :argmax, :argmin, :sum,  :first, :last)
    @eval begin
        Base.$op(a::MHWrapper) = $op(a.anom)
        Base.$op(a::MCWrapper) = $op(a.anom)
    end
end

for op = (:mean, :std)
    @eval begin
        $op(a::MHWrapper) = $op(a.anom)
        $op(a::MCWrapper) = $op(a.anom)
    end
end

mhcsminimax(a::MHWrapper) = maximum(a)
mhcsminimax(a::MCWrapper) = minimum(a)

mhcsarg(a::MHWrapper) = argmax(a)
mhcsarg(a::MCWrapper) = argmin(a)

# we expect Events to currently return Vector{Vector{T}} so the following function is targeting the field we want.

function mhcsminimax(x::EventsHW)
    return @inbounds [maximum(ia) for ia in x.minimaxes]
end

function mhcsminimax(x::EventsCS)
    return @inbounds [minimum(ia) for ia in x.minimaxes]
end

function timeindices(sstdate::Date, evtdate::Date)
    si::TI = findfirst(isequal(first(evtdate)), sstdate)
    ei::TI = findfirst(isequal(last(evtdate)), sstdate)
    return range(si, ei)
end

seamask(sst::AbstractArray, dims) = dropdims(count(!isnan, sst, dims=dims) .> 0, dims=dims)

function _subtemp(sst::AbstractArray, mhwix::Range, clmix::Range)
    N = ndims(sst)
    N ≤ 0 && error("0-dimensional data just can't work.")

    xz = seamask(sst, N)
    !xz[1] && error("Please check your data. It appears to be all `NaN` or `missing`.")

    if N > 1
        CIx = CartesianIndices(xz)[xz]
        @views begin
            msst = permutedims(sst[CIx, mhwix])
            csst = permutedims(sst[CIx, clmix])
        end
        x, y = Base.front(size(sst))
        nCIx = setdiff(CartesianIndices((x, y)), CIx)
        return msst, csst, (CIx, nCIx, x, y)
    elseif N == 1
        CIx = TI(xz[1])
        @views begin
            msst = sst[mhwix]
            csst = sst[clmix]
        end
        return msst, csst, CIx
    end
end


