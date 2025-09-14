
import Statistics : mean, quantile, std
import Distributions: cdf, FDist, TDist
using Dates
using NCDatasets
using SparseArrays

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
        Base.$op(a::EventsHW) = $op(a.minimaxes)
        Base.$op(a::EventsCS) = $op(a.minimaxes)
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

mhcsminimax(a::EventHW) = maximum(a)
mhcsminimax(a::EventCS) = minimum(a)

mhcsarg(a::MHWrapper) = argmax(a)
mhcsarg(a::MCWrapper) = argmin(a)



