
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

struct MHW{T<:AbstractVecOrMat{<:AbstractFloat}}
    temp::T
    clim::T
    thresh::T
    exceeds::Array{Bool}
end

function MHW(temp, clim, thresh)
    exc = _excess(temp, thresh)
    MHW{typeof(temp), typeof(clim), typeof(thresh), typeof(exc)}(temp, clim, thresh, exc)
end

function MCS(temp, clim, thresh)
    exc = _excess(thresh, temp)
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

MarineHW(args...)::MHW = MHW(subtemp(args...)...)
MarineCS(args...)::MCS = MCS(subtemp(args...)...)

leapyearday(mts::Date)::TI = dayofyear(mts) > 59 && !isleapyear(mts) ? dayofyear(mts) + 1 : dayofyear(mts)

_excess(tp, th) = tp .> th

function subtemp(sst, sstdate::Date, mhwdate::Date, clmdate::Date)
    mhwix = timeindices(sstdate, mhwdate)
    clmix = timeindices(sstdate, clmdate)
    mlyd = leapyearday.(sstdate[mhwix])
    clyd = leapyearday.(sstdate[clmix])
    mhtemp, ctemp, CIs = _subtemp(sst, mhwix, clmix)
    clima, thresh = climthresh(cltemp, clyd, mlyd, window, smoothwindow, threshold)

    return mhtemp, clima, thresh
end

function climthresh(cmst, clyd, mlyd, win::TI, swindow::TI, threshold::T)
    # T =  eltype(cmst)
    clim = Matrix{T}(undef, length(mlyd), size(cmst, 2))
    thresh = similar(clim)
    uranges = urange(clyd, win)
    cv = Vector{T}(undef, 366)
    th = similar(cv)
    for (n, vec) in enumerate(eachcol(cmst))
        for (m, ur) in enumerate(uranges)
            cvc = Iterators.flatten([vec[i] for i in ur])
            cv[m] = mean(cvc)
            th[m] = quantile(cvc, threshold)
        end
        cv[60] = 0.5(cv[59] + cv[61])
        th[60] = 0.5(th[59] + th[61])
        clim[:, n] = moving_means(cv, swindow, mlyd)
        thresh[:, n] = moving_means(th, swindow, mlyd)
    end
    clim, thresh
end

function moving_means(mt::Vector, pwidth, lyd; wrap=true)
    out = Vector{eltype(mt)}(undef, length(lyd))
    ins = moving_mean(mt, pwidth, wrap)
    for (i, t) in enumerate(lyd)
        @views out[i] = ins[t]
    end
    out
end

function moving_mean(A::AbstractVector, m::TI, wrap::Bool) 
    out = similar(A)
    R = LinearIndices(A)
    m2 = trunc(TI, 0.5m)
    Ifirst, Ilast = first(R), last(R)
    I1 = m2 * oneunit(Ifirst)
    if wrap
        R = axes(A, 1) .+ m2
        A = vcat(last(A, m2), A, first(A, m2))
        Ilast = lastindex(A)
    end
    for (I, J) in pairs(R)
        n, s = 0, zero(eltype(out))
        for K in max(Ifirst, J - I1):min(Ilast, J + I1)
            s += A[K]
            n += 1
        end
        out[I] = s / n
    end
    return out
end

