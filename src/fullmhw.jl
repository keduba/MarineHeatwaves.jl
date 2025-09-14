
import Statistics : mean, quantile, std
import Distributions: cdf, FDist, TDist
using Dates
using NCDatasets
using SparseArrays

const T = Float32
const TI = Int16

abstract type MExtreme end
abstract type MWrapper end

struct MHWrapper{T} <: MWrapper 
    anom::Array{T}
    threshanom::Array{T}
    onset::T
    decline::T
end

struct MCWrapper{T} <: MWrapper 
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

struct MCS{T<:AbstractVecOrMat{<:AbstractFloat}} <: MExtreme
    temp::T
    clim::T
    thresh::T
    exceeds::Array{Bool}
end

struct MHW{T<:AbstractVecOrMat{<:AbstractFloat}} <: MExtreme
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

excess(ms::MExtreme) = ms.exceeds

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

function urange(clyd::Vector{TI}, win::TI)
    out = [[] for _ in 1:366]
    for n in 1:366
           for (x, ) in Iterators.filter(p -> isequal(n, p.second), pairs(clyd))
               push!(out[n], max(1, x-win):min(length(clyd), x+win))
           end
    end
    out
end

function _mylabel(ms::MExtreme, mindur::TI, maxgap::TI)
    sty = excess(ms)
    stb = sparse(diff(sty, dims=1))
    cstt = Vector{Vector{TI}}(undef,  size(sty, 2))
    csee = Vector{Vector{TI}}(undef,  size(sty, 2))
    cols = TI[]
    for c in axes(stb, 2)
        cst = TI[] 
        cse = TI[]
        for i in nzrange(stb, c)
            if isequal(nonzeros(stb)[i], -1)
                push!(cse, rowvals(stb)[i])
            else
                push!(cst, rowvals(stb)[i])#+1)
            end
        end
        first(sty[:, c]) ? pushfirst!(cst, 1) : nothing
        last(sty[:, c]) ? push!(cse, lastindex(sty, 1)) : nothing
        length(cse) == length(cst)
        dur = cse - cst
        stss = cst[dur .≥ mindur]
        enss = cse[dur .≥ mindur]
        dft = stss[2:end] - enss[1:end-1] .> maxgap
        isempty(dft) && continue
        keepat!(stss, [true; dft])
        keepat!(enss, [dft; true])
        cstt[c] = stss
        csee[c] = enss
        push!(cols, c)
    end
    cstt, csee, cols
end

function _anomsa(M::DataType, sst::Vector{T}, clim::Vector{T}, thsh::Vector{T}, st::TI, se::TI, lm::TI)
    @views begin
    ant = sst[st:se] - clim[st:se]
    tht = thsh[st:se] - clim[st:se]
    ont = sst[max(1, st - 1)] - clim[max(1, st - 1)]
    dnt = sst[min(lm, se + 1)] - clim[min(lm, se + 1)]
    end
    return M(ant, tht, ont, dnt)
end

function _onset2(atod::MWrapper, mst)
    fan = first(atod) 
    nmx = mhcsminimax(atod) 
    ngx = mhcsarg(atod) 
    lnmx = nmx - 0.5(fan + atod.onset)
    snmx = nmx - fan
    mst > 1 ? /(lnmx, (ngx + 0.5)) : /(snmx, ngx)
end

function _decline2(atod::MWrapper, mse, lnx) 
    lan = last(atod)
    nmx = mhcsminimax(atod) 
    ngx = mhcsarg(atod) 
    wsnx::T = length(atod) - ngx
    lnmx = nmx - 0.5(lan + atod.decline)
    snmx = nmx - lan
    mse < lnx ? /(lnmx, (wsnx + 0.5)) : ngx == lnx ? snmx : /(snmx, wsnx)
end

function anomsa(m::MExtreme, evst, indices)
    MW, Ev = typeof(m) == MHW{Matrix{T}} ? (MHWrapper, EventsHW) : (MCWrapper, EventsCS)
    CIx, nCIx, x, y = indices
    mst, mse, cols = evst
    lm::TI = size(m.temp, 1)
    mt::TI = 6
    outemp, outhsh, outcat = ntuple(_ -> Array{T, 3}(undef, x, y, lm), 3)
    onsan, decan, means, cums, maxes, durs = ntuple(_ -> [Vector{T}(undef,m) for m in length.(mst)], mt)
    for (c, cst, cse) in zip(cols, mst, mse)
        for (d, (st, se)) in enumerate(zip(cst, cse))
            atod = _anomsa2(MW, m.temp[:, c], m.clim[:, c], m.thresh[:, c], st, se, lm)
            outemp[CIx[c], st:se] = atod.anom
            outhsh[CIx[c], st:se] = atod.threshanom
            outcat[CIx[c], st:se] .= categorys(atod)
            onsan[c][d] = _onset2(atod, st)
            decan[c][d] = _decline2(atod, se, lm)
            means[c][d], cums[c][d], maxes[c][d], durs[c][d] = _events2(atod)# vars[c][d]
        end
    end
    for bl in (outemp, outhsh, outcat)
        bl[nCIx, :] .= NaN
    end
    # if using mutable struct, we can also wrap outemp, outhsh, outcat
    return outemp, outhsh, outcat, Ev(onsan, decan, means, cums, maxes, durs)
end

_categorys(anom::Vector{T}, thsd::Vector{T}) = min(4, maximum(fld.(anom, thsd))) 

categorys(a::MWrapper) =  _categorys(a.anom, a.threshanom)

function _events(anom::MWrapper)
    # Per Event Metrics
    means = mean(anom)
    cums = sum(anom)
    maxes = mhcsminimax(anom)
    durs = length(anom)
    return means, cums, maxes, durs
end

