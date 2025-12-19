import Statistics: mean, quantile, std
import Distributions: cdf, FDist, TDist
using Dates
using NCDatasets
using SparseArrays

const T = Float32
const TI = Int16

const metrics = Dict(:means => 1,
                     :maxes => 2,
                     :sums => 3,
                     :duration => 4,
                     :frequency => 5,
                     :days => 6,
                     :onset => 7,
                     :decline => 8,
                     :categorys => 9)



abstract type MWrapper end
abstract type MEvents{TD} end 
abstract type MExtreme{TD<:AbstractVecOrMat{<:AbstractFloat}, V<:BitArray} end 


struct MCWrapper{T} <: MWrapper
    anom::Vector{T}
    threshanom::Vector{T}
    onset::T
    decline::T
end

struct MHWrapper{T} <: MWrapper
    anom::Vector{T}
    threshanom::Vector{T}
    onset::T
    decline::T
end
#
# struct Events{TE} <: MEvents{TE}
#     means::TE
#     minimaxes::TE
#     onset::TE
#     decline::TE
#     duration::TE
#     sums::TE
#     categorys::TE
#     dtype::Type{<:MEvents}
# end

struct EventsFull{TE<:AbstractVector, Ti<:AbstractFloat, N} #<: MEvents{TE}
    means::TE
    minimaxes::TE
    onset::TE
    decline::TE
    duration::TE
    sums::TE
    categorys::TE
    tpanom::Array{Ti, N}
    thanom::Array{Ti, N}
    categoryarr::Array{Ti, N}
    dtype::Type{<:MEvents}
end

struct EventHW{TA} <: MEvents{TA}
    maxes::TA
end

struct EventCS{TA} <: MEvents{TA}
    maxes::TA
end

struct MCS{TA,V} <: MExtreme{TA,V}
    temp::TA
    clim::TA
    thresh::TA
    exceeds::V
    edtype::Type{EventCS}
end

struct MHW{TA, V} <: MExtreme{TA,V}
    temp::TA
    clim::TA
    thresh::TA
    exceeds::V
    edtype::Type{EventHW}
end

function MHW(temp, clim, thresh)
    temp = convert(typeof(clim), temp)
    exc = _excess(temp, thresh)
    MHW(temp, clim, thresh, exc, EventHW)
end

function MCS(temp, clim, thresh)
    temp = convert(typeof(clim), temp)
    exc = _excess(thresh, temp)
    MCS{typeof(temp), typeof(clim), typeof(thresh), typeof(exc)}(temp, clim, thresh, exc, EventCS)
end

struct MHCMetrics{T<:AbstractFloat, Nm, N, TD<:AbstractArray{T, N}} 
    annuals::NTuple{Nm, TD}
    means::TD
    coeffs::TD
    errors::TD
    rsquared::TD
    intercepts::TD
    pvalues::TD
end

struct MHWCSO{T<:AbstractArray{<:AbstractFloat}}
    # outtempanom::T
    # outcats::T
    # outthreshanom::T
    # outclim::T
    # outhresh::T
    # outexceeds::Array{Bool}
    outmeans::T
    outannuals::T
    outpvalues::T
    outcoeff::T
    outrsquared::T
end

for op = (:length, :maximum, :minimum, :argmax, :argmin, :sum,  :first, :last)
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

mhcsminimax(x::EventHW{Vector{T}}) = maximum(x.maxes)
mhcsminimax(x::EventCS{Vector{T}}) = minimum(x.maxes)

function mhcsminimax(x::EventHW{Vector{Vector{T}}})
    return @inbounds [maximum(ia) for ia in x.maxes]
end

function mhcsminimax(x::EventCS{Vector{Vector{T}}})
    return @inbounds [minimum(ia) for ia in x.maxes]
end

"Return the range of date B in date A."
function timeindices(sstdate::StepRange{Date, Day}, evtdate::StepRange{Date, Day})
    si::TI = findfirst(isequal(first(evtdate)), sstdate)
    ei::TI = findfirst(isequal(last(evtdate)), sstdate)
    return range(si, ei)
end

seamask(sst::AbstractArray, dims) = dropdims(count(!isnan, sst, dims=dims) .> 0, dims=dims)

"""
    Calculate the subset of the input array. For n-dimensional array input where n > 1, return n-1 dimension. Return the subset and indices of mask.

"""
function _subtemp(sst::AbstractArray, mhwix::UnitRange, clmix::UnitRange)
    N = ndims(sst)
    N ≤ 0 && error("0-dimensional data just can't work.")

    xz = seamask(sst, N)
    # TODO: condition on N==1 
    if N == 1
        !xz[1] && error("Please check your data. It appears to be all `NaN` or `missing`.")
    end
    N == 1 && xz[1] || error("Please check your data. It appears to be all `NaN` or `missing`.")

    if N == 3
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


leapyearday(mts::Date)::TI = dayofyear(mts) > 59 && !isleapyear(mts) ? dayofyear(mts) + 1 : dayofyear(mts)

_excess(tp, th) = tp .> th

excess(ms::MExtreme) = ms.exceeds

const events = Dict(:mhw => (MHW, 0.9),
                    :mcs => (MCS, 0.1))

"""
    Return a MExtreme (MHW or MCS).

"""
function mextreme(sst, sstdate::StepRange, mhwdate::StepRange, clmdate::StepRange, event=:mhw; window=5, smoothwindow=31, threshold=nothing)
    window::TI = convert(TI, window)
    smoothwindow = convert(TI, smoothwindow)
    # convert(typeof)
    in(event, keys(events)) || error(event," is not a valid event. Try `:mhw` or `:mcs`")
    ME, mthreshold = get(events, event, :mhw)
    threshold::T = isnothing(threshold) ? convert(T, mthreshold) : convert(T, threshold)
    mhwix = timeindices(sstdate, mhwdate)
    clmix = timeindices(sstdate, clmdate)
    mlyd = leapyearday.(sstdate[mhwix])
    clyd = leapyearday.(sstdate[clmix])
    mhtemp, ctemp, CIs = _subtemp(sst, mhwix, clmix)
    clima, thresh = climthresh(ctemp, clyd, mlyd, window, smoothwindow, threshold)
    return ME(mhtemp, clima, thresh), CIs
end

"""
    Calculate and return the climatology mean and quantile.

"""
function climthresh(cmst::Matrix{T}, clyd, mlyd, win::TI, swindow::TI, threshold::T; wrap::Bool=true)
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
        clim[:, n] = moving_means(cv, swindow, mlyd, wrap=wrap)
        thresh[:, n] = moving_means(th, swindow, mlyd, wrap=wrap)
    end
    clim, thresh
end

function climthresh(cmst::Vector{T}, clyd, mlyd, win::TI, swindow::TI, threshold::T; wrap::Bool=true)
    uranges = urange(clyd, win)
    cv = Vector{T}(undef, 366)
    th = similar(cv)
    for (m, ur) in enumerate(uranges)
        cvc = Iterators.flatten([cmst[i] for i in ur])
        cv[m] = mean(cvc)
        th[m] = quantile(cvc, threshold)
    end
    cv[60] = 0.5(cv[59] + cv[61])
    th[60] = 0.5(th[59] + th[61])
    clim = moving_means(cv, swindow, mlyd, wrap=wrap)
    thresh = moving_means(th, swindow, mlyd, wrap=wrap)
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

"Return a vector of all indices of 366 day from the leapyearday vector"
function urange(clyd::Vector{TI}, win::TI)
    out = [[] for _ in 1:366]
    for n in 1:366
       for (x, ) in Iterators.filter(p -> isequal(n, p.second), pairs(clyd))
           push!(out[n], max(1, x-win):min(length(clyd), x+win))
       end
    end
    out
end

# TODO: remove the type signature from the mindur and maxgap. Enforce it inside if necessary.
# umbrella function should call mylabel and the anomsa in one step.
# Exceedance can be stored as a sparse array.

"""
    Return the start and end positions  and pixel columns where MExtremes were detected.

"""
function _mylabel(ms::MExtreme{Matrix{T}}, mindur::TI, maxgap::TI)
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
                push!(cst, rowvals(stb)[i]+1)
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

function _mylabel(ms::MExtreme{Vector{T}}, mindur::TI, maxgap::TI) 
    sty = excess(ms)
    stb = sparse(diff(sty))
    cst = TI[] 
    cse = TI[]
    for i in nzrange(stb, c)
        if isequal(nonzeros(stb)[i], -1)
            push!(cse, rowvals(stb)[i])
        else
            push!(cst, rowvals(stb)[i]+1)
        end
    end
    first(sty) ? pushfirst!(cst, 1) : nothing
    last(sty) ? push!(cse, lastindex(sty, 1)) : nothing
    length(cse) == length(cst)
    dur = cse - cst 
    stss = cst[dur .≥ mindur]
    enss = cse[dur .≥ mindur]
    dft = stss[2:end] - enss[1:end-1] .> maxgap
    isempty(dft) && error("No event was detected")
    keepat!(stss, [true; dft])
    keepat!(enss, [dft; true])
    stss, enss, 1
end

function _anomsa(sst::Vector{T}, clim::Vector{T}, thsh::Vector{T}, st::TI, se::TI, lm::TI)
    @views begin
    ant = sst[st:se] - clim[st:se]
    tht = thsh[st:se] - clim[st:se]
    ont = sst[max(1, st - 1)] - clim[max(1, st - 1)]
    dnt = sst[min(lm, se + 1)] - clim[min(lm, se + 1)]
    end
    return ant, tht, ont, dnt
end

function anomsa(m::MHW{Matrix{T}}, c, st, en, ls) 
    outs = _anomsa(m.temp[:,c], m.clim[:,c], m.thresh[:,c], st, en, ls)
    return MHWrapper(outs...)
end

function anomsa(m::MCS{Matrix{T}}, c, st, en, ls) 
    outs = _anomsa(m.temp[:,c], m.clim[:,c], m.thresh[:,c], st, en, ls)
    return MCWrapper(outs...)
end

function anomsa(m::MCS{Vector{T}}, st, en, ls) 
    outs = _anomsa(m.temp, m.clim, m.thresh, st, en, ls)
    return MCWrapper(outs...)
end

function anomsa(m::MHW{Vector{T}}, st, en, ls) 
    outs = _anomsa(m.temp, m.clim, m.thresh, st, en, ls)
    return MHWrapper(outs...)
end

"""
    Return the rate of onset for each event.

"""
function _onset(atod::MWrapper, mst)
    fan = first(atod) 
    nmx = mhcsminimax(atod) 
    ngx = mhcsarg(atod) 
    lnmx = nmx - 0.5(fan + atod.onset)
    snmx = nmx - fan
    mst > 1 ? /(lnmx, (ngx + 0.5)) : /(snmx, ngx)
end

"""
    Return the rate of decline for each event.

"""
function _decline(atod::MWrapper, mse, lnx)
    lan = last(atod)
    nmx = mhcsminimax(atod) 
    ngx = mhcsarg(atod) 
    wsnx::T = length(atod) - ngx
    lnmx = nmx - 0.5(lan + atod.decline)
    snmx = nmx - lan
    mse < lnx ? /(lnmx, (wsnx + 0.5)) : ngx == lnx ? snmx : /(snmx, wsnx)
end

"""
    Compute the anomaly of the events in each pixel and return the Event.

"""
function anomsa(m::MExtreme{Matrix{T}}, evst, indices) 
    CIx, nCIx, x, y = indices
    mst, mse, cols = evst
    lm::TI = size(m.temp, 1)
    mt::TI = 7
    outemp, outhsh, outcat = ntuple(_ -> Array{T, 3}(undef, x, y, lm), 3)
    onsan, decan, means, cums, maxes, durs, catso = ntuple(_ -> [Vector{T}(undef,m) for m in length.(mst)], mt)
    for (c, cst, cse) in zip(cols, mst, mse)
        for (d, (st, se)) in enumerate(zip(cst, cse))
            atod = anomsa(m, c, st, se, lm)
            outemp[CIx[c], st:se] = atod.anom
            outhsh[CIx[c], st:se] = atod.threshanom
            cats = categorys(atod)
            outcat[CIx[c], st:se] .= cats
            catso[c][d] = cats
            onsan[c][d] = _onset(atod, st)
            decan[c][d] = _decline(atod, se, lm)
            means[c][d], cums[c][d], maxes[c][d], durs[c][d] = _events(atod)
        end
    end
    for bl in (outemp, outhsh, outcat)
        bl[nCIx, :] .= NaN
    end
    # return outemp, outhsh, outcat, Events(means, maxes, onsan, decan, durs, cums, catso, m.edtype)
    return EventsFull(means, maxes, onsan, decan, durs, cums, catso, outemp, outhsh, outcat, m.edtype)
end

function anomsa(m::MExtreme{Vector{T}}, evst, indices)
    mst, mse = evst
    lm::TI = size(m.temp, 1)
    mt::TI = 7
    outemp, outhsh, outcat = ntuple(_ -> Vector{T}(undef, lm), 3)
    onsan, decan, means, cums, maxes, durs, catso = ntuple(_ -> Vector{T}(undef,  length(mst)), mt)
    for (d, (st, se)) in enumerate(zip(mst, mse))
        atod = _anomsa(m.wdtype, m.temp, m.clim, m.thresh, st, se, lm)
        outemp[st:se] = atod.anom
        outhsh[st:se] = atod.threshanom
            cats = categorys(atod)
            outcat[CIx[c], st:se] .= cats
            catso[c][d] = cats
        onsan[c][d] = _onset(atod, st)
        decan[c][d] = _decline(atod, se, lm)
        means[c][d], cums[c][d], maxes[c][d], durs[c][d] = _events(atod)
    end
    # if using mutable struct, we can also wrap outemp, outhsh, outcat
    # return outemp, outhsh, outcat, Events(means, maxes, onsan, decan, durs, cums, catso, m.edtype)
    return EventsFull(means, maxes, onsan, decan, durs, cums, catso, outemp, outhsh, outcat, m.edtype)
end

_categorys(anom::Vector{T}, thsd::Vector{T}) = min(4, maximum(fld.(anom, thsd))) 

categorys(a::MWrapper) = _categorys(a.anom, a.threshanom)

function _events(anom::MWrapper)
    # Per Event Metrics
    means = mean(anom)
    cums = sum(anom)
    maxes = mhcsminimax(anom)
    durs = length(anom)
    return means, cums, maxes, durs
end

"""
    Return the mean of the metrics in each pixel.

"""
function meanmetrics(ev::MEvents{Vector{Vector{T}}}, indices, mdate)
    CIx, nCIx, x, y = indices
    lfy = (length ∘ unique)(year.(mdate))
    z = 8 #length(metrics)
    outmean = ntuple(_ -> Matrix{T}(undef, x, y), z) 
    for i in eachindex(outmean)
        outmean[i][nCIx] .= T(NaN)
    end
    # check the metrics dictionary to ensure order of variables
    for fm in (:means, :sums, :onset, :decline, :duration)
        idx = metrics[fm]
        outmean[idx][CIx]  = mean.(getfield(ev, fm))
    end
    E = ev.dtype
    outmean[metrics[:maxes]][CIx] = mhcsminimax(E(ev.minimaxes)) 
    outmean[metrics[:frequency]][CIx] = @. length(ev.duration) / lfy
    outmean[metrics[:days]][CIx] = @. sum(ev.duration) / lfy 
    outmean
end

# potentially
"""
    Return the mean of the metrics in each pixel.
"""
function meanmetrics2(ev::MEvents{Vector{Vector{T}}}, indices, mdate)
    CIx, nCIx, x, y = indices
    lfy = (length ∘ unique)(year.(mdate))
    z = 8# length(metrics)
    outmean = Array{T, 3}(undef, x, y, z)
    outmean[nCIx, :] .= T(NaN)
    # check the metrics dictionary to ensure order of variables
    for fm in (:means, :sums, :onset, :decline, :duration)
        idx = metrics[fm]
        outmean[CIx, idx]  = mean.(getfield(ev, fm))
    end
    E = ev.dtype
    outmean[CIx, metrics[:maxes]] = mhcsminimax(E(ev.minimaxes)) 
    outmean[CIx, metrics[:frequency]] = @. length(ev.duration) / lfy
    outmean[CIx, metrics[:days]] = @. sum(ev.duration) / lfy 
    outmean
end

function meanmetrics(ev::MEvents{Vector{T}}, mdate)
    # Vector version
    lfy = (length ∘ unique)(year.(mdate))
    z = 8 #length(metrics)
    outmean = Vector{T}(undef, z) 
    # check the metrics dictionary to ensure order of variables
    for fm in (:means, :sums, :onset, :decline, :duration)
        idx = metrics[fm]
        outmean[idx]  = mean(getfield(ev, fm))
    end
    E = ev.dtype
    outmean[metrics[:maxes]] = mhcsminimax(E(ev.minimaxes)) 
    outmean[metrics[:frequency]] = length(ev.duration) / lfy
    outmean[metrics[:days]] = sum(ev.duration) / lfy 
    outmean
end

"""
    Helper for the annualmetrics to deal with the year start and end ranges.
"""
function _yrdate(mdate::StepRange{Date}, evst)
    cst, cse, cols = evst
    mcste = map((x, y) -> unique(year.(vcat(x,y))), cst, cse)
    myr = map((x, y) -> unique(year.(vcat(mdate[x], mdate[y]))), cst, cse)
    myst = map(x -> year.(mdate[x]), cst)
    myse = map(x -> year.(mdate[x]), cse)
    return mcste, myr, myst, myse
end

"""
    Compute the metrics for each year in the desired period.
"""
function annualmetrics(ev::MEvents{Vector{Vector{T}}}, indices, mdate, evst)
    mapcste, mapyr, mapyst, mapyse = _yrdate(mdate, evst)
    _, _, cols = evst
    lfy = (length ∘ unique)(year.(mdate))
    CIx, nCIx, x, y = indices
    z =  8 #length(metrics)
    E = ev.dtype
    outannual = ntuple(_ -> Array{T, 3}(undef, x, y, lfy), z)
    for i in eachindex(outannual)
        outannual[i][nCIx, :] .= T(NaN)
    end
    for (h, cx, mz, my, mt, me) in zip(cols, CIx, mapcste, mapyr, mapyst, mapyse)
        for (i, yr) in zip(mz, my)
            yx = findall(yr .== mt) 
            if isempty(yx) 
                yx = findall(yr .== me)
            end
            for fm in (:means, :sums, :onset, :decline, :duration)
                idx = metrics[fm]
                outannual[idx][cx, i]  = mean(getfield(ev, fm)[h][yx])
            end
            outannual[metrics[:maxes]][cx, i] = mhcsminimax(E(ev.minimaxes[h][yx]))
            outannual[metrics[:frequency]][cx, i] = convert(T, length(yx))
            outannual[metrics[:days]][cx, i] = convert(T, length(ev.duration[h][yx]))
        end
    end
    outannual
end

# vector version
function annualmetrics(ev::MEvents{Vector{T}}, indices, mdate, evst)
    mapcste, mapyr, mapyst, mapyse = _yrdate(mdate, evst)
    lfy = (length ∘ unique)(year.(mdate))
    z =  8 #length(metrics)
    E = ev.dtype
    outannual = ntuple(_ -> Vector{T}(undef, lfy), z)
    for (mz, my, mt, me) in zip(mapcste, mapyr, mapyst, mapyse)
        for (i, yr) in zip(mz, my)
            yx = findall(yr .== mt) 
            if isempty(yx) 
                yx = findall(yr .== me)
            end
            # NOTE: Added
            for fm in (:means, :sums, :onset, :decline, :duration)
                idx = metrics[fm]
                outannual[idx][i]  = mean(getfield(ev, fm)[yx])
            end
            outannual[metrics[:maxes]][i] = mhcsminimax(E(getfield(ev, :minimaxes)[yx]))
            outannual[metrics[:frequency]][i] = convert(T, length(yx))
            outannual[metrics[:days]][i] = length(getfield(ev, :duration)[yx])
        end
    end
    outannual
end

"""
    Return the coefficient, error, intercept, R² and p-value.

"""
function trend(outannual::NTuple{N, Array{T, 3}}, indices) where N
    CIx, nCIx, x, y = indices
    X = 1:size(first(outannual), 3)
    z = length(outannual)
    outpvalue = ntuple(_ -> Matrix{T}(undef, x, y), z)
    outcoeff = ntuple(_ -> Matrix{T}(undef, x, y), z)
    outrsqd = ntuple(_ -> Matrix{T}(undef, x, y), z)
    outerror_coeff = ntuple(_ -> Matrix{T}(undef, x, y), z)
    outintercept = ntuple(_ -> Matrix{T}(undef, x, y), z)
    # outpvalue = Array{T, 3}(undef, x, y, z)
    # outcoeff = similar(outpvalue)
    # outrsqd = similar(outpvalue)
    for i in 1:z
        for ci in CIx
            # outcoeff[i][ci],
            # outrsqd[i][ci], outpvalue[i][ci], = linreg(X, outannual[i][ci,:])
            # V2 if linreg returns without pvalues
            outlg = linreg(X, outannual[i][ci,:])
            outpvalue[i][ci], _ = _pvalue(outlg)
            outcoeff[i][ci] = getindex(outlg, 1)
            outintercept[i][ci] = getindex(outlg, 2)
            outrsqd[i][ci] = getindex(outlg, 3)
            outerror_coeff[i][ci] = getindex(outlg, 4)
        end
        outcoeff[i][nCIx] .= T(NaN)
        outpvalue[i][nCIx] .= T(NaN)
        outrsqd[i][nCIx] .= T(NaN)
    end
    outcoeff, outerror_coeff, outrsqd, outintercept, outpvalue 
end

function trend2(outannual::NTuple{N, Array{T, 3}}, indices) where N
    CIx, nCIx, x, y = indices
    X = 1:size(first(outannual), 3)
    z = length(outannual)
    outpvalue = Array{T, 3}(undef, x, y, z)
    outcoeff = similar(outpvalue)
    outrsqd = similar(outpvalue)
    outerror_coeff = similar(outpvalue)
    outintercept = similar(outpvalue)
    for i in 1:z
        for ci in CIx
            # V2 if linreg returns without pvalues
            outlg = linreg(X, outannual[i][ci,:])
            outpvalue[ci, i], _ = _pvalue(outlg)
            outcoeff[ci, i] = getindex(outlg, 1)
            outintercept[ci, i] = getindex(outlg, 2)
            outrsqd[ci, i] = getindex(outlg, 3)
            outerror_coeff[ci, i] = getindex(outlg, 4)
        end
        outcoeff[nCIx, i] .= T(NaN)
        outpvalue[nCIx, i] .= T(NaN)
        outrsqd[nCIx, i] .= T(NaN)
        outerror_coeff[nCIx, i] .= T(NaN)
        outintercept[nCIx, i] .= T(NaN)
    end
    outcoeff, outerror_coeff, outrsqd, outintercept, outpvalue 
end

# Vectorversion Output: each metric is a tuple of vectors
function trend(outannual::NTuple{N, Vector{T}}, indices) where N
    # CIx, nCIx, x, y = indices
    X = 1:size(first(outannual), 1)
    z = length(outannual)
    outpvalue = Vector{T}(undef, z)
    outcoeff = Vector{T}(undef, z)
    outrsqd = Vector{T}(undef, z)
    outintercept = Vector{T}(undef, z)
    outerror_coeff = Vector{T}(undef, z)
    for i in 1:z
        # outcoeff[i],
        # outrsqd[i], outpvalue[i], = linreg(X, outannual[i])
        # V2 if linreg returns without pvalues
        outlg = linreg(X, outannual[i])
        outpvalue[i], _ = _pvalue(outlg)
        outcoeff[i] = getindex(outlg, 1)
        outintercept[i] = getindex(outlg, 2)
        outrsqd[i] = getindex(outlg, 3)
        outerror_coeff[i] = getindex(outlg, 4)
    end
    outcoeff, outerror_coeff, outrsqd, outintercept, outpvalue 
end


function linreg(x::AbstractVector, y::AbstractVector)
    # Alexander Barth github.com/AlexanderBarth

    # remove NaNs
    ind = .!isnan.(x) .& .!isnan.(y)
    x = x[ind]
    y = y[ind]

    xm = mean(x)
    xa = x .- xm
    ym = mean(y)
    ya = y .- ym

    ssxy = xa' * ya
    ssxx = xa' * xa
    ssyy = ya' * ya

    b = ssxy / ssxx
    a = ym - b * xm
    r2 = ssxy^2 / (ssxx * ssyy)

    n = length(x)

    sigma_e = sqrt(max((ssyy - b^2 * ssxx) / (n - 2), 0))
    sigma_a = sigma_e * sqrt(sum(x .^ 2) / (n * ssxx))
    sigma_b = sigma_e / sqrt(ssxx)


    #= Calculate the F-statistic
    f_stat_a = (ssxy^2 / ssxx) / (ssyy - ssxy^2 / ssxx) * (n - 2)

    # Calculate the p-value for the F-statistic
    p_value_fa = 1 - cdf(FDist(1, n - 2), f_stat_a)

    # Calculate the t-statistic for the slope (b)
    t_stat_b = b / sigma_b

    # Calculate the p-value for the t-statistic
    p_value_ta = 2 * (1 - cdf(TDist(n - 2), abs(t_stat_b)))

    # Calculate T-Distribution p-value
    t_stat = b / sigma_b
    p_value_tb = 2 * (1 - cdf(TDist(n - 2), abs(t_stat)))

    # Calculate F-Distribution p-value
    f_stat_b = r2 / (1 - r2) * (n - 2) / 1
    p_value_f = 1 - cdf(FDist(1, n - 2), f_stat_b)
    =#
    return a, b, r2, sigma_a, sigma_b, sigma_e, convert(typeof(b), n)
end

function _pvalue(olg::NTuple{N, <:AbstractFloat}) where N
    b = olg[2]
    r2 = olg[3]
    sigma_b = olg[5]
    n = convert(Int, olg[7])

    # Calculate F-Distribution p-value
    f_stat = r2 / (1 - r2) * (n - 2) / 1
    p_value_f = 1 - cdf(FDist(1, n - 2), f_stat)

    # Calculate T-Distribution p-value
    t_stat = b / sigma_b
    p_value_t = 2 * (1 - cdf(TDist(n - 2), abs(t_stat)))
    
    return p_value_f, p_value_t
end

####
# Interfaces 
####
for f in (:pvalues, :coeffs, :rsquared)
    @eval begin
        function $f(am::MHCMetrics, metric)
            idx = get(metrics, metric, metric)
            <:(typeof(idx), Integer) || throw(KeyError(idx))
            xx = getfield(am, Symbol($f))
            ndims(xx) == 1 ? getindex(xx, idx) : getindex(xx, :,:,idx)
        end
    end
end

function mymetric(ev::EventsFull, metric)
    fd = first(propertynames(ev), 7)
    if metric in fd
        return reduce(vcat, getfield(ev, metric))
    else
        throw(KeyError(metric))
        @info "Perhaps you mean one of ", print(keys(metric))
    end
end

function mymetric(ev::EventsFull)
    # Return a vector of vectors
    vps = first(propertynames(ev), 7)
    [reduce(vcat, getfield(ev, t)) for t in vps]
end

    # I think you could stack it outside as in stack(mymetric(ev)) as (stack ∘ mymetrics)(ev)
 
# NOTE: We should also return the start and end dates 
# TODO: Return the anomalies and categories
# Change MEvents to EventsFull mostly globally
function mymetric(ev::EventsFull, indices)
    ids = mymetric(ev)
    cix = getindex(indices, 1)
    # the number of events per pixel
    drs = length.(getfield(ev, 1))
    ix = TI[]
    iy = TI[]
    for (q, s) in zip(drs, cix)
        append!(ix, repeat([Tuple(s)[1]], q))
        append!(iy, repeat([Tuple(s)[2]], q))
    end
    @assert length(ix) == length(iy) == sum(drs)
    stack([ids..., ix, iy])
end

function mymetric(mm::MExtreme, indices, arg::Symbol)
    # return clim, thresh or exceedance as 3-D arrays
    # mapping the argument to the fieldnames of the mm

    arg in fieldnames(typeof(mm)) || throw(KeyError(arg))
    gh = getfield(mm, arg)

    # NOTE: for clim and threshold, should we return 366 or leave it to the user?
    return _outarrays(gh, indices)
end

function mymetric(mm::MHCMetrics, field, metric)
    field in (:means, :annuals) || throw(error(field, " not a valid option. Try `:means` or `:annuals`"))
    idx = get(metrics, metric, metric)
    <:(typeof(idx), Integer) || throw(KeyError(idx))
    means = getfield(mm, field)
    outm = ndims(means) == 1 || isa(means, NTuple) ? getindex(means, idx) : getindex(means, :,:,idx)
    return outm 
end

"""
    Helper to return Matrix as 3-D array. 
"""
function _outarrays(gg::Matrix{T}, indices)
    CIx, nCIx, x, y = indices
    # Each column in gg is a pixel, and each row is a date. 
    oclim = Array{T, 3}(undef, x, y, size(gg, 1))
    for (col, cx) in enumerate(CIx)
        oclim[cx, :] = gg[:, col]
    end
    oclim[nCIx, :] .= NaN
    oclim
end

"""
    Return the Events and the labels (starts and end positions) of the event.
"""
function mevents(ms, mindur, maxgap, indices)
    mindur, maxgap = convert(Vector{TI}, [mindur, maxgap])
    evst = _mylabel(ms, mindur, maxgap)
    ev = anomsa(ms, evst, indices)
    return ev, evst
end

"""
    Return the computed metrics - means, annuals and linear regression outputs as MHCMetrics.
"""
function mmetrics(ev, evst, indices, mdate)
    mm = meanmetrics2(ev::EventsFull{Vector{Vector{T}}}, indices, mdate)
    am = annualmetrics(ev::EventsFull{Vector{Vector{T}}}, indices, mdate, evst)
    M = length(am)
    N = ndims(mm)
    tr = trend2(am::NTuple{M, Array{T, N}}, indices)
    MHCMetrics{T, M, N, typeof(mm)}(am, mm, tr...)
end

#= STEPS
    # 1. calculate the MHW/MCS
    ms, indices = mextreme(sst, sstdate::StepRange, mhwdate::StepRange, clmdate::StepRange, event=:mhw; window=5, smoothwindow=31, threshold=nothing)
    # 2. label + events
    ev, lb = mevents(ms, mindur, maxgap, indices)
    # 3. pixel, annual and trend
    mht = mmetrics(ev, lb, indices, mdate)
=#
