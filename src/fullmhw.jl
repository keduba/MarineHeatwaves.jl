
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
                     :decline => 8)

# abstract type MExtreme end
abstract type MWrapper end
abstract type MEvents{TD<:AbstractVector} end #<:AbstractFloat}} end
abstract type MExtreme{TD<:AbstractVecOrMat{<:AbstractFloat}, V<:BitArray} end # MCS/mhw


struct MCWrapper{T}
    anom::Vector{T}
    threshanom::Vector{T}
    onset::T
    decline::T
end

struct MHWrapper{T}
    anom::Vector{T}
    threshanom::Vector{T}
    onset::T
    decline::T
end

struct Events{TE} <: MEvents{TE}
    means::TE
    minimaxes::TE
    onset::TE
    decline::TE
    duration::TE
    sums::TE
    dtype::DataType
end

# # FIX: to remove
# struct EventsCS{T} <: MEvents
#     means::Array{T, 1}
#     minimaxes::Array{T, 1}
#     onset::Array{T, 1}
#     decline::Array{T, 1}
#     duration::Array{T, 1}
#     sums::Array{T, 1}
#     dtype::DataType
# end

# struct Eventx{TA} <: MEv{TA}
#     where TA <: AbstractVector
#     sy::TA{T}
# end

struct EventHW{TA} <: MEvents{TA}
    # where TA <: AbstractVector{T}
    maxes::TA
end

struct EventCS{TA} <: MEvents{TA}
    # where TA <: AbstractVector{T}
    maxes::TA
end
 
# takes care of mhw/mcs as sample
# struct MHWn{TA,V} <: MHExt{TA,V}
#     where TA <: AbstractVecOrMat
#     temp::TA{T}
#     clim::TA{T}
#     thresh::TA{T}
#     exceeds::V
#     edtype::DataType
# end

struct MCS{TA,V} <: MExtreme{TA,V}
    # where TA <: AbstractVecOrMat
    temp::TA
    clim::TA
    thresh::TA
    exceeds::V
    edtype::Type{EventCS}
end

struct MHW{TA, V} <: MExtreme{TA,V}
    # where TA <: AbstractVecOrMat
    temp::TA
    clim::TA
    thresh::TA
    exceeds::V
    edtype::Type{EventHW}
end

function MHW(temp, clim, thresh)
    exc = _excess(temp, thresh)
    MHW(temp, clim, thresh, exc, EventHW)
end

function MCS(temp, clim, thresh)
    exc = _excess(thresh, temp)
    MCS{typeof(temp), typeof(clim), typeof(thresh), typeof(exc)}(temp, clim, thresh, exc, EventCS)
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

# we expect Events to currently return Vector{Vector{T}} so the following function is targeting the field we want.
mhcsminimax(x::EventHW{Vector{T}}) = maximum(x.maxes)
mhcsminimax(x::EventCS{Vector{T}}) = minimum(x.maxes)

function mhcsminimax(x::EventHW{Vector{Vector{T}}})
    return @inbounds [maximum(ia) for ia in x.maxes]
end

function mhcsminimax(x::EventCS{Vector{Vector{T}}})
    return @inbounds [minimum(ia) for ia in x.maxes]
end

function timeindices(sstdate::StepRange{Date, Day}, evtdate::StepRange{Date, Day})
    si::TI = findfirst(isequal(first(evtdate)), sstdate)
    ei::TI = findfirst(isequal(last(evtdate)), sstdate)
    return range(si, ei)
end

seamask(sst::AbstractArray, dims) = dropdims(count(!isnan, sst, dims=dims) .> 0, dims=dims)

function _subtemp(sst::AbstractArray, mhwix::UnitRange, clmix::UnitRange)
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

# FIX: Wrapper for event
# 1. Base entry function that allows to select if it's mhw or mcs
# 2. From this entry point, I return the appropriate wrapper for the event
# 3. This should allow me to return the indices and the wrapper
#
# Works for MarineHW tested. okay.
# MarineHW(args) = MHW(subtemp(args...)...)
# MarineCS(args)::MCS = MCS(subtemp(args...)...)

leapyearday(mts::Date)::TI = dayofyear(mts) > 59 && !isleapyear(mts) ? dayofyear(mts) + 1 : dayofyear(mts)

_excess(tp, th) = tp .> th

excess(ms::MExtreme) = ms.exceeds

const events = Dict(:mhw => (MHW, 0.9),
                    :mcs => (MCS, 0.1))

# TODO: this object returns the clim and thresh by returning a MHW or MCS object so should reflect the name
function subtemp(sst, sstdate::StepRange, mhwdate::StepRange, clmdate::StepRange, event=:mhw; window=5, smoothwindow=31, threshold=nothing)
    window::TI = convert(TI, window)
    smoothwindow = convert(TI, smoothwindow)
    ME, mthreshold = events[event]
    threshold::T = isnothing(threshold) ? convert(T, mthreshold) : convert(T, threshold)
    mhwix = timeindices(sstdate, mhwdate)
    clmix = timeindices(sstdate, clmdate)
    mlyd = leapyearday.(sstdate[mhwix])
    clyd = leapyearday.(sstdate[clmix])
    mhtemp, ctemp, CIs = _subtemp(sst, mhwix, clmix)
    clima, thresh = climthresh(ctemp, clyd, mlyd, window, smoothwindow, threshold)
    return ME(mhtemp, clima, thresh), CIs
end

function climthresh(cmst::Matrix{T}, clyd, mlyd, win::TI, swindow::TI, threshold::T) # wrap?
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

function climthresh(cmst::Vector{T}, clyd, mlyd, win::TI, swindow::TI, threshold::T)
    # T =  eltype(cmst)
    # clim = Vector{T}(undef, length(mlyd))
    # thresh = similar(clim)
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
    clim = moving_means(cv, swindow, mlyd)
    thresh = moving_means(th, swindow, mlyd)
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

function _mylabel(ms::MExtreme{Vector{T}}, mindur, maxgap) 
    sty = excess(ms)
    stb = sparse(diff(sty))
    cst = TI[] 
    cse = TI[]
    for i in nzrange(stb, c)
        if isequal(nonzeros(stb)[i], -1)
            push!(cse, rowvals(stb)[i])
        else
            push!(cst, rowvals(stb)[i])
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

function anomsa(m::MHW{Matrix{T}}, c, args...) 
    outs = _anomsa(m.temp[:,c], m.clim[:,c], m.thresh[:,c], args...)
    return MHWrapper(outs...)
end

function anomsa(m::MCS{Matrix{T}}, c, args...) 
    outs = _anomsa(m.temp[:,c], m.clim[:,c], m.thresh[:,c], args...)
    return MCWrapper(outs...)
end

function anomsa(m::MCS{Vector{T}}, args...) 
    outs = _anomsa(m.temp, m.clim, m.thresh, args...)
    return MCWrapper(outs...)
end

function anomsa(m::MHW{Vector{T}}, args...) 
    outs = _anomsa(m.temp, m.clim, m.thresh, args...)
    return MHWrapper(outs...)
end

function _onset(atod::MWrapper, mst)
    fan = first(atod) 
    nmx = mhcsminimax(atod) 
    ngx = mhcsarg(atod) 
    lnmx = nmx - 0.5(fan + atod.onset)
    snmx = nmx - fan
    mst > 1 ? /(lnmx, (ngx + 0.5)) : /(snmx, ngx)
end

function _decline(atod::MWrapper, mse, lnx)
    lan = last(atod)
    nmx = mhcsminimax(atod) 
    ngx = mhcsarg(atod) 
    wsnx::T = length(atod) - ngx
    lnmx = nmx - 0.5(lan + atod.decline)
    snmx = nmx - lan
    mse < lnx ? /(lnmx, (wsnx + 0.5)) : ngx == lnx ? snmx : /(snmx, wsnx)
end

function anomsa(m::MExtreme{Matrix{T}}, evst, indices) where M<:AbstractMatrix{T}
    # MW, Ev = typeof(m) == MHW{Matrix{T}} ? (MHWrapper, EventHW) : (MCWrapper, EventCS)
    CIx, nCIx, x, y = indices
    mst, mse, cols = evst
    lm::TI = size(m.temp, 1)
    mt::TI = 6
    outemp, outhsh, outcat = ntuple(_ -> Array{T, 3}(undef, x, y, lm), 3)
    onsan, decan, means, cums, maxes, durs = ntuple(_ -> [Vector{T}(undef,m) for m in length.(mst)], mt)
    for (c, cst, cse) in zip(cols, mst, mse)
        for (d, (st, se)) in enumerate(zip(cst, cse))
            atod = anomsa(m, c, st, se, lm)
            # atod = _anomsa(m.temp[:, c], m.clim[:, c], m.thresh[:, c], st, se, lm)
            outemp[CIx[c], st:se] = atod.anom
            outhsh[CIx[c], st:se] = atod.threshanom
            outcat[CIx[c], st:se] .= categorys(atod)
            onsan[c][d] = _onset(atod, st)
            decan[c][d] = _decline(atod, se, lm)
            means[c][d], cums[c][d], maxes[c][d], durs[c][d] = _events(atod)# vars[c][d]
        end
    end
    for bl in (outemp, outhsh, outcat)
        bl[nCIx, :] .= NaN
    end
    # if using mutable struct, we can also wrap outemp, outhsh, outcat
    return outemp, outhsh, outcat, Events(onsan, decan, means, cums, maxes, durs, m.edtype)
end

function anomsa(m::MExtreme, evst, indices)
    # MW, Ev = typeof(m) == MHW{Vector{T}} ? (MHWrapper, EventHW) : (MCWrapper, EventCS)
    mst, mse = evst
    lm::TI = size(m.temp, 1)
    mt::TI = 6
    outemp, outhsh, outcat = ntuple(_ -> Vector{T}(undef, lm), 3)
    onsan, decan, means, cums, maxes, durs = ntuple(_ -> Vector{T}(undef,  length(mst)), mt)
    # for (c, cst, cse) in zip(cols, mst, mse)
        for (d, (st, se)) in enumerate(zip(mst, mse))
            atod = _anomsa(m.wdtype, m.temp, m.clim, m.thresh, st, se, lm)
            outemp[st:se] = atod.anom
            outhsh[st:se] = atod.threshanom
            outcat[st:se] .= categorys(atod)
            onsan[c][d] = _onset(atod, st)
            decan[c][d] = _decline(atod, se, lm)
            means[c][d], cums[c][d], maxes[c][d], durs[c][d] = _events(atod)# vars[c][d]
        end
    # end
    # for bl in (outemp, outhsh, outcat)
    #     bl[nCIx, :] .= NaN
    # end
    # if using mutable struct, we can also wrap outemp, outhsh, outcat
    return outemp, outhsh, outcat, Events(onsan, decan, means, cums, maxes, durs, m.edtype)
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

function meanmetrics3(ev::MEvents, indices, mdate)
    CIx, nCIx, x, y = indices
    lfy = (length ∘ unique)(year.(mdate))
    z = length(metrics)
    outmean = ntuple(_ -> Matrix{T}(undef, x, y), z) 
    for i in eachindex(outmean)
        outmean[i][nCIx] .= T(NaN)
    end
    # check the metrics dictionary to ensure order of variables
    for fm in (:means, :sums, :onset, :decline, :duration)
        idx = metrics[fm]
        outmean[idx][CIx]  = mean.(ev.fm)
        setindex!(outmean[idx], CIx, mean.(ev.fm))
    end
    E = ev.dtype
    # em = ev.dtype{Vector{T}}(ev.minimaxes)
    # outmean[6][CIx] = mhcsminimax(em)
    outmean[6][CIx] = mhcsminimax(E(ev.minimaxes)) 
    outmean[7][CIx] = @. length(ev.durs) / lfy # frequency
    outmean[z][CIx] = @. sum(ev.durs) / lfy #days
    outmean
end

function meanmetrics3(ev::MEvents, indices, mdate)
    # CIx, nCIx, x, y = indices
    # VEctor version
    lfy = (length ∘ unique)(year.(mdate))
    z = length(metrics)
    outmean = Vector{T}(undef, z) 
    # for i in eachindex(outmean)
    #     outmean[i][nCIx] .= T(NaN)
    # end
    # check the metrics dictionary to ensure order of variables
    for fm in (:means, :sums, :onset, :decline, :duration)
        idx = metrics[fm]
        outmean[idx]  = mean(ev.$fm)
        setindex!(outmean[idx], mean(ev.$fm))
    end
    E = ev.dtype
    outmean[6] = mhcsminimax(E(ev.minimaxes)) # mhcminimax(ev.maxes)
    outmean[7] = length(ev.durs) / lfy # frequency
    outmean[z] = sum(ev.durs) / lfy #days
    outmean
end

function anumets4(ev::MEvents, indices, mdate, evst)
    cst, cse, cols = evst
    # f = isa(EventsHW, ev) ? maximum : minimum 
    mapcste = map((x, y) -> unique(year.(vcat(x,y))), cst, cse)
    mapyr = map((x, y) -> unique(year.(vcat(mdate[x], mdate[y]))), cst, cse)
    mapyst = map(x -> year.(mdate[x]), cst)
    mapyse = map(x -> year.(mdate[x]), cse)
    lfy = (length ∘ unique)(year.(mdate))
    CIx, nCIx, x, y = indices
    z =  length(metrics)
    E = ev.dtype
    outannual = ntuple(_ -> Array{T, 3}(undef, x, y, lfy), z)
    for i in eachindex(outannual)
        outannual[i][nCIx, :] .= T(NaN)
    end
    # for j in eachindex(outannual)# To remove this outer loop
        for (h, cx, mz, my, mt, me) in zip(cols, CIx, mapcste, mapyr, mapyst, mapyse)
            for (i, yr) in zip(mz, my)
                yx = findall(yr .== mt) 
                if isempty(yx) 
                    yx = findall(yr .== me)
                end
                # NOTE: Added
                for fm in (:means, :sums, :onset, :decline, :duration)
                    idx = metrics[fm]
                    # outmean[idx][CIx]  = mean.(ev.fm)
                    setindex!(outannual[idx], cx, i, mean(ev.fm[h][yx]))
                end
                setindex!(outannual[6], cx, i, mhcsminimax(E(ev.minimaxes[h][yx])))
                setindex!(outannual[7], cx, i, convert(T, length(yx)))
                setindex!(outannual[z], cx, i, length(ev.durations[h][yx])) 
                # if j == 6 
                #     # TODO: Fix mini-maximum and type passed to allow indexing
                #     # create a wrapper that's a type of array
                #     # outannual[j][cx, i] = maximum(ev[j][h][yx]) # WARN: change maximum to minimax
                #     outannual[j][cx, i] = f(ev.maxes[h][yx]) # WARN: change maximum to minimax
                # elseif j == 7 # frequency maybe
                #     outannual[j][cx, i] = convert(T, length(yx))
                # elseif j == 8 # durations
                #     outannual[j][cx, i] = length(ev.durations[h][yx])
                # else
                #     outannual[j][cx, i] = mean(ev[j][h][yx])
                # end
            end
        end
    # end
    outannual
end


# vector version
function anumets4(ev::MEvents, indices, mdate, evst)
    cst, cse = evst
    mapcste = map((x, y) -> unique(year.(vcat(x,y))), cst, cse)
    mapyr = map((x, y) -> unique(year.(vcat(mdate[x], mdate[y]))), cst, cse)
    mapyst = map(x -> year.(mdate[x]), cst)
    mapyse = map(x -> year.(mdate[x]), cse)
    lfy = (length ∘ unique)(year.(mdate))
    z =  length(metrics)
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
                    # outmean[idx][CIx]  = mean.(ev.fm)
                    setindex!(outannual[idx], i, mean(ev.fm[yx]))
                end
                setindex!(outannual[6], i, mhcsminimax(E(ev.minimaxes[yx])))
                setindex!(outannual[7], i, convert(T, length(yx)))
                setindex!(outannual[z], i, length(ev.durations[yx])) 
            end
        end
    outannual
end

# Input: outannual is a tuple of 3-D arrays
# Output: each metric is a tuple of Matrices
function trend(outannual::NTuple{8, Array{T, 3}}, indices)
    CIx, nCIx, x, y = indices
    X = 1:size(first(outannual), 3)
    z = 8 # no of metrics length(outannual)
    outpvalue = ntuple(_ -> Matrix{T}(undef, x, y), z)
    outcoeff = ntuple(_ -> Matrix{T}(undef, x, y), z)
    outrsqd = ntuple(_ -> Matrix{T}(undef, x, y), z)
    for i in 1:z
        for ci in CIx
            outcoeff[i][ci],
            outrsqd[i][ci], outpvalue[i][ci], = linreg(X, outannual[i][ci,:])
            # V2 if linreg returns without pvalues
            # outlg = linreg(X, outannual[i][ci,:])
            # outpvalue[i][ci] = _pvalue(outlg)
            # outcoeff[i][ci] = getindex(outlg, 2)
            # outrsqd[i][ci] = getindex(outlg, 3)
        end
    end
    outpvalue, outcoeff, outrsqd
end

# Vectorversion Output: each metric is a tuple of vectors
function trend(outannual::NTuple{8, Vector{T}}, indices)
    # CIx, nCIx, x, y = indices
    X = 1:size(first(outannual), 1)
    z = 8 # no of metrics length(outannual)
    outpvalue = Vector{T}(undef, z)
    outcoeff = Vector{T}(undef, z)
    outrsqd = Vector{T}(undef, z)
    for i in 1:z
        # for ci in CIx
            outcoeff[i],
            outrsqd[i], outpvalue[i], = linreg(X, outannual[i])
            # V2 if linreg returns without pvalues
            # outlg = linreg(X, outannual[i])
            # outpvalue[i] = _pvalue(outlg)
            # outcoeff[i] = getindex(outlg, 2)
            # outrsqd[i] = getindex(outlg, 3)
        # end
    end
    outpvalue, outcoeff, outrsqd
end

####
# Here begin the Interfaces to access the results and outputs
# NOTE: the internal variable `indices` has to be stored somewhere in one of the outward structs if it is to be used to return 2D to 3D outputs
####

for (f, fn) in ((:pvalues, :outpvalue), (:coeffs, :outcoeff), (:rsquared, :outrsquared))
    @eval begin
        function $f(am::MHWCSO, metric)
            idx = get(metrics, metric, metric)
            if typeof(idx) == Int
                return getindex(am.$fn, idx)
            else
                throw(KeyError(idx))
            end
        end
    end
end

function mymetric(ev::MEvents, metric)
    fd = fieldnames(ev)
    if metric in fd
        return reduce(vcat, ev.$metric)
    else
        throw(Keyerror(metric))
        @info "Perhaps you mean one of ", print(keys(metric))
    end
end


function mymetric(ev::MEvents)
    # Return a vector of vectors
    [reduce(vcat, ev.t) for t in fieldnames(ev)]
end

    # I think you could stack it outside as in stack(mymetric(ev)) as (stack ∘ mymetrics)(ev)
 
function mymetric(ev::MEvents, indices)
    # Return a the metrics as vector 
    mymetric(ev)
    cix = getindex(indices, 1)
    # the number of events per pixel
    drs = length.(ev.fieldnames(ev)[1])
    ix = TI[]
    iy = TI[]
    for (q, s) in zip(drs, cix)
        append!(ix, repeat([Tuple(s)[1]], q))
        append!(iy, repeat([Tuple(s)[2]], q))
    end
    @assert length(ix) == length(iy) == sum(drs)
end

function mymetric(mm::MExtreme, indices, arg::Symbol=:anom)
    # return clim, thresh or exceedance as 3-D arrays
    # Matrix to Array conversion: require indices 
    # mapping the argument to the fieldnames of the mm
    Dict(:anom => :anoma,
         :thresh => :threshold,
         :exceeds => :exceedance)

    arg in fieldnames(typeof(mm)) || throw(KeyError(arg))
    gh = getfield(mm, arg)
    return _outarrays(gh, indices)
end

function mymetric(mm::MHWCSO, arg::Symbol=:anom)
    # return clim, thresh or exceedance as 3-D arrays
    # Matrix to Array conversion: require indices 
    # mapping the argument to the fieldnames of the mm
    Dict(:anom => :anoma,
         :thresh => :threshold,
         :exceeds => :exceedance)

    arg in fieldnames(typeof(mm)) || throw(KeyError(arg))
    return gh = getfield(mm, arg)
end

function _outarrays(gg::Matrix{T}, indices)
    CIx, nCIx, x, y = indices
    oclim = Array{T, 3}(undef, x, y, size(gg, 1))
    for (col, cx) in enumerate(CIx)
        oclim[cx, :] = gg[:, col]
    end
    oclim[nCIx, :] .= NaN
    oclim
end

