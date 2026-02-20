import Statistics: mean, quantile, std
import Distributions: cdf, FDist, TDist
using Dates
using NCDatasets
using SparseArrays

const T = Float32
const TI = Int16

const SymOrString = Union{Symbol, String}
const metrics = Dict(:means => 1,
                     :maxes => 2,
                     :sums => 3,
                     :duration => 4,
                     :frequency => 5,
                     :days => 6,
                     :onset => 7,
                     :decline => 8)

# The Abstract types

abstract type MWrapper end
abstract type MEvents{TE <: AbstractVector} end
abstract type MExtreme{TA <: AbstractVecOrMat{<: AbstractFloat}, V <: BitArray} end


# Input structs

struct EventHW{TE} <: MEvents{TE}
    maxes::TE
end

struct EventCS{TE} <: MEvents{TE}
    maxes::TE
end

struct MCS{TA, V} <: MExtreme{TA, V}
    temp::TA
    clim::TA
    thresh::TA
    exceeds::V
    edtype::Type{EventCS}
end

function MCS(temp::AT, clim::AT, thresh::AT) where AT<:VecOrMat{<:AbstractFloat}
    temp = convert(typeof(clim), temp)
    exc = _excess(thresh, temp)
    MCS{typeof(temp), typeof(exc)}(temp, clim, thresh, exc, EventCS)
end

struct MHW{TA, V} <: MExtreme{TA, V}
    temp::TA
    clim::TA
    thresh::TA
    exceeds::V
    edtype::Type{EventHW}
end

function MHW(temp::AT, clim::AT, thresh::AT) where AT<:VecOrMat{<:AbstractFloat}
    temp = convert(typeof(clim), temp)
    exc = _excess(temp, thresh)
    MHW{typeof(temp), typeof(exc)}(temp, clim, thresh, exc, EventHW)
end

const events = Dict(:mhw => (MHW, 0.9),
                    :mcs => (MCS, 0.1))

# Internal structs

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

# struct Events{TE} <: MEvents{TE}
# end

struct Events{TE} <: MEvents{TE}#<:AbstractVector, Ti<:AbstractVecOrMat}
    means::TE
    minimaxes::TE
    onset::TE
    decline::TE
    duration::TE
    sums::TE
    categorys::TE
    tpanom::VecOrMat{T}
    thanom::VecOrMat{T}
    dtype::Type{<:MEvents}
end

# Output structs

struct MHCMetricsm{TA, TV}
    annuals::TA
    means::TV
    coeffs::TV
    errors::TV
    rsquared::TV
    intercepts::TV
    pvalues::TV
end

# Base implementations for internal structs

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

"""
Return the range of date B in date A.

timeindices(dateA, dateB)::UnitRange
"""
function timeindices(sstdate::StepRange{Date, Day}, evtdate::StepRange{Date, Day})::UnitRange
    si = findfirst(isequal(first(evtdate)), sstdate)
    ei = findfirst(isequal(last(evtdate)), sstdate)
    (si == nothing || ei == nothing) && error("The date $(evtdate) is outside $(sstdate)")
    si, ei = TI.([si, ei])
    return range(si, ei)
end

seamask(sst::AbstractArray, dims) = dropdims(count(!isnan, sst, dims=dims) .> 0, dims=dims)

"""
Calculate the subset of the input array. For n-dimensional array input where n > 1, return n-1 dimension. 

Return the subset and indices of mask.
"""
function _subtemp(sst::AbstractArray, mhwix::UnitRange, clmix::UnitRange)
    N = ndims(sst)
    N ≤ 0 && error("0-dimensional data just can't work.")
    in(N, [1, 3]) || error("Expected 1 or 3 dimensions. Got $N dimensions.")
    xz = seamask(sst, N)

    if N == 1
        xz[1] || error("Please check your data. It appears to be all `NaN` or `missing`.")
    end

    if N == 3
        CIx = CartesianIndices(xz)[xz]
        @views begin
            msst = permutedims(sst[CIx, mhwix])
            csst = permutedims(sst[CIx, clmix])
        end
        x, y = Base.front(size(sst))
        # nCIx = setdiff(CartesianIndices((x, y)), CIx)
        return msst, csst, (CIx, x, y)
        # return msst, csst, (CIx, nCIx, x, y)
    else
        CIx = TI(xz[1])
        msst = sst[mhwix]
        csst = sst[clmix]
        return msst, csst, CIx
    end
end


leapyearday(mts::Date)::TI = dayofyear(mts) > 59 && !isleapyear(mts) ? dayofyear(mts) + 1 : dayofyear(mts)

_excess(tp, th) = tp .> th

excess(ms::MExtreme) = ms.exceeds

"""
Return a MExtreme (MHW or MCS).

"""
function mextreme(sst, sstdate::StepRange, mhwdate::StepRange, clmdate::StepRange; event=:mhw, window=5, smoothwindow=31, threshold=nothing, wrap=true)
    in(event, keys(events)) || error("`:$event` is not a valid event. Try `:mhw` or `:mcs`")
    ME, mthreshold = get(events, event, :mhw)
    threshold = isnothing(threshold) ? convert(T, mthreshold) : convert(T, threshold)
    window, smoothwindow = convert.(TI, [window, smoothwindow])
    mhwix = timeindices(sstdate, mhwdate)
    clmix = timeindices(sstdate, clmdate)
    mlyd = leapyearday.(sstdate[mhwix])
    clyd = leapyearday.(sstdate[clmix])
    mhtemp, ctemp, CIs = _subtemp(sst, mhwix, clmix)
    clima, thresh = climthresh(ctemp, clyd, mlyd, window, smoothwindow, threshold, wrap=wrap)
    return ME(mhtemp, clima, thresh), CIs
end

"""
Calculate and return the climatology mean and quantile.

"""
function climthresh(cmst::Matrix{T}, clyd::Vector{TI}, mlyd::Vector{TI}, window::Integer, smoothwindow::Integer, threshold::T; wrap::Bool=false)
    window::TI = convert(TI, window)
    smoothwindow = convert(TI, smoothwindow)
    clim = Matrix{T}(undef, length(mlyd), size(cmst, 2))
    thresh = similar(clim)
    uranges = urange(clyd, window)
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
        clim[:, n] = moving_means(cv, smoothwindow, mlyd, wrap=wrap)
        thresh[:, n] = moving_means(th, smoothwindow, mlyd, wrap=wrap)
    end
    clim, thresh
end

function climthresh(cmst::Vector{T}, clyd::Vector{TI}, mlyd::Vector{TI}, window::Integer, smoothwindow::Integer, threshold::T; wrap::Bool=true)
    window::TI = convert(TI, window)
    smoothwindow = convert(TI, smoothwindow)
    uranges = urange(clyd, window)
    cv = Vector{T}(undef, 366)
    th = similar(cv)
    for (m, ur) in enumerate(uranges)
        cvc = Iterators.flatten([cmst[i] for i in ur])
        cv[m] = mean(cvc)
        th[m] = quantile(cvc, threshold)
    end
    cv[60] = 0.5(cv[59] + cv[61])
    th[60] = 0.5(th[59] + th[61])
    clim = moving_means(cv, smoothwindow, mlyd, wrap=wrap)
    thresh = moving_means(th, smoothwindow, mlyd, wrap=wrap)
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

"""
Return the start and stop positions and pixel columns at which MExtremes were detected.

mylabel(me::MExtreme, minimum\\_duration, maximum\\_gap)

"""
function mylabel(ms::MExtreme{Matrix{T}}, mindur::TI, maxgap::TI)
    (mindur <= 0 || maxgap <=0) && error("The minimum duration or maximum gap cannot be less than 1.")
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

function mylabel(ms::MExtreme{Vector{T}}, mindur::TI, maxgap::TI)
    (mindur <= 0 || maxgap <= 0) && error("The minimum duration or maximum gap cannot be less than 1.")
    sty = excess(ms)
    stb = sparse(diff(sty))
    cst = TI[] 
    cse = TI[]
    for i in nzrange(stb, 1)
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

"""
Return the rate of onset for each event.

"""
function onset(atod::MWrapper, mst)
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
function decline(atod::MWrapper, mse, lnx)
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
function anomsam(m::MExtreme{Matrix{T}}, evst::Tuple)
    MW = m.edtype == Type{EventHW} ? MHWrapper : MCWrapper
    mst, mse, cols = evst
    lm::TI = size(m.temp, 1)
    mt::TI = 7
    outemp, outhsh = ntuple(_ -> Matrix{T}(undef, lm, length(cols)), 2)
    # outemp .= NaN; outhsh .= NaN
    onsan, decan, means, cums, maxes, durs, catso = ntuple(_ -> [Vector{T}(undef, m) for m in length.(mst)], mt)
    for (c, cst, cse) in zip(cols, mst, mse)
        for (d, (st, se)) in enumerate(zip(cst, cse))
            atod = MW(_anomsa(m.temp[:,c], m.clim[:,c], m.thresh[:,c], st, se, lm)...)
            outemp[st:se, c] = atod.anom
            outhsh[st:se, c] = atod.threshanom
            catso[c][d] = categorys(atod)
            onsan[c][d] = onset(atod, st)
            decan[c][d] = decline(atod, se, lm)
            means[c][d], cums[c][d], maxes[c][d], durs[c][d] = mevents(atod)
        end
    end
    return Events(means, maxes, onsan, decan, durs, cums, catso, outemp, outhsh, m.edtype)
end


function anomsav(m::MExtreme{Vector{T}}, evst::Tuple)
    MW = m.edtype == Type{EventHW} ? MHWrapper : MCWrapper
    mst, mse = evst
    lm::TI = size(m.temp, 1)
    mt::TI = 7
    outemp, outhsh = ntuple(_ -> Vector{T}(undef, lm), 2)
    # outemp .= NaN; outhsh .= NaN
    onsan, decan, means, cums, maxes, durs, catso = ntuple(_ -> Vector{T}(undef, length(mst)), mt)
    for (d, (st, se)) in enumerate(zip(mst, mse))
        atod = MW(_anomsa(m.temp, m.clim, m.thresh, st, se, lm)...)
        outemp[st:se] = atod.anom
        outhsh[st:se] = atod.threshanom
        cats = categorys(atod)
        # outcat[st:se] .= cats
        catso[d] = cats
        onsan[d] = onset(atod, st)
        decan[d] = decline(atod, se, lm)
        means[d], cums[d], maxes[d], durs[d] = mevents(atod)
    end
    return Events(means, maxes, onsan, decan, durs, cums, catso, outemp, outhsh, m.edtype)
end

_categorys(anom::Vector{T}, thsd::Vector{T}) = min(4, maximum(fld.(anom, thsd)))

categorys(a::MWrapper) = _categorys(a.anom, a.threshanom)

function mevents(anom::MWrapper)
    # Per Event Metrics
    means = mean(anom)
    cums = sum(anom)
    maxes = mhcsminimax(anom)
    durs = length(anom)
    return means, cums, maxes, durs
end


# potentially
"""
Return the mean of the metrics in each pixel.
"""
function meanmetricsm(ev::Events{Vector{Vector{T}}}, mdate::StepRange{Date, Day})
    lfy = (length ∘ unique)(year.(mdate))
    z = length(metrics)
    y = length(getfield(ev, :means))
    outmean = Matrix{T}(undef, y, z)
    # outmean .= NaN
    for fm in (:means, :sums, :onset, :decline, :duration)
        idx = metrics[fm]
        outmean[:, idx]  = mean.(getfield(ev, fm))
    end
    E = ev.dtype
    outmean[:, metrics[:maxes]] = mhcsminimax(E(ev.minimaxes))
    outmean[:, metrics[:frequency]] = @. length(ev.duration)/lfy
    outmean[:, metrics[:days]] = @. sum(ev.duration)/lfy
    outmean
end


function meanmetricsv(ev::Events{Vector{T}}, mdate::StepRange{Date, Day})
    # Vector version
    lfy = (length ∘ unique)(year.(mdate))
    z = length(metrics)
    outmean = Vector{T}(undef, z)
    # check the metrics dictionary to ensure order of variables
    for fm in (:means, :sums, :onset, :decline, :duration)
        idx = metrics[fm]
        outmean[idx]  = mean(getfield(ev, fm))
    end
    E = ev.dtype
    outmean[metrics[:maxes]] = mhcsminimax(E(ev.minimaxes))
    outmean[metrics[:frequency]] = length(ev.duration)/lfy
    outmean[metrics[:days]] = sum(ev.duration)/lfy
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
function annualmetricsm(ev::Events{Vector{Vector{T}}}, mdate::StepRange{Date, Day}, evst::Tuple)
    mapcste, mapyr, mapyst, mapyse = _yrdate(mdate, evst)
    _, _, cols = evst
    lfy = (length ∘ unique)(year.(mdate))
    y = length(getfield(ev, :means))
    z = length(metrics)
    E = ev.dtype
    # no of years * no of pixels * no of metrics
    # outannual = Array{T, 3}(undef, y, lfy, z)
    outannual = Array{T, 3}(undef, lfy, y, z)
    # outannual .= NaN
    for (c, mz, my, mt, me) in zip(cols, mapcste, mapyr, mapyst, mapyse)
        for (i, yr) in zip(mz, my)
            yx = findall(yr .== mt)
            if isempty(yx)
                yx = findall(yr .== me)
            end
            for fm in (:means, :sums, :duration, :onset, :decline)
                idx = metrics[fm]
                outannual[i, c, idx]  = mean(getfield(ev, fm)[c][yx])
            end
            outannual[i, c, metrics[:maxes]] = mhcsminimax(E(ev.minimaxes[c][yx]))
            outannual[i, c, metrics[:frequency]] = convert(T, length(yx))
            outannual[i, c, metrics[:days]] = convert(T, sum(ev.duration[c][yx]))
        end
    end
    outannual
end


# vector version
function annualmetricsv(ev::Events{Vector{T}}, mdate::StepRange{Date}, evst::Tuple)
    mapcste, mapyr, mapyst, mapyse = _yrdate(mdate, evst)
    lfy = (length ∘ unique)(year.(mdate))
    z =  length(metrics)
    E = ev.dtype
    # outannual = ntuple(_ -> Vector{T}(undef, lfy), z)
    # no of years * no of metrics
    outannual = Matrix{T}(undef, lfy, z)
    # outannual .= NaN
    for (mz, my, mt, me) in zip(mapcste, mapyr, mapyst, mapyse)
        for (i, yr) in zip(mz, my)
            yx = findall(yr .== mt)
            if isempty(yx)
                yx = findall(yr .== me)
            end
            for fm in (:means, :sums, :onset, :decline, :duration)
                idx = metrics[fm]
                outannual[i, idx]  = mean(getfield(ev, fm)[yx])
            end
            outannual[i, metrics[:maxes]] = mhcsminimax(E(getfield(ev, :minimaxes)[yx]))
            outannual[i, metrics[:frequency]] = convert(T, length(yx))
            outannual[i, metrics[:days]] = sum(getfield(ev, :duration)[yx])
        end
    end
    outannual
end


function trendm(outannual::Array{T, 3})
    X = 1:size(outannual, 1) # the number of years
    z = length(metrics)
    sz = Base.tail(size(outannual))
    outpvalue = Matrix{T}(undef, sz)
    outcoeff = similar(outpvalue)
    outrsqd = similar(outpvalue)
    outerror_coeff = similar(outpvalue)
    outintercept = similar(outpvalue)
    for j in axes(outannual, 3)
        for i in axes(outannual, 2)
            outlg = linreg(X, outannual[:, i, j])
            outpvalue[i, j] = _pvalue(outlg)
            outcoeff[i, j] = getindex(outlg, 1)
            outintercept[i, j] = getindex(outlg, 2)
            outrsqd[i, j] = getindex(outlg, 3)
            outerror_coeff[i, j] = getindex(outlg, 4)
        end
    end
    outcoeff, outerror_coeff, outrsqd, outintercept, outpvalue
end

# Vector version Output: 
function trendv(outannual::Matrix{T})
    X = 1:size(outannual, 1)
    z = length(metrics)
    outpvalue = Vector{T}(undef, z)
    outcoeff = Vector{T}(undef, z)
    outrsqd = Vector{T}(undef, z)
    outintercept = Vector{T}(undef, z)
    outerror_coeff = Vector{T}(undef, z)
    for i in axes(outannual, 2)
        outlg = linreg(X, outannual[:, i])
        outpvalue[i] = _pvalue(outlg)
        outcoeff[i] = getindex(outlg, 1)
        outintercept[i] = getindex(outlg, 2)
        outrsqd[i] = getindex(outlg, 3)
        outerror_coeff[i] = getindex(outlg, 4)
    end
    outcoeff, outerror_coeff, outrsqd, outintercept, outpvalue
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
    return p_value_f #, p_value_t
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
            ndims(xx) == 1 ? getindex(xx, idx) : getindex(xx, :, :, idx)
        end
    end
end

"""
Return the Eventsfull of a single pixel as a matrix

"""
function mymetric(ev::Events{Vector{T}})
    # Return a vector of vectors
    vps = first(propertynames(ev), 7)
    stack([reduce(vcat, getfield(ev, t)) for t in vps])
end

"""
Return all the values of the events of a given metric as a vector.

mymetric(eve, "means")::Vector{T}
"""
function mymetric(ev::Events, metric::SymOrString)
    metric = typeof(metric) == String ? Symbol(metric) : metric
    fd = first(propertynames(ev), 7)
    if metric in fd
        return reduce(vcat, getfield(ev, metric))
    else
        @info "Perhaps you mean one of ", keys(metrics)
        throw(KeyError(metric))
    end
end
 
# NOTE: We should also return the start and end dates

"""
Return the Events be used as a table (or dataframe).

mymetric(evt, indices, startends)::Matrix
"""
function mymetric(ev::Events{Vector{Vector{T}}}, indices::Tuple, evst::Tuple)
    vps = first(propertynames(ev), 7)
    evs = [reduce(vcat, getfield(ev, t)) for t in vps]
    cix = getindex(indices, 1)
    # the number of events per pixel
    drs = length.(getfield(ev, 1))
    ix = TI[]
    iy = TI[]
    for (q, s) in zip(drs, cix)
        append!(ix, repeat([Tuple(s)[1]], q))
        append!(iy, repeat([Tuple(s)[2]], q))
    end
    stdates = [mdate[i] for i in Iterators.flatten(evst[1]) |> collect]
    endates = [mdate[i] for i in Iterators.flatten(evst[2]) |> collect]
    @assert length(ix) == length(iy) == sum(drs)
    stack([evs..., stdates, endates, ix, iy])
end


"""
Return the climatology mean and threshold as input array.

ma = mymetric(mm, arg, indices)
"""
function mymetric(mm::MExtreme, field::SymOrString, indices)
    # return clim, thresh, annual, tempanom, threshanom or exceedance as 3-D arrays
    # mapping the argument to the fieldnames of the mm
    field = typeof(field) == String ? Symbol(field) : field
    field in fieldnames(typeof(mm)) || throw(error(field, " not a valid option. Try any of `:clim`, `:thresh`, `:exceeds`")) 
    CIx, x, y = indices
    # CIx, nCIx, x, y = indices
    gg = getfield(mm, field)
    outs = eltype(gg) == Bool ? falses(x, y, size(gg, 1)) : fill(T(NaN), x, y, size(gg, 1)) # Array{T, 3}(undef, x, y, size(gg, 1))
    # outs = fill(T(NaN), x, y, size(gg, 1))
    outs[CIx, :] = permutedims(gg)
    # outs[nCIx, :] .= NaN
    return outs
end

function mymetricm(mm::MHCMetricsm,
                   field::SymOrString,
                   metric::Union{SymOrString, Tuple{Vararg{SymOrString}}},
                   indices)

    field = typeof(field) == String ? Symbol(field) : field
    field in propertynames(mm) || throw(KeyError(field))
    gg = getfield(mm, field)
    # CIx, nCIx, x, y = indices
    CIx, x, y = indices
        # first for a single metric
    if metric isa SymOrString
        metric = typeof(metric) == String ? Symbol(metric) : metric
        metric in keys(metrics) || throw(KeyError(metric))
        idx = getindex(metrics, metric)
        # :annuals is 1-D larger than the others
        if field == :annuals
            gg = getindex(gg, :, :, idx)
            outs = fill(T(NaN), x, y, size(gg, 1))
            # outs = Array{T, 3}(undef, x, y, size(gg, 1))
            outs[CIx, :] = permutedims(gg)
            # outs[nCIx, :] .= NaN
            return outs
        end
        # outs = Matrix{T}(undef, x, y)
        outs = fill(T(NaN), x, y)
        outs[CIx] = gg[:, idx]
        # outs[nCIx] .= NaN
        return outs
    else # metric is a tuple
        M = length(metric)
        nmet = typeof(metric) == NTuple{M, Symbol} ? metric : Symbol.(metric)
        idx = TI[]
        for mt in nmet
            mt in keys(metrics) || throw(KeyError(mt))
            push!(idx, getindex(metrics, mt))
        end
        @assert M == length(idx)
        if field == :annuals
            # outs = ntuple(_ -> Array{T, 3}(undef, x, y, size(gg, 1)), M)
            outs = ntuple(_ -> fill(T(NaN), x, y, size(gg, 1)), M)
            for i in 1:M
                og = getindex(gg, :, :, getindex(idx, i))
                outs[i][CIx, :] = permutedims(og)
                # outs[i][nCIx, :] .= NaN
            end
            return outs
        end
        # outs = ntuple(_ -> Matrix{T}(undef, x, y), M)
        outs = ntuple(_ -> fill(T(NaN), x, y), M)
        for i in 1:M
            outs[i][CIx] = getindex(gg, :, getindex(idx, i))
            # outs[i][nCIx] .= NaN
        end
        return outs
    end
end


"""
Return the Events and the labels (starts and end positions) of the event.
"""

function meventsm(ms::MExtreme, mindur::Integer, maxgap::Integer)
    mindur, maxgap = convert(Vector{TI}, [mindur, maxgap])
    evst = mylabel(ms, mindur, maxgap)
    ev = anomsam(ms, evst)
    return ev, evst
end


"""
Return the computed metrics - means, annuals and linear regression outputs as MHCMetrics.
"""

function mmetricsm(ev::Events, evst, mdate, indices)
    mm = meanmetricsm(ev, mdate)
    am = annualmetricsm(ev, mdate, evst)
    tr = trendm(am)
    MHCMetricsm{typeof(am), typeof(mm)}(am, mm, tr...)
end

#= STEPS
    # 1. calculate the MHW/MCS
    ms, indices = mextreme(sst, sstdate::StepRange, mhwdate::StepRange, clmdate::StepRange, event=:mhw; window=5, smoothwindow=31, threshold=nothing)
    # 2. label + events
    ev, lb = mevents(ms, mindur, maxgap, indices)
    # 3. pixel, annual and trend
    mht = mmetrics(ev, lb, indices, mdate)
=#
export seamask, timeindices, _subtemp, MHW, MCS, EventHW, EventCS, anomsa, Events, trendm
