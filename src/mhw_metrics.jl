"""
    This function calculates the anomalies for the temperature array and the clim/threshold.
- anom: the anomaly
- thsd: the threshold difference
- anbx (index): for the rate of onset
- anfx (index): for the rate of decline

- anbf (value): for the rate of onset
- anft (value): for the rate of decline
"""

function _anoms(sst, clim, thsh, lyd, mst, mse)
    anom = sst[mst:mse] - clim[lyd[mst:mse]]
    thsd = thsh[lyd[mst:mse]] - clim[lyd[mst:mse]]
    lnxt = lastindex(sst)
    anbx, anfx = max(1, mst - 1), min(lnxt, mse + 1)
    anbf = sst[anbx] - clim[lyd[anbx]]
    anft = sst[anfx] - clim[lyd[anfx]]
    return (; anom, anbf, anft, lnxt, thsd)
end

# TODO: DOUBLE CHECK THIS FOR THE COLD SPELL
# INFO: This is not where the issue lies.

_categorys(anoms) = min(4, maximum(fld.(anoms.anom, anoms.thsd)))

function ronset(anoms, mst, anomfn, argfn)
    anom, anbf = anoms.anom, anoms.anbf
    nmx, fan, ngx = anomfn(anom), first(anom), argfn(anom)
    lnmx, snmx = nmx - mean((fan, anbf)), nmx - fan
    ron = mst > 1 ? /(lnmx, (ngx + 0.5)) : /(snmx, ngx)
    return ron
end

function rdecline(anoms, mse, anomfn, argfn)
    anft, lnx, nmx, lan, ngx, lnt = (anoms.anft, anoms.lnxt, anomfn(anoms.anom), last(anoms.anom), argfn(anoms.anom), length(anoms.anom))
    lnmx, snmx = (nmx - mean((lan, anft)), nmx - lan)
    rdc = mse < lnx ? /(lnmx, (lnt - ngx + 0.5)) : ngx == lnx ? /(snmx, 1.0) : /(snmx, (lnt - ngx))
    return rdc
end

"""
`_eventmetrics` is the base function for calculating the events in each pixel. Since each pixel has potentially many events, we create vectors to hold each output. It returns:
- category
- rate of onset
- rate of decline
- event anomaly
- start date
- end date
- the marine heatwave/coldspell vector (mhwout)
- the category vector (catout)


"""
function _eventmetrics(sst, clim, thrs, lyd, anomfn, argfn, evdate, mstarts, mends)
    l = length(mstarts)
    cats, rons, rdcs = ntuple(_ -> Vector{eltype(sst)}(undef, l), 3)
    evanom = Vector{Vector{eltype(sst)}}(undef, l)
    stdate, endate = ntuple(_ -> Vector{Date}(undef, l), 2)
    mhwout, catout = (zeros(size(sst)) for _ in 1:2)
    for (m, (mst, mse)) in enumerate(zip(mstarts, mends))
        eanoms = _anoms(sst, clim, thrs, lyd, mst, mse)
        evanom[m] = eanoms.anom
        cats[m] = _categorys(eanoms)
        rons[m] = ronset(eanoms, mst, anomfn, argfn)
        rdcs[m] = rdecline(eanoms, mse, anomfn, argfn)
        stdate[m] = evdate[mst]
        endate[m] = evdate[mse]
        mhwout[mst:mse] = eanoms.anom
        catout[mst:mse] .= cats[m]
    end
    return evanom, rons, rdcs, cats, stdate, endate, mhwout, catout
end


"""
    Version for single vector which could be part of a loop if we do the values of one pixel all the way to the end.

    eventmetrics(sst, climthresh, startendsindices) -> Vector{Vector{T}}[events, mhwout, mhwcat]
"""
eventmetrics(msst::Union{MCTemp,MHTemp}, mstartsxs, mendsxs) = _eventmetrics(msst.temp, msst.clima, msst.thresh, msst.lyday, msst.anomfn, msst.argfn, msst.dates, mstartsxs, mendsxs)

function meanmetrics(evanom, rons, rdec, fullyears)
    lfy = length(fullyears)
    # Per Event Metrics
    meanint = mean.(evanom)
    cumint = sum.(evanom)
    maxint = maximum.(evanom)
    duration = length.(evanom)
    varint = std.(evanom)
    evmets = (; meanint, cumint, maxint, duration, varint)

    # Overall Mean Metrics
    meanint = mean(Iterators.flatten(evanom))
    cumint = mean(cumint)
    maxint = maximum(Iterators.flatten(evanom))
    ronset = mean(rons)
    rdecline = mean(rdec)
    frequency = length(duration) / lfy
    days = sum(duration) / lfy
    duration = mean(length.(evanom))

    meanmets = (; meanint, cumint, maxint,
        ronset, rdecline, duration,
        days, frequency)

    return meanmets, evmets
end

function annualmetrics(evanom, ronset, rdecline, stdate, endate, fullyears)
    lfy = length(fullyears)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maxint, :days, :frequency)
    evyears = collect(Iterators.flatten(map((x, y) -> x:y, stdate, endate))) .|> year
    styr, enyr = year.(stdate), year.(endate)
    evanomf = collect(Iterators.flatten(evanom))
    annuals = ntuple(_ -> zeros(lfy), length(metrics))
    for (ey, yr) in enumerate(fullyears)
        yearixs = findall(isequal(yr), evyears)
        yrix = [i for (i, (a, b)) in enumerate(zip(styr, enyr)) if a == yr || b == yr]
        if isempty(yearixs)
            continue
        end
        annuals[1][ey] = mean(evanomf[yearixs])      # meanint
        annuals[2][ey] = mean(sum.(evanom[yrix]))    # cumint
        annuals[3][ey] = mean(ronset[yrix])
        annuals[4][ey] = mean(rdecline[yrix])
        annuals[5][ey] = mean(length.(evanom[yrix])) # duration
        annuals[6][ey] = maximum(evanomf[yearixs])   # maxint
        annuals[7][ey] = length(evanomf[yearixs])    # days
        annuals[8][ey] = length(yrix)                # frequency
    end
    return NamedTuple{metrics}(annuals)
end

function __trendv(metric)
    trend, pval, pmet = ntuple(_ -> ones(1), 3)
    m = size(metric, 1)
    X = 1:m
    y = replace(metric, NaN => missing)
    data = DataFrame(X=X, Y=y)
    lr = lm(@formula(Y ~ X), data)
    trend[1] = coef(lr)[2]
    pval[1] = coeftable(lr).cols[4][2]
    pmet[1] = Float64(prod(confint(lr)[2, :]) ≥ 0)
    return trend[1], pval[1], pmet[1]
end

function trends(annualmetrics)
    trds = (:trends, :pvalues, :pmetrics)
    m = length(annualmetrics)
    tss = ntuple(_ -> ones(m), 3)
    for (m, metric) in enumerate(annualmetrics)
        tss[1][m], tss[2][m], tss[3][m] = __trendv(metric)
    end
    tss = (NamedTuple{keys(annualmetrics)}(x) for x in tss)
    return (; zip(trds, tss)...)
end


