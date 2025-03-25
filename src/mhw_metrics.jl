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
eventmetrics(
    msst::Union{MCTemp,MHTemp}, mstartsxs, mendsxs
) = _eventmetrics(msst.temp, msst.clima, msst.thresh, msst.lyday, msst.anomfn, msst.argfn, msst.dates, mstartsxs, mendsxs)

#=

function matevents(msms::MarineHW, msst::MarineHeatWave, mexc::BitMatrix, mstartsxs, mendsxs)
    mni, cmi, mxi, drt, cats, rons, rdcs, vri, sdt, edt, cls, rws = ntuple(_ -> Vector(undef, length(mstartsxs)), 12)
    fullyears = unique(year.(msst.dates))
    mask = msst.mask
    x, y = size(msms.temp)
    climout, threshout = ntuple(_ -> similar(msms.temp, x, y, 366), 2)
    msms.exceed[mask, :] = mexc
    for i in eachindex(mstartsxs)
        evmeta = _eventmetrics(eachcol(msst.temp)[i], eachcol(msst.clima)[i], eachcol(msst.thresh)[i], msst.lyday, msst.anomfn, msst.argfn, msst.dates, mstartsxs[i], mendsxs[i])
        meanmets, evmets = meanmetrics(evmeta[1], evmeta[2], evmeta[3], fullyears, msst.anomfn)
        cats[i] = evmeta[4]
        anmets = annualmetrics(evmeta[1], evmeta[2], evmeta[3], evmeta[5], evmeta[6], fullyears)
        msms.temp[mask[i], :] = evmeta[7]
        msms.category[mask[i], :] = evmeta[8]
        climout[mask[i], :] = eachcol(msst.clima)[i]
        threshout[mask[i], :] = eachcol(msst.thresh)[i]
        rons[i], rdcs[i], sdt[i], edt[i] = evmeta[2], evmeta[3], evmeta[5], evmeta[6]
        rws[i] = repeat([mask[i][1]], length(evmeta[2]))
        cls[i] = repeat([mask[i][2]], length(evmeta[2]))
        mni[i], cmi[i], mxi[i], drt[i], vri[i] = evmets
        trmets = trends(anmets)

        for m in keys(msms.annuals)
            msms.annuals[m][mask[i], :] = anmets[m]
            msms.means[m][mask[i]] = meanmets[m]
            msms.trends[m][mask[i]] = trmets.trends[m]
            msms.pvalues[m][mask[i]] = trmets.pvalues[m]
            msms.pmetrics[m][mask[i]] = trmets.pmetrics[m]
        end

    end

    edf = (MeanInt=reduce(vcat, mni), CumInt=reduce(vcat, cmi), MaxInt=reduce(vcat, mxi), Duration=reduce(vcat, drt), Category=reduce(vcat, cats), ROnset=reduce(vcat, rons), RDecline=reduce(vcat, rdcs), VarInt=reduce(vcat, vri), StartDate=reduce(vcat, sdt), EndDate=reduce(vcat, edt), Xind=reduce(vcat, rws), Yind=reduce(vcat, cls))

    return msms, edf, climout, threshout
end


# TODO: a more appropriate name for this function.
- is it more safe memory wise to save all the vector space in memory or should we automatically chunk it? How do we achieve this?

=#

function meanmetrics(evanom, rons, rdec, fullyears, anomfn)
    lfy = length(fullyears)
    # Per Event Metrics
    meanint = mean.(evanom)
    cumint = sum.(evanom)
    maxint = anomfn.(evanom)
    duration = length.(evanom)
    varint = std.(evanom)
    evmets = (; meanint, cumint, maxint, duration, varint)

    # Overall Mean Metrics
    meanint = mean(Iterators.flatten(evanom))
    cumint = mean(cumint)
    maxint = anomfn(Iterators.flatten(evanom))
    ronset = mean(rons)
    rdecline = mean(rdec)
    frequency = length(duration) / lfy
    days = sum(duration) / lfy
    duration = mean(length.(evanom))

    meanmets = (; meanint, cumint, maxint, ronset, rdecline, duration, days, frequency)

    return meanmets, evmets
end

function annualmetrics(evanom, ronset, rdecline, stdate, endate, fullyears, anomfn)
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
        annuals[6][ey] = anomfn(evanomf[yearixs])   # maxint
        annuals[7][ey] = length(evanomf[yearixs])    # days
        annuals[8][ey] = length(yrix)                # frequency
    end
    return NamedTuple{metrics}(annuals)
end

"""
helper function to compute the trend. It takes in a specific annual metric and calculates the linear regression on it.

- metric : a vector that potentially has NaN values, its length should be the number of years in the annual metrics.
- 
"""
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

"""
this function is supposed to wrap the four processes internal of obtaining the values of an event.

- Essentially, we'll wrap this up for a vector and a matrix.
- So the vector version will loop through the starts and ends keeping the sst, clim and threshold constant.
- Perhaps we can use the lyd to create a vector for the clim and thresh to be as long as the sst.
- it returns the key metrics that will be incorporated in the output without presenting the final outputs themselves.
- metrics returned are: vector - anomaly, scalars - category, rate of onset and rate of decline
"""

# NOTE: A matrix version of the inner functions is not feasible because each pixel represented by a column vector has different starts and ends.

function newmets(sst, clim, thrs, lyd, mst, mse, anomfn, argfn)
    eanoms = _anoms(sst, clim, thrs, lyd, mst, mse)
    cats = _categorys(eanoms)
    rons = ronset(eanoms, mst, anomfn, argfn)
    rdcs = rdecline(eanoms, mse, anomfn, argfn)
    return eanoms.anom, cats, rons, rdcs
end

"""
We return the mhw anomalies and categories in the `mhwout` and `catout` in the same format as the initial sst data. So it has the same dimensions with non-mhw events as `0` and mhw events having `non-zero` values.
"""
function newmetsvector(msstobject::Vector, mstarts, mends)

    mhwout, catout = (zeros(size(msstobject.temp)) for _ in 1:2)
    ll = length(mstarts)
    cats, rons, rdcs = ntuple(_ -> Vector{eltype(msstobject.temp)}(undef, ll), 3)
    mdurations = mends - mstarts .+ 1
    # evanoms = [Vector{eltype(msstobject.temp)}(undef, mstarts[l]) for l in eachindex(mstarts)] 
    evanoms = [Vector{eltype(msstobject.temp)}(undef, l) for l in mdurations]

    for (i, (mst, mse)) in enumerate(zip(mstarts, mends))
        eanoms = newmets(msstobject.temp, msstobject.clim, msstobject.thresh, msstobject.lyd, mst, mse, msstobject.anomfn, msstobject.argfn)
        cats[i], rons[i], rdcs[i] = eanoms.cats, eanoms.rons, eanoms.rdcs
        evanoms[i] = eanoms.anom
        mhwout[mst:mse] = eanoms.anom
        catout[mst:mse] .= eanoms.cats
    end
    return mhwout, catout, evanoms, rons, rdcs, cats
    # NOTE: are we not returning the `eanoms` for further work in the metrics?
    # Answer: we are in the evanoms variable.
end

# bear in mind that in the matrix version, the starts and ends are vectors of vectors, meaning that we are going into it at two levels before performing the operation.
# can we use the vector version, probably not as it targets the object itself not the array under.

function newmetsmatrix(msstobject::Matrix, mstarts, mends)

    mhwout, catout = (zeros(size(msstobject.temp)) for _ in 1:2)
    ll = length.(mstarts)
    cats, rons, rdcs = ntuple(_ -> [Vector{eltype(msstobject.temp)}(undef, l) for l in ll], 3)
    mdurations = mends - mstarts .+ 1
    # here mdurations is a vector of vectors
    evanoms = [[Vector{eltype(msstobject.temp)}(undef, l) for l in mdurations[x]] for x in eachindex(mdurations)]
    # @assert length.(evanoms) == ll
    # @assert length.(evanoms[1]) == mstarts[1]
    # evanoms = Vector{Vector{eltype(sst)}}(undef, ll)
    # not yet completed
    for j in axes(msstobject.temp, 2)
        for (i, (mst, mse)) in enumerate(zip(mstarts[j], mends[j]))
            eanoms = newmets(msstobject.temp[:, j], msstobject.clim[:, j], msstobject.thresh[:, j], msstobject.lyd, mst, mse, msstobject.anomfn, msstobject.argfn)
            cats[j][i], rons[j][i], rdcs[j][i] = eanoms.cats, eanoms.rons, eanoms.rdcs
            evanoms[j][i] = eanoms.anom
            mhwout[mst:mse, j] = eanoms.anom
            catout[mst:mse, j] .= eanoms.cats
        end
    end
    # the mhwout and catout should be transposed
    return mhwout', catout', evanoms, rons, rdcs, cats
end



#= For the mean and annual metrics

evyears = vcat(map((x, y) -> x:y, stdate, endate)...) .|> year

NOTE: for the meanmetrics we can pass the yearobject and the anomfn. so no need to pass the msstobject.

function meanmetrics(evanom, rons, rdec, yearobject, anomfn)
    # fullyears = unique(year.(msstobject.dates))
    lfy = length(yearobject.fullyears)
    # Per Event Metrics
    meanint = mean.(evanom)
    cumint = sum.(evanom)
    maxint = anomfn.(evanom)
    duration = length.(evanom)
    varint = std.(evanom)
    evmets = (; meanint, cumint, maxint, duration, varint)

    # Overall Mean Metrics
    meanint = mean(Iterators.flatten(evanom))
    cumint = mean(cumint)
    maxint = anomfn(Iterators.flatten(evanom))
    ronset = mean(rons)
    rdecline = mean(rdec)
    frequency = length(duration) / lfy
    days = sum(duration) / lfy
    duration = mean(length.(evanom))

    meanmets = (; meanint, cumint, maxint, ronset, rdecline, duration, days, frequency)

    return meanmets, evmets
end



NOTE: In the function we are only using the full years, event years and start and end years so essentially that's what we should return ⌣ No need for the start/end dates themselves.
What we actually need is a function that takes the object, the start and end indices and returns the following:
- event years
- full years
- start years
- end years

# can broadcast either outside or inside at year()

this would then return the start/end year as either scalar or vector,
the eventyear will always be a vector regardless of input, fullyears is always a vector. Can return as a named tuple so we sub what we want by name.

- we could call getyears inside the annual metrics but we also need it for mean metrics so it makes more sense to compute it before the two metrics functions.
- good place to call it would be at the vector level since what we're passing to mean and annual metrics are vectors.

=#
function getyears(msstobject, mst, mse)

    # this would work for a scalar or vector of scalars
    startyear = year.(msstobject.dates[mst])
    endyear = year.(msstobject.dates[mse])
    fullyears = unique(year.(msstobject.dates)) # could be ignored though

    # base function below
    gevyear(dateobject, mst, mse) = year.(dateobject[mst:mse])

    # generalises to scalars and vectors
    # eventyears = map((x, y) -> gevyear(msstobject.dates, x, y), mst, mse)

    # flattened using vcat and splat
    eventyears = vcat(map((x, y) -> gevyear(msstobject.dates, x, y), mst, mse)...)

    return (; startyear, endyear, eventyears, fullyears)
end


function _ntrendv(metric)
    x = 1:length(metric)
    lr = linreg(x, metric)
    return lr
end

function ntrends(annualmetrics)
    trds = (:trends, :pvalues, :pmetrics)
    m = length(annualmetrics)
    tss = ntuple(_ -> ones(m), 3)
    for (m, metric) in enumerate(annualmetrics)
        tss[1][m], tss[2][m], tss[3][m] = _ntrendv(metric)
    end
    tss = (NamedTuple{keys(annualmetrics)}(x) for x in tss)
    return (; zip(trds, tss)...)
end

function meanmetsv(evanom, rons, rdec, yearobject, anomfn)
    # Default metrics
    lfy = length(yearobject.fullyears)
    # Overall Mean Metrics
    meanint = mean(Iterators.flatten(evanom))
    cumint = mean(sum.(evanom))
    maxint = anomfn(Iterators.flatten(evanom))
    ronset = mean(rons)
    rdecline = mean(rdec)
    frequency = length(length.(evanom)) / lfy # length.(evanom) can be made a temp variable since it is used thrice. Same as Iter.flatten(evanom) above?
    days = sum(length.(evanom)) / lfy
    duration = mean(length.(evanom))
    return  (; meanint, cumint, maxint, ronset, rdecline, duration, days, frequency)
end

function meanmetsm(evanom, rons, rdec, yearobject, anomfn)
    # Default metrics
    meansmets = map((x, y, z) -> meanmetsv(x, y, z, yearobject, anomfn), evanom, rons, rdec)
    @assert length(meansmets) == length(evanom)
    return meansmets # vector of named tuples. length 
end

# Example usage
# meansmatrixout = Matrix(undef, size(sst, 1), size(sst, 2))
# meansmatrixout[maskci] = meansmets
# getindex.(meansmets, 1) # == vec(meanint)

function eventmets(evanom, anomfn)
    
    # Per Event Metrics
    meanint = mean.(evanom)
    cumint = sum.(evanom)
    maxint = anomfn.(evanom)
    duration = length.(evanom)
    varint = std.(evanom)
    (; meanint, cumint, maxint, duration, varint)
end

function annualmetricsv(evanom, ronset, rdecline, yearobject, anomfn)
    lfy = length(yearobject.fullyears)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maxint, :days, :frequency)
    evanomf = collect(Iterators.flatten(evanom))
    ametrics = NamedTuple{metrics}([zeros(lfy) for _ in metrics])
    for (ey, yr) in enumerate(yearobject.fullyears)
        yearixs = findall(isequal(yr), yearobject.eventyears)
        yrix = [i for (i, (a, b)) in enumerate(zip(yearobject.startyear, yearobject.endyear)) if a == yr || b == yr]
        if isempty(yearixs)
            continue
        end
        ametrics.meanint[ey] = mean(evanomf[yearixs])      # meanint
        ametrics.cumint[ey] = mean(sum.(evanom[yrix]))    # cumint
        ametrics.ronset[ey] = mean(ronset[yrix])
        ametrics.rdecline[ey] = mean(rdecline[yrix])
        ametrics.duration[ey] = mean(length.(evanom[yrix])) # duration
        ametrics.maxint[ey] = anomfn(evanomf[yearixs])   # maxint
        ametrics.days[ey] = length(evanomf[yearixs])    # days
        ametrics.frequency[ey] = length(yrix)                # frequency
    end
    return ametrics
end

function annualmetricsv(evanom, ronset, rdecline, yearobject, anomfn)
    lfy = length(yearobject.fullyears)
    metrics = (:meanint, :cumint, :onset, :decline, :duration, :maxint, :days, :frequency)
    evanomf = collect(Iterators.flatten(evanom))
    ametrics = NamedTuple{metrics}([zeros(lfy) for _ in metrics])
    for (ey, yr) in enumerate(yearobject.fullyears)
        yearixs = findall(isequal(yr), yearobject.eventyears)
        yrix = [i for (i, (a, b)) in enumerate(zip(yearobject.startyear, yearobject.endyear)) if a == yr || b == yr]
        if isempty(yearixs)
            continue
        end

        if :meanint in metrics
            ametrics.meanint[ey] = mean(evanomf[yearixs])      # meanint
        end
        if :cumint in metrics
            ametrics.cumint[ey] = mean(sum.(evanom[yrix]))    # cumint
        end
        if :onset in metrics
            ametrics.onset[ey] = mean(ronset[yrix])
        end
        if :decline in metrics
            ametrics.decline[ey] = mean(rdecline[yrix])
        end
        if :duration in metrics
            ametrics.duration[ey] = mean(length.(evanom[yrix])) # duration
        end
        if :maxint in metrics # or minint for coldspells
            ametrics.maxint[ey] = anomfn(evanomf[yearixs])   # maxint
        end
        if :days in metrics
            ametrics.days[ey] = length(evanomf[yearixs])    # days
        end
        if :frequency in metrics
            ametrics.frequency[ey] = length(yrix)                # frequency
        end
    end
    return ametrics
end


function annualmetricsm(evanom, ronset, rdecline, yearobject, anomfn)
    
    annmets = map((x, y, z) -> annualmetricsv(x, y, z, yearobject, anomfn), evanom, ronset, rdecline)
end
