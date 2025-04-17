"""
`_eventmetrics` is the base function for calculating the events in each pixel. Since each pixel has potentially many events, we create vectors to hold each output. It returns:
- category
using Dates: max_width
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
        anmets = annualmetrics(evmeta[1], evmeta[2], evmeta[3], evmeta[5], evmeta[6], fullyears, msst.anomfn)
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
# - is it more safe memory wise to save all the vector space in memory or should we automatically chunk it? How do we achieve this?



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

function _newmets(sst, clim, thrs, lyd, mst, mse, anomfn, argfn)
    eanoms = _anoms(sst, clim, thrs, lyd, mst, mse)
    cats = _categorys(eanoms)
    rons = ronset(eanoms, mst, anomfn, argfn)
    rdcs = rdecline(eanoms, mse, anomfn, argfn)
    return eanoms.anom, cats, rons, rdcs
end

function _anomsnew(sst, clim, thsh, mst, mse)
    anom = sst[mst:mse] - clim[mst:mse]
    thsd = thsh[mst:mse] - clim[mst:mse]
    onsan = sst[max(1, mst - 1)] - clim[max(1, mst - 1)]
    decan = sst[min(lastindex(sst), mse + 1)] - clim[min(lastindex(sst), mse + 1)]
    return anom, onsan, decan, thsd
end

anoms(mhw::MHTemp, mst, mse) = WVH(_anomsnew(mhw.temp, mhw.clima, mhw.thresh, mst, mse))

anoms(mhw::MCTemp, mst, mse) = WVC(_anomsnew(mhw.temp, mhw.clima, mhw.thresh, mst, mse))

function _newmetsnew(mhw, mst, mse)
    eanoms = anoms(mhw, mst, mse)
    cats = _categorys(eanoms)
    rons = ronsetnew(eanoms, mst)
    rdcs = rdeclinenew(eanoms, mse, lnx)
    return eanoms.anom, cats, rons, rdcs
end
"""
We return the mhw anomalies and categories in the `mhwout` and `catout` in the same format as the initial sst data. So it has the same dimensions with non-mhw events as `0` and mhw events having `non-zero` values.
"""
function newmets(outarray::AbstractArray, msstobject, mstarts::Vector{Int}, mends::Vector{Int})

    mhwout, catout = (zeros(size(msstobject.temp)) for _ in 1:2)
    ll = length(mstarts)
    cats, rons, rdcs = ntuple(_ -> Vector{eltype(msstobject.temp)}(undef, ll), 3)
    mdurations = mends - mstarts .+ 1
    evanoms = [Vector{eltype(msstobject.temp)}(undef, l) for l in mdurations]

    for (i, (mst, mse)) in enumerate(zip(mstarts, mends))
        eanoms = _newmets(msstobject.temp, msstobject.clim, msstobject.thresh, msstobject.lyd, mst, mse, msstobject.anomfn, msstobject.argfn)
        cats[i], rons[i], rdcs[i] = eanoms.cats, eanoms.rons, eanoms.rdcs
        evanoms[i] = eanoms.anom
        outarray[1][mst:mse] = eanoms.anom
        outarray[2][mst:mse] .= eanoms.cats
        # mhwout[mst:mse] = eanoms.anom
        # catout[mst:mse] .= eanoms.cats
    end
    # return mhwout, catout, evanoms, rons, rdcs, cats
    return outarray, evanoms, rons, rdcs, cats
end

# bear in mind that in the matrix version, the starts and ends are vectors of vectors, meaning that we are going into it at two levels before performing the operation.

function newmets(msstobject::Matrix, mstarts, mends)
    # TODO: the mhwout and catout should be the final one so we don't need to write a different function to return it in its output state
    mhwout, catout = (zeros(size(msstobject.temp)) for _ in 1:2)
    ll = length.(mstarts)
    cats, rons, rdcs = ntuple(_ -> [Vector{eltype(msstobject.temp)}(undef, l) for l in ll], 3)
    mdurations = mends - mstarts .+ 1
    # here mdurations is a vector of vectors
    evanoms = [[Vector{eltype(msstobject.temp)}(undef, l) for l in mdurations[x]] for x in eachindex(mdurations)]
    # @assert length.(evanoms) == ll
    # @assert length.(evanoms[1]) == mstarts[1]
    for j in axes(msstobject.temp, 2)
        for (i, (mst, mse)) in enumerate(zip(mstarts[j], mends[j]))
            eanoms = _newmets(msstobject.temp[:, j], msstobject.clim[:, j], msstobject.thresh[:, j], msstobject.lyd, mst, mse, msstobject.anomfn, msstobject.argfn)
            cats[j][i], rons[j][i], rdcs[j][i] = eanoms.cats, eanoms.rons, eanoms.rdcs
            evanoms[j][i] = eanoms.anom
            mhwout[mst:mse, j] = eanoms.anom
            catout[mst:mse, j] .= eanoms.cats
        end
    end
    # the mhwout and catout should be transposed
    return mhwout', catout', evanoms, rons, rdcs, cats
end

# complementary function that returns the output array from the previous method
# function mhwcatout(outputarray, first(newmetsmatrix, 2))
# run a for loop
# outputarray
# end

# Version that takes in the output array
function newmets(outputarray::AbstractArray, msstobject::Matrix, mstarts, mends)
    # TODO: the mhwout and catout should be the final one so we don't need to write a different function to return it in its output state

    # mhwout, catout = (zeros(size(msstobject.temp)) for _ in 1:2)
    ll = length.(mstarts)
    cats, rons, rdcs = ntuple(_ -> [Vector{eltype(msstobject.temp)}(undef, l) for l in ll], 3)
    mdurations = mends - mstarts .+ 1
    # here mdurations is a vector of vectors
    evanoms = [[Vector{eltype(msstobject.temp)}(undef, l) for l in mdurations[x]] for x in eachindex(mdurations)]
    # @assert length.(evanoms) == ll
    # @assert length.(evanoms[1]) == mstarts[1]
    for j in axes(msstobject.temp, 2)
        for (i, (mst, mse)) in enumerate(zip(mstarts[j], mends[j]))
            eanoms = _newmets(msstobject.temp[:, j], msstobject.clim[:, j], msstobject.thresh[:, j], msstobject.lyd, mst, mse, msstobject.anomfn, msstobject.argfn)
            cats[j][i], rons[j][i], rdcs[j][i] = eanoms.cats, eanoms.rons, eanoms.rdcs
            evanoms[j][i] = eanoms.anom
            # mhwout[mst:mse, j] = eanoms.anom
            # catout[mst:mse, j] .= eanoms.cats
            outputarray[1][msstobject.mask[j], mst:mse] = eanoms.anom
            outputarray[2][msstobject.mask[j], mst:mse] .= eanoms.cats
        end
    end
    return outputarray, evanoms, rons, rdcs, cats
end
# what is the fill value for the contents of the outputarrays where there is no mhw? NaN or 0?


#= For the mean and annual metrics



NOTE: In the function we are only using the full years, event years and start and end years so essentially that's what we should return \:smiley: No need for the start/end dates themselves.
What we actually need is a function that takes the object, the start and end indices and returns the following:
- event years
- full years
- start years
- end years


this would then return the start/end year as either scalar or vector,
the eventyear will always be a vector regardless of input, fullyears is always a vector. Can return as a named tuple so we sub what we want by name.

- we could call getyears inside the annual metrics but we also need it for mean metrics so it makes more sense to compute it before the two metrics functions.
- good place to call it would be at the vector level since what we're passing to mean and annual metrics are vectors.
=#
# TODO: ensure that getyears works well with the annual and mean metrics: DONE 

function getyears(dateobject, mst::Vector{Int}, mse::Vector{Int})
    # this would work for a scalar or vector of scalars 
    startyear = year.(dateobject[mst])
    endyear = year.(dateobject[mse])

    # base function below
    gevyear(dateobject, mst, mse) = year.(dateobject[mst:mse])

    # flattened using vcat and splat
    eventyears = vcat(map((x, y) -> gevyear(dateobject, x, y), mst, mse)...)

    return (; startyear, endyear, eventyears)
end

getyears(dateobject, mst::Vector{Vector{Int}}, mse::Vector{Vector{Int}}) = map((x, y) -> getyears(dateobject, x, y), mst, mse)

#=
test getyears:
tdate = Date(2020):Date(2025);
ym = MarineHeatwaves.getyears(tdate, [[4, 12, 29],[4, 10, 19]] , [[9, 20, 36], [10, 19, 24]])
getindex.(ym, :startyear) [2020, 2020, 2020]
 [2020, 2020, 2020]
getindex.(ym, :eventyears)
ym[1]
(startyear = [2020, 2020, 2020], endyear = [2020, 2020, 2020], eventyears = [2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020  …  2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020])
=#


# TODO: modify the linreg function to also compute the pvalue
# What does linreg return : a, b, r2, sigma_a, sigma_b, sigma_e

function _meanmets(evanom::Vector{AbstractFloat}, rons, rdec, fullyears, anomfn)
    lfy = length(fullyears)
    # Overall Mean Metrics
    meanint = mean(Iterators.flatten(evanom))
    cumint = mean(sum.(evanom))
    maxint = anomfn(Iterators.flatten(evanom))
    ronset = mean(rons)
    rdecline = mean(rdec)
    frequency = length(length.(evanom)) / lfy # length.(evanom) can be made a temp variable since it is used thrice. Same as Iter.flatten(evanom) above?
    days = sum(length.(evanom)) / lfy
    duration = mean(length.(evanom))
    return meanint, cumint, maxint, ronset, rdecline, duration, days, frequency
end

meanmets(evanom::Vector{AbstractFloat}, rons, rdcs, fyears, anomfn) = _meanmets(evanom, rons, rdcs, fyears, anomfn)

function meanmets(evanom::Vector{Vector{AbstractFloat}}, rons, rdec, fullyears, anomfn)
    # Default metrics
    outmatrix = [zeros(8) for _ in eachindex(evanom)]
    for (i, eva) in pairs(evanom)
        outmatrix[i] = _meanmets(eva, rons[i], rdec[i], fullyears, anomfn)
    end
    # @assert length(meansmets) == length(evanom)
    return outmatrix # vector of tuples. length 
end

# Example usage

# meansmatrixout = Matrix(undef, size(sst, 1), size(sst, 2))
# meansmatrixout[maskci]= meansmets
# getindex.(meansmets, 1) # == vec(meanint)
function meanmetsout(outarray::VecOrMat, meanmets, mask)
    if isa(meanmets, Vector{Vector{AbstractFloat}})
        for (i, mk) in pairs(mask)
            outarray[mk] = meanmets[i]
        end
    else
        outarray[:] = meanmets
    end
    return outarray
end

function eventmets(evanom::Vector, anomfn)
    # Per Event Metrics
    meanint = mean.(evanom)
    cumint = sum.(evanom)
    maxint = anomfn.(evanom)
    duration = length.(evanom)
    varint = std.(evanom)
    (; meanint, cumint, maxint, duration, varint)
end

eventmets(evanom::Vector{Vector}, anomfn) = map(ev -> eventmets(ev, anomfn), evanom)
# usage:: Matrix
# getindex.(eventmets, :meanint) etc.
# getindex.(eventmets, 1)
#
"""

This doesn't return tuples. It starts with tuples and returns a vector of vectors with dimension: length(fullyears) * length(metrics)
"""
function annualmets(evanom::Vector{AbstractFloat}, ronset, rdecline, eventyears, startyear, endyear, fullyears, anomfn)
    # lfy = length(fullyears)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maxint, :days, :frequency)
    evanomf = collect(Iterators.flatten(evanom))
    ametrics = NamedTuple(zip(metrics, ntuple(_ -> zeros(length(fullyears)), length(metrics))))
    for (ey, yr) in enumerate(fullyears)
        yearixs = findall(isequal(yr), eventyears)
        yrix = [i for (i, (a, b)) in enumerate(zip(startyear, endyear)) if a == yr || b == yr]
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
    # return ametrics
    return stack(ametrics) |> transpose |> eachcol
    # this splits the metrics into vectors of metrics for each year
end

# TODO: Should probably be able to pass the metrics as an argument somehow

# function annualmetricsv(evanom, ronset, rdecline,  eventyears, startyear, endyear, fullyears, anomfn)
#     lfy = length(fullyears)
#     metrics = (:meanint, :cumint, :onset, :decline, :duration, :maxint, :days, :frequency)
#     evanomf = collect(Iterators.flatten(evanom))
# ametrics = NamedTuple(zip(metrics, ntuple(_ -> zeros(length(fullyears)), length(metrics))))
#     for (ey, yr) in enumerate(fullyears)
#         yearixs = findall(isequal(yr), eventyears)
#         yrix = [i for (i, (a, b)) in enumerate(zip(startyear, endyear)) if a == yr || b == yr]
#
#         if isempty(yearixs)
#             continue
#         end

#         if :meanint in metrics
#             ametrics.meanint[ey] = mean(evanomf[yearixs])      # meanint
#         end
#         if :cumint in metrics
#             ametrics.cumint[ey] = mean(sum.(evanom[yrix]))    # cumint
#         end
#         if :onset in metrics
#             ametrics.onset[ey] = mean(ronset[yrix])
#         end
#         if :decline in metrics
#             ametrics.decline[ey] = mean(rdecline[yrix])
#         end
#         if :duration in metrics
#             ametrics.duration[ey] = mean(length.(evanom[yrix])) # duration
#         end
#         if :maxint in metrics # or minint for coldspells
#             ametrics.maxint[ey] = anomfn(evanomf[yearixs])   # maxint
#         end
#         if :days in metrics
#             ametrics.days[ey] = length(evanomf[yearixs])    # days
#         end
#         if :frequency in metrics
#             ametrics.frequency[ey] = length(yrix)                # frequency
#         end
#     end
#     return ametrics
# return stack(ametrics) |> transpose |> eachcol
# end

# TODO:
# Implement a version of annualmetricsv that:
# 1. doesn't use namedtuples
# 2. is a vector of vectors (length of years * number of metrics)
# 3. will be indexed by number
# 4. each metric has a fixed position
# 5. [zeros(length(metrics)) for _ in eachindex(fullyears)]
# 6. This currently works for a fixed number of metrics.
# 7. How to implement it for a variable number of metrics passed by the user??

function annualmets2(evanom, ronset, rdecline, eventyears, startyear, endyear, fullyears, anomfn)
    # lfy = length(fullyears)
    metrics = (:meanint, :cumint, :ronset, :rdecline, :duration, :maxint, :days, :frequency)
    evanomf = collect(Iterators.flatten(evanom))
    ametrics = [zeros(length(metrics)) for _ in eachindex(fullyears)]
    for (ey, yr) in enumerate(fullyears)
        yearixs = findall(isequal(yr), eventyears)
        yrix = [i for (i, (a, b)) in enumerate(zip(startyear, endyear)) if a == yr || b == yr]
        if isempty(yearixs)
            continue
        end
        ametrics[ey][1] = mean(evanomf[yearixs])      # meanint
        ametrics[ey][2] = mean(sum.(evanom[yrix]))    # cumint
        ametrics[ey][3] = mean(ronset[yrix])
        ametrics[ey][4] = mean(rdecline[yrix])
        ametrics[ey][5] = mean(length.(evanom[yrix])) # duration
        ametrics[ey][6] = anomfn(evanomf[yearixs])   # maxint
        ametrics[ey][7] = length(evanomf[yearixs])    # days
        ametrics[ey][8] = length(yrix)                # frequency
    end
    return ametrics
end

function annualmets3(evanom, ronset, rdecline, eventyears, startyear, endyear, fullyears, anomfn)
    # lfy = length(fullyears)
    metrics = (:meanint, :cumint, :onset, :decline, :duration, :maxint, :days, :frequency)
    evanomf = collect(Iterators.flatten(evanom))
    # ametrics = NamedTuple(zip(metrics, ntuple(_ -> zeros(length(fullyears)), length(metrics))))
    ametrics = [zeros(length(metrics)) for _ in eachindex(fullyears)]
    for (ey, yr) in enumerate(fullyears)
        yearixs = findall(isequal(yr), eventyears)
        yrix = [i for (i, (a, b)) in enumerate(zip(startyear, endyear)) if a == yr || b == yr]
        if isempty(yearixs)
            continue
        end
        for (m, metric) in pairs(metrics)
            if metric == :meanint #in metrics
                ametrics[ey][m] = mean(evanomf[yearixs])      # meanint
            end
            if metric == :cumint# in metrics
                ametrics[ey][m] = mean(sum.(evanom[yrix]))    # cumint
            end
            if metric == :onset #in metrics
                ametrics[ey][m] = mean(ronset[yrix])
            end
            if metric == :decline #in metrics
                ametrics[ey][m] = mean(rdecline[yrix])
            end
            if metric == :duration# in metrics
                ametrics[ey][m] = mean(length.(evanom[yrix])) # duration
            end
            if metric in (:maxint, :minint) #in metrics # minint for coldspells
                ametrics[ey][m] = anomfn(evanomf[yearixs])   # maxint
            end
            if metric == :days #in metrics
                ametrics[ey][m] = length(evanomf[yearixs])    # days
            end
            if metric == :frequency #in metrics
                ametrics[ey][m] = length(yrix)                # frequency
            end
        end
    end
    return ametrics
    # return stack(ametrics) |> transpose |> eachcol
end

function annualmets(evanom::Vector{Vector{AbstractFloat}}, ronset, rdecline, eventyears, startyear, endyear, fullyears, anomfn)
    # the fixed variables are: fullyear, anomfn
    annmets = map((an, ons, dec, evys, sy, ey) -> annualmets2(an, ons, dec, evys, sy, ey, fullyears, anomfn), evanom, ronset, rdecline, eventyears, startyear, endyear)
    return annmets
end

function trendout(outmatrixes::AbstractMatrix, annualmetrics; mask=mask)
    X = eachindex(first(annualmetrics))
    for (p, pix) in pairs(annualmetrics)
        for n in eachindex(first(pix))
            outmatrixes[1][mask][p][n],
            outmatrixes[2][mask][p][n],
            outmatrixes[3][mask][p][n] = linreg(X, getindex.(pix, n))
        end
    end
end

function trendout(outvectors::AbstractVector, annualmetrics)
    X = eachindex(first(annualmetrics))
    for n in eachindex(first(annualmetrics))
        outvectors[1][n],
        outvectors[2][n],
        outvectors[3][n] = linreg(X, getindex.(annualmetrics, n))
    end
end

function annualmetricsout(outarray::AbstractArray, annualmetrics, mask)
    if isa(annualmetrics, Vector{Vector{Vector{AbstractFloat}}})
        for (p, pix) in pairs(annualmetrics)
            outarray[mask[p], :] = pix
        end
        return outarray
    else
        return outarray[:] = annualmetrics
    end
end
# TODO: 
# INFO:
# NOTE:
# WARN:
# PERF:
# TEST:
# FIX:
# HACK:
# 
function newdetect(msst)
    mstarts, mends = mylabeling(exceed(msst), min_dur, max_gap)
    outarray, evanoms, rons, rdcs, cats = newmets(outarray, msst, mstarts, mends)
    yrobj = getyears(msst.dates, mstarts, mends)
    fyears = unique(year.(msst.dates))
    meanmet = meanmets(evanoms, rons, rdcs, fyears, msst.anomfn)
    evmets = eventmets(evanoms, msst.anomfn)
    annmets = annualmets(evanoms, rons, rdcs, yrobj.eventyears, yrobj.startyear, yrobj.endyear, fyears, msst.anomfn)
    outmeans = meanmetsout(outmean, meanmet, mask)
    outtrend = trendout(outtrend, annmets)
    outannualmets = annualmetricsout(outarray, annmets, mask)

    return evmets, outtrend, outmeans, outannualmets
end

#=
FIX:
if ndims(sst) == 3
    x, y, z = size(sst)
else
    x = size(sst)
end

outarray = zeros(x, y, length(msst.date)) # matrix -> array
outarray = zeros(length(msst.date)) # vector 

- For the vector we have two choices: 
1. return a vector of vectors one vector for each of mean, coeff, pvalue and r_sq containing the means of the metrics (about 8 of them) so Vector{Vector{Vector, 8},4}


outmeansv1a = [zeros(size(metrics)) for _ in 1:4] # vector * count(mean, coeff, pvalue, r_squared) @btime [zeros(5) for _ in 1:4];139.946 ns (10 allocations: 480 bytes)

outmeansv1b = fill(zeros(size(metrics)), 4) # @btime fill(zeros(5), 4);56.918 ns (4 allocations: 192 bytes) NOTE: preferred.
- accessing: 
    a. getindex(outmeansv1b[1], 1) == the 1st metric of the first of coeff, mean, pval, rsquared.
    b. getindex.(outmeansv1b, 1) == the coeff, mean, pval, rsquared of the first metric
    c. getindex(outmeansv1b, 1) == the coeff for all metrics


2. return a matrix in which the rows correspond to the metrics and the columns correspond to the 4 variable outputs.

outmeansv2 = zeros(size(metrics), 4) # vector * length((mean, coeff, pvalue, r_squared)) in this case each row is a metric and each column is one of mean, coeff, pvalue and r_squared @btime zeros(5, 4); 40.205 ns (2 allocations: 240 bytes)
- accessing: 
    a. outmeansv2[:, 1] == all metrics in first column (e.g, coeff) so all coeffs.
    b. outmeansv2[1, :] == all of coeff, means, pval, rsq for first metric
    b. outmeansv2[1, 1] == first metric of first coeff (for example)

- Or we create an ntuple?

outmeans, outcoeff, outpvalue, outrsquared = ntuple(_ -> zeros(size(metrics)), 4) # @btime a, b,c,d = (zeros(5) for _ in 1:4);66.167 ns (4 allocations: 256 bytes), @btime a, b,c,d = ntuple(_ -> zeros(5), 4); 107.929 ns (8 allocations: 384 bytes)

outannuals1 = zeros(length(fullyears), size(metrics)) # vector input, matrix output column for each metric 20 years and 8 metrics = 20×8 Matrix # @btime zeros(20, 8); 164.210 ns (2 allocations: 1.38 KiB)

outannuals2 = fill(zeros(size(metrics)),length(fullyears)) # vector input, @btime fill(zeros(8),20); 82.647 ns (4 allocations: 352 bytes) NOTE: preferred.

- For the array input

# mean metrics and others (+3)
outmeans1 = fill(zeros(size(metrics)), x,y,z);  @btime fill(zeros(5), 4,5,3); 141.133 ns (4 allocations: 688 bytes), NOTE: preferred.
    - for accessing: x,y are dimensions of original, z is no of (coeff, mean, pvalue and rsquared)
    - for the inner metrics: getindex.(outmeans[1], 1) # outmeans[1] == coeff  and 1 is the meanint (for example)

outmeans2 = [zeros(size(metrics)) for _ in 1:x, _ in 1:y, _ in 1:z] # @btime [zeros(3,) for _ in 1:4, _ in 1:5, _ in 1:3]; 1.586 μs (122 allocations: 5.27 KiB)

# annual metrics
outannuals1 = fill(zeros(size(metrics)), x,y,z); z == no of fullyears PERF: preferred

outannuals = [zeros(size(metrics)) for _ in 1:x, _ in 1:y, _ in eachindex(fullyears)]
=#
