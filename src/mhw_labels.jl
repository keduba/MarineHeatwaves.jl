"""
this is a helper function to calculate where the temperature exceeds the threshold. lyd is the leapyearday to create a longer vector so that threshdata is the same length as tempdata (temperature)

"""
_exceed(exfn, tempdata::Vector, threshdata::Vector, lyd) = exfn.(tempdata, threshdata[lyd])

_exceed(exfn, tempdata::Matrix, threshdata::Matrix, lyd) = exfn.(tempdata, threshdata[:, lyd])

"""
`exceed` is used to apply the exceedance function to the mhw or mcs types (`MCTemp` or `MHTemp`). It returns an array (matrix or vector). Should be a bitmatrix or bitvector or array of booleans.

exceed(x::MCTemp) -> exceedance::VecOrMat
Label the starts and ends using the exceedance and their differences.
"""
function exceed(x::Union{MCTemp,MHTemp})
    return _exceed(x.excfn, x.temp, x.thresh, x.lyday)
end
# TODO: review how the labeling function works from start to end and refactor.
function endlabel(mexd)
    menders = Int[i for (i, m) in enumerate(mexd) if isequal(-1, m)]
    return menders
end

function startlabel(mexd)
    mstarts = Int[i for (i, m) in enumerate(mexd) if isone(m)]
    return mstarts
end
#=function _mylabeling(mexc::BitVector)
    mexd = diff(mexc)
    ss = ifelse(isone(first(mexc)), count(==(1), mexd) + 1, count(==(1), mexd))
    es = ifelse(isone(last(mexc)), count(==(-1), mexd) + 1, count(==(-1), mexd))
    ss == es || throw("Tengo frio! Something's off with the starts and ends in the `_mylabeling`")
    mstarts = startlabel(mexd)
    menders = endlabel(mexd)
    mstts, mends, hna = _indices(mstarts, menders, mindur, maxgap)
    mstartxs = startindices(mstts, hna)
    mendsxs = endindices(mends, hna)
    return mstartxs, mendsxs
end

function _mylabeling(mexc::BitMatrix) 
    mexd = diff(mexc, dims=2)
    f(mexc, mexd) = ifelse(isone(first(mexc)), count(==(1), mexd) + 1, count(==(1), mexd))
    f2(mexc, mexd) = ifelse(isone(last(mexc)), count(==(-1), mexd) + 1, count(==(-1), mexd))
    # ss = [f(rwc, rwd) for (rwc, rwd) in eachrow(mexd)]
    # es = [count(==(-1), row) for row in eachrow(mexd)]
    ss = map((x, y) -> f(x, y), eachrow(mexc), eachrow(mexd))
    es = map((x, y) -> f2(x, y), eachrow(mexc), eachrow(mexd))
    ss == es || throw("Tengo frio! Something's off with the starts and ends in the `_mylabeling`")

    mstarts, menders = ntuple(_ -> [ones(Int, e) for e in es], 2)
    for (m, mxd) in enumerate(eachrow(mexd))
        mstarts[m] = startlabel(mxd)
        menders[m] = endlabel(mxd)
    end

    mstartxs, mendsxs = ntuple(_ -> typeof(mstarts)(undef, size(mexc, 1)), 2)
    for s in eachindex(mstarts, menders)
        mstts, mends, hna = _indices(mstarts[s], menders[s], mindur, maxgap)
        mstartxs[s] = startindices(mstts, hna)
        mendsxs[s] = endindices(mends, hna)
    end
    mstartxs, mendsxs
end=#

function __mylabeling(mexc, mexcd)
    S = ifelse(isone(first(mexc)), count(==(1), mexcd) + 1, count(==(1), mexcd))
    E = ifelse(isone(last(mexc)), count(==(-1), mexcd) + 1, count(==(-1), mexcd))
    S == E || error("Something's not right with the starts and ends in the exceedance vector.")
    mstarts, menders = ntuple(_ -> ones(Int, S), 2)
    menders[end] = lastindex(mexc)
    for (i, (k,)) in enumerate(Iterators.filter(p -> isone(p.second), pairs(mexcd)))
        isone(first(mexc)) ? mstarts[i+1] = k : mstarts[i] = k
    end
    for (i, (k,)) in enumerate(Iterators.filter(p -> isequal(-1, p.second), pairs(mexcd)))
        menders[i] = k
    end
    return mstarts, menders
end

function _mylabeling(mexc::BitVector)
    mexcd = diff(mexc)
    mstarts, menders = __mylabeling(mexc, mexcd)
    mstts, mends, hna = _indices(mstarts, menders, mindur, maxgap)
    mstartxs = startindices(mstts, hna)
    mendsxs = endindices(mends, hna)
    return mstartxs, mendsxs
end

function _mylabeling(mexc::BitMatrix)
    mexcd = diff(mexc, dims=2)
    mstarts, menders = ([Int[] for _ in axes(mexc, 1)] for _ in 1:2)
    for (m, (colc, cold)) in enumerate(zip(eachrow(mexc), eachrow(mexcd)))
        mstarts[m], menders[m] = __mylabeling(colc, cold)
    end
    mstartxs, mendsxs = ntuple(_ -> typeof(mstarts)(undef, size(mexc, 1)), 2)
    for s in eachindex(mstarts, menders)
        mstts, mends, hna = _indices(mstarts[s], menders[s], mindur, maxgap)
        mstartxs[s] = startindices(mstts, hna)
        mendsxs[s] = endindices(mends, hna)
    end
    return mstartxs, mendsxs
end

"""
    The indices in the vector/matrix where the events begin and end.

_indices(sts, ends, minduration, maximumgap) -> (stixs, enixs)::Tuple{2, Vector{Integer}}
"""
function _indices(mstarts, menders, mindur, maxgap)
    oldurations = menders - mstarts
    mstts, mends = (ix[oldurations.≥mindur] for ix in (mstarts, menders))
    hna = mstts[2:end] - mends[1:end-1] .> maxgap
    return (mstts, mends, hna)
end

function startindices(mstts, hna)
    mstartsxs = ones(Int, sum(hna) + 1)
    mstartsxs[1] = first(mstts)
    mstartsxs[2:end] = mstts[2:end][hna] .+ 1
    return mstartsxs
end

function endindices(mends, hna)
    mendsxs = ones(Int, sum(hna) + 1)
    mendsxs[end] = last(mends)
    mendsxs[begin:end-1] = mends[begin:end-1][hna]
    return mendsxs
end


