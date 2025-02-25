"""
this is a helper function to calculate where the temperature exceeds the threshold. lyd is the leapyearday to create a longer vector so that threshdata is the same length as tempdata (temperature)

"""
_exceed(exfn, tempdata::Vector, threshdata::Vector, lyd) = exfn.(tempdata, threshdata[lyd])

_exceed(exfn, tempdata::Matrix, threshdata::Matrix, lyd) = exfn.(tempdata, threshdata[lyd, :])

"""
`exceed` is used to apply the exceedance function to the mhw or mcs types (`MCTemp` or `MHTemp`). It returns an array (matrix or vector). Should be a bitmatrix or bitvector or array of booleans.

exceed(x::MCTemp) -> exceedance::VecOrMat
Label the starts and ends using the exceedance and their differences.
"""
function exceed(x::Union{MCTemp,MHTemp})
    return _exceed(x.excfn, x.temp, x.thresh, x.lyday)
end
# TODO: review how the labeling function works from start to end and refactor.

# The endlabel and startlabel functions are vector based. They take a column of vectors of the exceedance and return the position.
# This may also be accomplished using the pointer in SparseColumns.

# TODO: how does the pointe in SparseCSC work? Could use it to replace this function.

function _endlabel(exdiff)# typeof(exceedance) == Vector
    menders = Int[i for (i, m) in enumerate(exdiff) if isequal(-1, m)]
    return menders
end

function _startlabel(exdiff)
    mstarts = Int[i for (i, m) in enumerate(exdiff) if isone(m)]
    return mstarts
end

# NOTE: HOW does the labeling work?
#
# First we obtain the points where the temperature is higher/less than the threshold using the `exceed` function. The result is a boolean array or bitarray (VecOrMat).
# Next, we get the pointer/index (SparseCSC or `startlabel/endlabel` function.) where the values meet the condition (> | <)
# Then, we take the difference within the column (so it is row difference). That means that the final row length should be n-1.
#

#=function _mylabeling(exceedance::BitVector)
    mexd = diff(exceedance)
    ss = ifelse(isone(first(exceedance)), count(==(1), mexd) + 1, count(==(1), mexd))
    es = ifelse(isone(last(exceedance)), count(==(-1), mexd) + 1, count(==(-1), mexd))
    ss == es || throw("Tengo frio! Something's off with the starts and ends in the `_mylabeling`")
    mstarts = startlabel(mexd)
    menders = endlabel(mexd)
    mstts, mends, hna = _indices(mstarts, menders, mindur, maxgap)
    mstartxs = startindices(mstts, hna)
    mendsxs = endindices(mends, hna)
    return mstartxs, mendsxs
end

function _mylabeling(exceedance::BitMatrix) 
    mexd = diff(exceedance, dims=2)
    f(exceedance, mexd) = ifelse(isone(first(exceedance)), count(==(1), mexd) + 1, count(==(1), mexd))
    f2(exceedance, mexd) = ifelse(isone(last(exceedance)), count(==(-1), mexd) + 1, count(==(-1), mexd))
    ss = map((x, y) -> f(x, y), eachrow(exceedance), eachrow(mexd))
    es = map((x, y) -> f2(x, y), eachrow(exceedance), eachrow(mexd))
    ss == es || throw("Tengo frio! Something's off with the starts and ends in the `_mylabeling`")

    mstarts, menders = ntuple(_ -> [ones(Int, e) for e in es], 2)
    for (m, mxd) in enumerate(eachrow(mexd))
        mstarts[m] = startlabel(mxd)
        menders[m] = endlabel(mxd)
    end

    mstartxs, mendsxs = ntuple(_ -> typeof(mstarts)(undef, size(exceedance, 1)), 2)
    for s in eachindex(mstarts, menders)
        mstts, mends, hna = _indices(mstarts[s], menders[s], mindur, maxgap)
        mstartxs[s] = startindices(mstts, hna)
        mendsxs[s] = endindices(mends, hna)
    end
    mstartxs, mendsxs
end=#

function __mylabeling(exceedance, exdiff)
    # first function checking that the length of starts and ends of events are similar
    S = ifelse(isone(first(exceedance)), count(==(1), exdiff) + 1, count(==(1), exdiff))
    E = ifelse(isone(last(exceedance)), count(==(-1), exdiff) + 1, count(==(-1), exdiff))
    S == E || error("Something's not right with the starts and ends in the exceedance vector.")
    # second function getting the starts and ends based on the length of the events.
    mstarts, menders = ntuple(_ -> ones(Int, S), 2)
    menders[end] = lastindex(exceedance)
    for (i, (k,)) in enumerate(Iterators.filter(p -> isone(p.second), pairs(exdiff)))
        isone(first(exceedance)) ? mstarts[i+1] = k : mstarts[i] = k
    end

    for (i, (k,)) in enumerate(Iterators.filter(p -> isequal(-1, p.second), pairs(exdiff)))
        menders[i] = k
    end
    return mstarts, menders
end

#    isone(first(exceedance)) 
#    ? mstarts[2:end] = startlabel(exdiff)
#    : mstarts[begin:end] = startlabel(exdiff)
# 
# WARN: No need for the following lines. Instead obtain the starts and ends and then check the length. if not equal, throw error. That's not all though.
#
# ss = isone(first(exceedance))
#  ? count(==(1), exdiff) + 1
#  : count(==(1), exdiff) 
#  
# es = isone(last(exceedance))
#  ? count(==(-1), exdiff) + 1
#  : count(==(-1), exdiff) 
# assert ss == es | "check the length of start and end"
#
function startlabel(exdiff, exceedance)
    mstarts = _startlabel(exdiff)
    if isone(first(exceedance))
        ss = count(==(1), exdiff) + 1
        mstarts = ones(Int, ss)
        mstarts[2:end] = _startlabel(exdiff)
    end
    return mstarts
end

function endlabel(exdiff, exceedance)
    menders = _endlabel(exdiff)
    if isone(last(exceedance))
        es = count(==(-1), exdiff) + 1
        menders = fill(lastindex(exceedance), es)
        menders[1:end-1] = _endlabel(exdiff)
    end
    return menders
end


#
# function _mylabeling(exceedance::BitVector)
#     exdiff = diff(exceedance)
#     mstarts, menders = __mylabeling(exceedance, exdiff)
#     mstts, mends, hna = _indices(mstarts, menders, mindur, maxgap)
#     mstartxs = startindices(mstts, hna)
#     mendsxs = endindices(mends, hna)
#     return mstartxs, mendsxs
# end

function _mylabeling(exceedance::BitVector, min_dur, max_gap)
    exdiff = diff(exceedance)
    mstarts = startlabel(exdiff, exceedance)
    menders = endlabel(exdiff, exceedance)
    mstts, mends, hna = _indices(mstarts, menders, min_dur, max_gap)
    mstartxs = startindices(mstts, hna)
    mendsxs = endindices(mends, hna)
    return mstartxs, mendsxs
end

# function _mylabeling(exceedance::BitMatrix, min_dur, max_gap)
#     exdiff = diff(exceedance, dims=1)
#     mstarts, menders = ([Int[] for _ in axes(exceedance, 2)] for _ in 1:2)
#     for (m, (colc, cold)) in enumerate(zip(eachrow(exceedance), eachrow(exdiff)))
#         mstarts[m], menders[m] = __mylabeling(colc, cold)
#     end
#     mstartxs, mendsxs = ntuple(_ -> typeof(mstarts)(undef, size(exceedance, 2)), 2)
#     for s in eachindex(mstarts, menders)
#         mstts, mends, hna = _indices(mstarts[s], menders[s], mindur, maxgap)
#         mstartxs[s] = startindices(mstts, hna)
#         mendsxs[s] = endindices(mends, hna)
#     end
#     return mstartxs, mendsxs
# end

function mylabeling(exceedance::Matrix)
    exdiff = diff(exceedance, dims=1)
    mstartsxs, mendersxs = ([Int[] for _ in axes(exceedance, 2)] for _ in 1:2)
    for x in axes(exdiff, 2)
        mstarts = startlabel(exdiff[x])
        mendser = endlabel(exdiff[x])

        mstts, mends, hna = _indices(mstarts, menders, mindur, maxgap)
        # can make the mindur, maxgap, thresholdand other exposables either a kwdef struct or a named tuple to make carrying them and modifying them 
        mstartxs[x] = startindices(mstts, hna)
        mendsxs[x] = endindices(mends, hna)
    end

end
"""
    The indices in the vector/matrix where the events begin and end.

_indices(sts, ends, minduration, maximumgap) -> (stixs, enixs)::Tuple{2, Vector{Integer}}
"""
function _indices(mstarts, menders, mindur, maxgap)
    durations = menders - mstarts
    # first the indices where the number of days > minimum duration
    mstts, mends = (ix[durations.≥mindur] for ix in (mstarts, menders))
    # then only those gap between starts and ends aremore than the specified maximum gap
    hna = (mstts[2:end] - mends[1:end-1]) .> maxgap
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
    mendsxs[1:end-1] = mends[1:end-1][hna]
    return mendsxs
end
