"""
tresh(indata::VecOrMat, daterange, threshold) -> climthresh::VecorMat

`tresh` calculates the climatological threshold of the `indata` to return the percentile at each point using the `threshold`. It uses `daterange` to select the specific days in the vector or matrix of the `indata` for which the desired percentile is returned.
"""

function tresh(input::Matrix, drange, thresh)
    inp = deepcopy(input)
    outthresh = [quantile!.(eachrow(inp[vcat(drange[d]...), :]), thresh) for d in eachindex(drange)]
    outthresh[60] = mean((outthresh[59], outthresh[61]))
    return outthresh
end

function tresh(input::Vector, drange, thresh)
    return [quantile(findices(input, dd), thresh) for dd in drange]
end

"""
clim(indata::VecOrMat, daterange) -> climamean::VecOrMat

`clim` calculates the climatological mean of the `indata` using `daterange` to select the specific days in the vector or matrix of the `indata` for which the desired mean is returned.
"""

function clim(input::Vector, drange)
    return [mean(findices(input, dd)) for dd in drange]
end

function clim(input::Matrix, drange)
    climamean = [vec(mean(input[vcat(drange[i]...), :], dims=2)) for i in eachindex(drange)]
    climamean[60] = mean((climamean[59], climamean[61]))
    return climamean
end

# TODO: add a way to target the smoothwindow
# ensure that the reduce(hcat,...) in the clthr function is not error prone
"""
clthr(indata::VecOrMat, daterange, threshold) -> climamean, climthresh

`clthr` combines the `clim` and `thresh` to return the climatological mean and thresholds after running `smoothdata!` to smooth the data for the given window `pctwidth` (default is 31).
"""

function clthr(input::Vector, drange, thresh, pctwidth)
    clima = clim(input, drange)
    climq = tresh(input, drange, thresh)
    _smoothdata!(clima, pctwidth)
    _smoothdata!(climq, pctwidth)
    return clima, climq
end

function clthr(input::Matrix, drange, thresh, pctwidth)
    clima = reduce(hcat, clim(input, drange))
    climq = reduce(hcat, tresh(input, drange, thresh))
    _smoothdata!.(eachcol(clima), pctwidth)
    _smoothdata!.(eachcol(climq), pctwidth)
    return clima, climq
end

"""
    _smoothdata!(ctarray::Union{Vector, SubArray}, pw) -> ctarray

    This function modifies `ctarray` in place. `ctarray` is the climatology/threshold array. It calculates and returns the moving mean. `pw` is the smoothing window. By default, this value is `31`.
"""
function _smoothdata!(ctarray, pw=pctwidth)
    B = repeat(ctarray, 3)
    N = length(ctarray)
    lastindex(B) == 3N || error("The smoothing failed. Please check the `_smoothdata!` function.")
    ctarray[begin:end] = view(movmean(B, pw), N+1:2N)
    return ctarray
end


