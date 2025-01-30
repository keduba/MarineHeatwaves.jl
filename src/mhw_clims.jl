

function tresh(input::Matrix, drange, thresh)
    inp = deepcopy(input)
    outq = [quantile!.(eachrow(inp[:, vcat(drange[d]...)]), thresh) for d in eachindex(drange)]
    outq[60] = mean((outq[59], outq[61]))
    return outq
end

function clim(input::Vector, drange)
    return [mean(findices(input, dd)) for dd in drange]
end

function tresh(input::Vector, drange, thresh)
    return [quantile(findices(input, dd), thresh) for dd in drange]
end

function clim(input::Matrix, drange)
    clima = [vec(mean(input[:, vcat(drange[i]...)], dims=2)) for i in eachindex(drange)]
    clima[60] = mean((clima[59], clima[61]))
    return clima
end

function clthr(input::Vector, drange, thresh)
    clima = clim(input, drange)
    climq = tresh(input, drange, thresh)
    _smoothdata!(clima)
    _smoothdata!(climq)
    return clima, climq
end

function clthr(input::Matrix, drange, thresh)
    clima = reduce(hcat, clim(input, drange))
    climq = reduce(hcat, tresh(input, drange, thresh))
    _smoothdata!.(eachrow(clima))
    _smoothdata!.(eachrow(climq))
    return clima, climq
end

"""
    _smoothdata!(ctarray::Union{Vector, SubArray}, pw) -> ctarray

    This function modifies `ctarray` in place. `ctarray` is the climatology/threshold array. It calculates and returns the moving mean. `pw` is the smoothing window. By default, this value is `31`.
    This is a test: lastindex(B) == 3N || error("The smoothing failed. Please check the `_smoothdata!` function.")
"""
function _smoothdata!(ctarray, pw=pctwidth)
    B = repeat(ctarray, 3)
    N = length(ctarray)
    ctarray[begin:end] = view(movmean(B, pw), N+1:2N)
    return ctarray
end


