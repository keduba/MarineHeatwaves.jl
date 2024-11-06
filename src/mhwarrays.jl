using Dates
using NaNStatistics

struct MCTemp{T <: AbstractFloat, Ti <: Integer}
    temp::VecOrMat{T}
    dates::StepRange{Date, Day}
    lyday::Vector{Ti}
    clima::VecOrMat
    thresh::VecOrMat
    pctile # = 0.1
    argfn  # = nanargmin
    anomfn # = nanminimum
    mask
end

struct MHTemp{T <: AbstractFloat, Ti <: Integer}
    temp::VecOrMat{T}
    dates::StepRange{Date, Day}
    lyday::Vector{Ti}
    clima::VecOrMat
    thresh::VecOrMat
    pctile::Float32 # = 0.9
    argfn::Function # = nanargmax
    anomfn::Function # = nanmaximum
    mask
end

struct MarineHW{T, N}
    temp::Array{T, N}
    category::Array{T, N}
    exceed::Array{Union{Missing, Bool}, N}
    annuals::NamedTuple
    means::NamedTuple
    trends::NamedTuple
    pvalues::NamedTuple
    pmetrics::NamedTuple
end


const winwidth = 5
const pctwidth = 31
const mindur = 5
const maxgap = 2

Base.length(fl::Base.Iterators.Flatten{<:AbstractVector{<:UnitRange}}) = sum(length, fl.it)

Base.IteratorSize(::Base.Iterators.Flatten{<:AbstractVector{<:UnitRange}}) = Base.HasLength()

"""
    seamask(sst) -> BitVecOrMat

    Returns only sea areas. Optional dimension. For arrays, expected dimension is the 3rd dimension.
"""
seamask(xy, dims=3) = dropdims(count(!isnan, xy, dims=dims) .> 0, dims=dims)


"""
    _subtemp(sst::Array{T, N}, sstindex) -> subsst, mask indices

    This function is used to subset the sst to the indices for the mhwdate and/or climatology date. It's not exported.
"""
function _subtemp(sst::Vector, mhwix) 
   xz = seamask(sst, 1)[1]
   if xz
       msst = view(sst, mhwix)
       return msst, xz
   else
       @warn "The temperature data is all NaNs. Please check again."
   end
end

function _subtemp(sst::Array, mhwix)
    xz = seamask(sst)
    xyz = CartesianIndices(xz)[xz]
    msst = view(sst, xyz, mhwix)
    return msst, xz, xyz
end

function subtemp(sst::Vector, sstdate, evdate)
    etix, elyd = timeindices(sstdate, evdate)
    evsst, mask = _subtemp(sst, etix)
    return evsst, mask, elyd
end

function findices(av, dr)
    (av[i] for i in flatten(dr))
end

function tresh(input::Matrix, drange, thresh)
    inp = deepcopy(input)
    outq = [quantile!.(eachrow(inp[:,vcat(drange[d]...)]), thresh) for d in eachindex(drange)]
    outq[60] = mean((outq[59], outq[61]))
    return outq
end

function clim(input::Matrix, drange)
    clima = [mean(input[:,vcat(drange[i]...)], dims = 2) for i in eachindex(drange)]
    clima[60] = mean((clima[59], clima[61]))
    return clima
end

clim(input::Vector, drange) = [mean(findices(input, dd)) for dd in drange]
tresh(input::Vector, drange, thresh) = [quantile(findices(input, dd), thresh) for dd in drange]

function clthr(input::Vector, drange, thresh)
    clima = clim(input, drange)
    climq = tresh(input, drange, thresh)
    _smoothdata!(clima)
    _smoothdata!(climq)
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
    ctarray .= view(movmean(B, pw), N+1:2N)
    return ctarray
end

function mhctemp(sst::Vector{T}, sstdate::StepRange{D,Dt}, mdate::StepRange{D,Dt}, cdate::StepRange{D, Dt}; threshold=0.1) where {T, D, Dt}
    mhsst, mask, mlyd = subtemp(sst, sstdate, mdate)
    clsst, mask, clyd = subtemp(sst, sstdate, cdate)
    dvec = daterange(clyd, winwidth)
    clima, climq = clthr(clsst, dvec, threshold)
    return mhsst, mdate, mlyd, mask, clima, climq
end






