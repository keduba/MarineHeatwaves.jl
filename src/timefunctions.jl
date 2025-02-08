


"""
   Obtain the dayofyear of a date as though it fell in a leapyear.
   
   leapyearday(dt::Date) -> y::Integer 
   
   Examples 
   --------

   julia> leapyearday(Date(2013,03,31)) # non-leapyear
   91

   julia> dayofyear(Date(2016,03,31)) # leapyear
   91   
   
   julia> dayofyear(Date(2013,03,31))
   90   
"""
leapyearday(mts) = dayofyear(mts) > 59 && !isleapyear(mts) ? dayofyear(mts) + 1 : dayofyear(mts)

"""
Internal function for returning the indices of the start and end dates for subsetting the climatology and marine heatwave data arrays from the input temperature array, and the leapyearday vector of the second input argument.

timeindices(sdt::Union{Vector, StepRange{Date}}, edt::Union{Vector, StepRange{Date}}) -> Range, Vector{Int}
"""
function timeindices(sstdate, evtdate)
    si = findfirst(isequal(first(evtdate)), sstdate)
    ei = findfirst(isequal(last(evtdate)), sstdate)
    senix = range(si, ei)
    elyd = leapyearday.(sstdate[senix])
    return senix, elyd
end


"""
    Return a vector of the dates ranges based on the leapyearday.

daterange(lyd::Vector{Integer}, window::Integer=5) -> Vector{UnitRange{Integer}}
"""
function daterange(lydd::Vector{T}, w::T) where {T<:Integer}
    nd = Vector{Vector{UnitRange{T}}}(undef, 366)
    for n in eachindex(nd)
        nd[n] = UnitRange{T}[max(1, x - w):min(lastindex(lydd), x + w) for (x, _) in Iterators.filter(p -> isequal(n, p.second), pairs(lydd))]
    end
    return nd
end

# Test suite for the date variables: temperature, marineheatwave and climatology

# Test suite for the temperature array and the date.
function testarrays(sst, sstdate)
    errmsgdims = "Dimensions of temperature array and time do not match: "
    size(sst, ndims(sst)) == size(sstdate, ndims(sstdate)) || throw(DimensionMismatch("$(errmsgdims) temperature: $(size(sst, ndims(sst))), time: $(size(sstdate, ndims(sstdate)))."))
end

function testdates(sstdate, mhwdate, climdate)
    errmsglen = "The end should be greater than the start: "
    errmsgform = "The date should be either a string as in 'yyyy-mm-dd' -> '1990-01-02' or a date as in `Date(yyyy,mm,dd)` -> `Date(1990,3,12)`."
    errmsgrange = "TimeError: Period is outside temperature date period: "
    # errmsglen2 = "The length of period is greater than the length of `sstdate`. Check: "

    # Format/type tests
    (first(sstdate) isa String || first(sstdate) isa Date) || errmsgform

    (first(climdate) isa String || first(climdate) isa Date) || errmsgform
    (first(mhwdate) isa String || first(mhwdate) isa Date) || errmsgform

    # Assure dates are date objects
    ymd = "yyyy-mm-dd"
    sstdate = eltype(sstdate) == String ? range(Date(first(sstdate), ymd), Date(last(sstdate), ymd)) : sstdate
    mhwdate = eltype(mhwdate) == String ? range(Date(first(mhwdate), ymd), Date(last(mhwdate), ymd)) : mhwdate
    climdate = eltype(climdate) == String ? range(Date(first(climdate), ymd), Date(last(climdate), ymd)) : climdate

    # Make sure the order of the dates is correct
    first(sstdate) < last(sstdate) || errmsglen
    first(mhwdate) < last(mhwdate) || errmsglen
    first(climdate) < last(climdate) || errmsglen

    # Range tests
    first(mhwdate) ≥ first(sstdate) && last(mhwdate) ≤ last(sstdate) || error(errmsgrange, mhwdate)

    first(climdate) ≥ first(sstdate) && last(climdate) ≤ last(sstdate) || error(errmsgrange, climdate)

    climdate = range(max(climdate[1] - Day(5), first(sstdate)), min(climdate[end] + Day(5), last(sstdate)))
    return sstdate, mhwdate, climdate
end

