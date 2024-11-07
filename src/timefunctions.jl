


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
function daterange(lydd::Vector{T}, w::T) where T <: Integer
    nd = Vector{Vector{UnitRange{T}}}(undef, 366);
    for n in eachindex(nd)
        nd[n] = UnitRange{T}[max(1, x-w):min(lastindex(lydd), x+w) for (x, _) in  Iterators.filter(p -> isequal(n, p.second), pairs(lydd))]
    end
    return nd
end
