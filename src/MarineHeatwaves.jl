module MarineHeatwaves

# Write your package code here.
using Dates
using Statistics
using NaNStatistics


include("timefunctions.jl")
include("mhwarrays.jl")
export leapyearday

end
