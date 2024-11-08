module MarineHeatwaves

# Write your package code here.
using Dates
using Statistics
using GLM
using NaNStatistics
using Base.Iterators: flatten
using DataFrames


include("timefunctions.jl")
include("mhwarrays.jl")
export leapyearday
export edetect

end
