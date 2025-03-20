module MarineHeatwaves

# Write your package code here.
using Dates
using Statistics
using GLM
using NaNStatistics
using Base.Iterators
using DataFrames


include("timefunctions.jl")
include("mhwarrays.jl")
include("mhw_clims.jl")
include("mhw_events.jl")
include("mhw_labels.jl")
include("mhw_metrics.jl")
include("linreg.jl")

export leapyearday
export edetect
export evtable
export evclim
export evthresh
export MarineCS
export MarineHW

end
