module MarineHeatwaves

# Write your package code here.
using Base.Iterators
import Statistics: mean, quantile, std
import Distributions: cdf, FDist, TDist
using Dates
# using NCDatasets
using SparseArrays


include("fullmhw.jl")
include("linreg.jl")
# include("timefunctions.jl")
# include("mhwarrays.jl")
# include("mhw_clims.jl")
# include("mhw_events.jl")
# include("mhw_labels.jl")
# include("mhw_metrics.jl")

export leapyearday
# export edetect
# export evtable
# export evclim
# export evthresh
export MarineCS
export MarineHW

end
