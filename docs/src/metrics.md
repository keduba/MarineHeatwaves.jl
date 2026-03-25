# Obtaining the metrics.

## Key metrics

The key metrics for MHWs or MCSs as defined by Hobday et al. are contained in a constant `const metrics` which is a dictionary of key-value pairs.

These metrics are calculated from the `Events` which is passed in as an argument to `meanmetrics` and `annualmetrics`.

Obtaining the trend from `trend` requires the output from `annualmetrics`.

```@docs
MarineHeatwaves.meanmetrics
MarineHeatwaves.annualmetrics
MarineHeatwaves.mhmetrics
MarineHeatwaves.mymetric
MarineHeatwaves.linreg
```
