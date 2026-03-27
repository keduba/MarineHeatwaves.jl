# Obtaining the events metrics

## Key metrics

The key metrics for MHWs or MCSs as defined by Hobday et al. are contained in a constant `const Metrics` which is a dictionary of key-value pairs.

These metrics are calculated from the `Events` which is passed in as an argument to `meanmetrics` and `annualmetrics`.

Obtaining the `trend` requires the output from `annualmetrics`.

## Calculating the metrics

```@docs
MarineHeatwaves.meanmetrics
MarineHeatwaves.annualmetrics
MarineHeatwaves.linreg
```

## Returning and viewing the metrics

```@docs
MarineHeatwaves.mhmetrics
MarineHeatwaves.mymetric
```
