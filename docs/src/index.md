```@meta
CurrentModule = MarineHeatwaves
``` 

# MarineHeatwaves

Documentation for [MarineHeatwaves](https://github.com/keduba/MarineHeatwaves.jl).

Marine heatwaves are extremes in the ocrean that have led to studies among scientists. Due to the importance, several packages have been created in the efforts to study and understand their occurrence.

MarineHeatwaves.jl is one such attempt to create a package in Julia that solves the problem, written in a high performance language and making use of native julia.

The functions and methods are also applicable to atmospheric heatwaves. The principal method used in the package is Hobday et al.(2015, 2018). Other newer methods are also going to be added in a standardised form.

## Installing

Well, first we have to make sure we have Julia on our computer.

1. [Download and install Julia](https://julialang.org/downloads). You can follow the instructions for your PC (macOS, Linux or Windows).


2. Launch Julia by typing `julia` in your terminal.

3. And then, install `MarineHeatwaves`:

```julia
]add MarineHeatwaves
```

It will take a few minutes to set up and compile the package and the other packages it needs.

4. When it's done, hit `backspace` to exit the package manager.

5. To use `MarineHeatwaves`:
```julia
using MarineHeatwaves
```
