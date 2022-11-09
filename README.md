# MD Analysis Toolkit

This is a library of simple functions written in Julia to analyze MD simulations
in post via LAMMPS dump files (see my [LammpsFiles](https://code.ornl.gov/8mj/LammpsFiles)
project, a dependency).

These functions should be called with the form

```julia
using MDAnalysis

timesteps = collect(10:10:10000)
dt = 0.002
times = dt * timesteps
filenames = ["dump.$(lpad(i, 9, '0')).txt" for i in timesteps]
msd = mean_square_displacement(filenames)
bins, rdf = radial_dist_func(filenames)
rg2 = square_Rg(filesnames)

using Plots
pyplot()
plot(times, msd)
plot(bins, rdf)
println(rg2)
```
