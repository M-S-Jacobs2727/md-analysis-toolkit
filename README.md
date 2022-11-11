# MD Analysis Toolkit

This is a library of simple functions written in Julia to analyze MD simulations in post via LAMMPS dump files (see my [LammpsFiles](https://code.ornl.gov/8mj/LammpsFiles) project, a dependency).

These functions should be called with the following form:

```julia
using MDAnalysis

timesteps = 10:10:10000
dt = 0.002
times = dt * timesteps
filenames = ["dump.$(lpad(i, 9, '0')).txt" for i in timesteps]
msd = mean_square_displacement(filenames)
bins, rdf = radial_dist_func(filenames, by_type=true)
rg2 = square_Rg(filenames)

using DelimitedFiles
input_data = readdlm("time_averaged_data.txt", comments=true)
dt_values, corr_data = autocorrelation(input_data)
dt_values *= dt * 10

using Plots
pyplot()
plot(times[1:100], msd)
plot(bins, rdf)
plot(times, rg2)
plot(dt_values, corr_data, xscale=:log10)
```
