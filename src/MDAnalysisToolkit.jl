module MDAnalysisToolkit
export square_Rg, autocorrelation, mean_square_displacement, radial_dist_func
include("mean_square_displacement.jl")
include("radial_dist_func.jl")
include("square_Rg.jl")
include("autocorrelation.jl")
end
