module MDAnalysisToolkit
export squareRg, autocorrelation, meanSquareDisplacement, radialDistFunc
include("mean_square_displacement.jl")
include("radial_dist_func.jl")
include("square_Rg.jl")
include("autocorrelation.jl")
end
