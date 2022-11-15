module MDAnalysisToolkit
export squareRg, autocorrelation, meanSquareDisplacement, radialDistFunc, computeBond
include("meanSquareDisplacement.jl")
include("radialDistFunc.jl")
include("squareRg.jl")
include("autocorrelation.jl")
include("computeBond.jl")
end
