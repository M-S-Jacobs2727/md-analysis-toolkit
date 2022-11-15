module MDAnalysisToolkit
export squareRg, autocorrelation, meanSquareDisplacement, radialDistFunc
include("meanSquareDisplacement.jl")
include("radialDistFunc.jl")
include("squareRg.jl")
include("autocorrelation.jl")
end
