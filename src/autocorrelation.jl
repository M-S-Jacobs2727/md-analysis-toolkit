"""
    autocorrelation(data::Matrix{<:Number}, dtvalues::Vector{<:Real}=nothing)

# Positional Arguments

`data` (Matrix of Numbers): NxC matrix where each of the C columns is
correlated independently along N values equally spaced in time.

# Keyword Arguments

`dtvalues` (Vector of Reals): The time increments over which to perform the
correlations. By default, these values scale step-wise exponentially from 0
to 245760. More precisely, the values will be 0-15, then `[2^(i+3), 2^(i+4))`
with steps of `2^i` for `i=1:14` (maximum value of `2^18-2^14 = 245760`).
If given, the values are sorted and de-duplicated, then returned in that form.

# Return Values

`dtvalues` (Vector of Floats): See above.

`acf` (Matrix of Numbers): The autocorrelations of `data`, taken independently
along each column, with time increments given by `dtvalues`. Size 
`(length(dtvalues), C)`.

Recommended: read the data from a text file with

    using DelimitedFiles
    data = readdlm("output_data.txt", comments=true)
    dtvalues, acf = autocorrelation(data)
    
"""
function autocorrelation(data::Matrix{<:Number}, dtvalues::Vector{<:Real}=nothing)
    if dtvalues === nothing
        dtvalues = zeros(Int, 128)
        dtvalues[1:16] = 0:15
        for i in 1:14
            dtvalues[i*8+9:i*8+16] = range(2^(i+3), step=2^i, length=8)
        end
    else
        dtvalues |> sort! |> unique!
    end

    acf = zeros(eltype(data), size(data, 2), length(dtvalues))
    for (i, dt) in enumerate(dtvalues)
        for (r1, r2) in zip(eachrow(data[1:end-dt]), eachrow(data[dt+1:end]))
            acf[ : , i] += r1 .* conj(r2)
        end
    end

    return dtvalues, transpose(acf)
end