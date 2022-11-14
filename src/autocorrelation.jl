"""
    autocorrelation(data)

Only works for real-valued data (not complex).

Recommended: read the data from a text file with

    using DelimitedFiles
    data = readdlm("output_data.txt", comments=true)
    acf = autocorrelation(data)

## Input

`data`: NxC matrix

`dtvalues` (optional): vector of length T. If not provided, the length will
be 128, and the values will be 0-15, then `[2^(i+3), 2^(i+4))` with steps of
`2^i` for `i=1:14` (maximum value of `2^18-2^14 = 245760`).

## Returns

`dtvalues`: vector of Ints (see above)

`acf`: TxC matrix

Note: 
    
"""
function autocorrelation(data, dtvalues=nothing)
    if dtvalues === nothing
        dtvalues = zeros(Int, 128)
        dtvalues[1:16] = 0:15
        for i in 1:14
            dtvalues[i*8+9:i*8+16] = range(2^(i+3), step=2^i, length=8)
        end
    else
        dtvalues |> sort! |> unique!
    end

    acf = zeros(Float64, size(data, 2), length(dtvalues))
    for (i, dt) in enumerate(dtvalues)
        for (r1, r2) in zip(eachrow(data[1:end-dt]), eachrow(data[dt+1:end]))
            acf[ : , i] += r1 .* r2
        end
    end

    return dtvalues, transpose(acf)
end