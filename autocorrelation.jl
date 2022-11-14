"""
    autocorrelation(data)

Only works for real-valued data (not complex).

Recommended: read the data from a text file with

    using DelimitedFiles
    data = readdlm("output_data.txt", comments=true)
    acf = autocorrelation(data)

## Input

`data`: NxC matrix

`dt_values` (optional): vector of length T. If not provided, the length will
be 128, and the values will be 0-15, then `[2^(i+3), 2^(i+4))` with steps of
`2^i` for `i=1:14` (maximum value of `2^18-2^14 = 245760`).

## Returns

`dt_values`: vector of Ints (see above)

`acf`: TxC matrix

Note: 
    
"""
function autocorrelation(data, dt_values=nothing)
    if dt_values === nothing
        dt_values = zeros(Int, 128)
        dt_values[1:16] = 0:15
        for i in 1:14
            dt_values[i*8+9:i*8+16] = range(2^(i+3), step=2^i, length=8)
        end
    else
        dt_values |> sort! |> unique!
    end

    acf = zeros(Float64, size(data, 2), length(dt_values))
    for (i, dt) in enumerate(dt_values)
        for (r1, r2) in zip(eachrow(data[1:end-dt]), eachrow(data[dt+1:end]))
            acf[ : , i] += r1 .* r2
        end
    end

    return dt_values, transpose(acf)
end