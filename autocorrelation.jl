function compute_acf(data, dt)
    result = 0.0
    for (r1, r2) in zip(eachrow(data[1:end-dt]), eachrow(data[dt+1:end]))
        result += r1 .* r2
    end
    return result
end

"""
    autocorrelation(data)

    Only works for real-valued data (not complex).

    Recommended: Read the data from a text file with

        using DelimitedFiles
        data = readdlm("output_data.txt", comments=true)
        acf = autocorrelation(data)

    Input:
        data: NxC matrix
    
    Output:
        dt_values: 128-length vector
        acf: 128xC matrix
    
"""
function autocorrelation(data)
    acf = zeros(eltype(data), size(data, 2), 128)
    dt_values = zeros(Int, 128)
    i = 1
    for dt in 0:15
        dt_values[i] = dt
        acf[i] = compute_acf(data, dt_values[i])
        i += 1
    end

    step = 2
    start = 16
    for _ in 1:14
        for j in 0:7
            dt_values[i] = start + step * j
            acf[ : , i] = compute_acf(data, dt_values[i])
            i += 1
        end
        step *= 2
        start *= 2
        (step*2 > size(data, 1)) && break
    end

    return dt_values, acf'
end