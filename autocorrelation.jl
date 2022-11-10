function compute_acf(data, dt)
    result = 0.0
    for (r1, r2) in zip(eachrow(data[1:end-dt]), eachrow(data[dt+1:end]))
        result += r1 .* r2
    end
    return result
end

"""
    autocorrelation(data)

    Only works for floating point data (not complex).

    Should be more than 2^18 = 262,144 rows

    Input:
        data: NxC matrix (N > 2^18)
    
    Output:
        acf: 128xC matrix
    
"""
function autocorrelation(data)
    acf = zeros(eltype(data), size(data, 2), 128)
    i = 1
    for dt in 0:15
        acf[i] = compute_acf(data, dt)
        i += 1
    end

    step = 2
    start = 16
    for _ in 1:14
        for j in 0:7
            dt = start + step * j
            acf[ : , i] = compute_acf(data, dt)
            i += 1
        end
        step *= 2
        start *= 2
    end

    return acf'
end