import LammpsFiles

"""
    radial_dist_func(filenames...; binwidth=0.05, bytype=false, maxdistance=nothing, ndim=3)

Compute the radial distribution function for a collection of LAMMPS dump
files (`filenames`). 

The width of the radial bins is set to `bin_width=0.05` by default.

If `bytype=false` (default) then the overall RDF is computed. If 
`true`, the RDF between each type is computed.

If given, the `maxdistance` dictates the position of the last bin. 
Otherwise, it is dictated by the box size (equivalent to `maxdistance=L`
where `L` is the box length).

The number of dimensions is set by `ndim` (default 3).
If `ndim=2`, then the 3rd dimension is not used to computed
periodic images.

## Returns
`binedges`: the maximum of each bin
`rdf`: the value of the RDF for each bin

"""
function radialDistFunc(filenames...; binwidth=0.05, bytype=false, maxdistance=nothing, ndim=3)
    ndim == 3 || ndim == 2 || throw(ArgumentError("argument ndim must be 2 or 3"))
    frame = LammpsFiles.read_dump(filenames[1])
    natoms = frame.natoms
    
    type_col = findfirst(x->x=="type", frame.properties)
    if bytype
        types = Int.(frame.atoms[type_col, : ])
        uniquetypes = unique(types)
        ntypes = length(uniquetypes)
        typecombos = [(i, j) for i in 1:ntypes for j in i:ntypes]
        typedict = Dict(type=>i+1 for (i, type) in enumerate(typecombos))
    else
        ntypes = 0
    end

    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["x", "y", "z"][1:ndim]]
    boxdims = frame.box[ : , 2] - frame.box[ : , 1]
    if maxdistance === nothing
        maxdistance = minimum(boxdims[1:ndim])
    end
    density = natoms / prod(boxdims[1:ndim])
    
    binedges = range(start=0, stop=maxdistance, step=binwidth)
    counts = zeros(Int, length(binedges), ntypes * (ntypes + 1) / 2 + 1)
    maxdistance = binedges[end]

    for filename in filenames
        frame = LammpsFiles.read_dump(filename)
        coords = frame.atoms[coord_cols, :]
        for i = 1 : natoms-1
            for j = i+1 : natoms
                differencevec = coords[ : , i] - coords[ : , j]
                for (d, l) in zip(differencevec, boxdims)
                    if (abs(d) > l / 2) differencevec -= l * sign(d) end
                end
                distance = hypot(differencevec)
                if distance >= maxdistance continue end
                bin = fld(distance, binwidth) + 1
                counts[bin, 1] += 1
                if bytype
                    index = typedict[minmax(types[i], types[j])]
                    counts[bin, index] += 1
                end
            end
        end
    end

    if ndim == 3
        rdf = counts ./ (density * 4pi/3 * (binedges.^3 - (binedges .- binwidth).^3) * length(filenames))
    else
        rdf = counts ./ (density * pi * (binedges.^2 - (binedges .- binwidth).^2) * length(filenames))
    end

    return (binedges, rdf)
end
