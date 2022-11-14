import LammpsFiles

"""
    radialDistFunc(filenames::Vector{<:AbstractString}; binwidth::Real=0.05, bytype::Bool=false, maxdistance::Real=nothing, ndim::Integer=3)

Compute the radial distribution function for a collection of LAMMPS dump
files.

## Positional Arguments

`filenames` (Vector of Strings): Paths leading to LAMMPS dump files.

## Keyword Arguments

`binwidth` (Real): Width of the radial bins. Default: 0.05.

`bytype` (Bool): If `false` (the default) then the overall RDF is computed. If 
`true`, the RDF between each type is computed.

`maxdistance` (Real): If given, indictates the position of the last bin. 
By default, it is dictated by the box size (equivalent to `maxdistance=L`
where `L` is the box length).

`ndim` (2 or 3): The number of spatial dimensions (default 3).
If `ndim=2`, then the 3rd (z) dimension is not used to compute
periodic images.

## Return Values

`binedges` (Vector of Reals): The edges of each bin (length N+1).

`columnnames` (Vector of Strings): The name of each column of `rdf`
(length C). The first is always 'all'. If `bytype` is `true`, then
each unordered pair of types '(m, n)' is given, corresponding to a
the respective column in `rdf`.

`rdf` (Matrix of Reals): The value of the RDF for each bin (size (N, C)).
Each column corresponds to the respective name in `columnnames`. Each row
`i` corresponds to the bin [`binedges[i]`, `binedges[i+1]`).

"""
function radialDistFunc(filenames::Vector{<:AbstractString}; binwidth::Real=0.05, bytype::Bool=false, maxdistance::Real=nothing, ndim::Integer=3)
    ndim == 3 || ndim == 2 || throw(ArgumentError("argument ndim must be 2 or 3"))
    frame = LammpsFiles.read_dump(filenames[1])
    natoms = frame.natoms
    
    type_col = findfirst(x->x=="type", frame.properties)
    typecombos = []
    ntypes = 0
    if bytype
        types = Int.(frame.atoms[type_col, : ])
        uniquetypes = unique(types)
        ntypes = length(uniquetypes)
        typecombos = [(i, j) for i in 1:ntypes for j in i:ntypes]
        typedict = Dict(type=>i+1 for (i, type) in enumerate(typecombos))
    end

    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["x", "y", "z"][1:ndim]]
    boxdims = frame.box[ : , 2] - frame.box[ : , 1]
    if maxdistance === nothing
        maxdistance = minimum(boxdims[1:ndim])
    end
    density = natoms / prod(boxdims[1:ndim])
    
    binedges = range(start=0, stop=maxdistance, step=binwidth)
    counts = zeros(Int, length(binedges)-1, length(typecombos) + 1)
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

    rdf = counts ./ (
        density
        * ((ndim==3) ? 4pi/3 : pi)
        * (binedges[2:end].^ndim - binedges[1:end-1].^ndim)
        * length(filenames)
    )

    columnnames::Vector{String} = ["all"; ["$a" for a in typecombos]]
    return (binedges, columnnames, rdf)
end
