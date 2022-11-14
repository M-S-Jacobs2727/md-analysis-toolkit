import LammpsFiles

"""
    radial_dist_func(filenames...; bin_width=0.05, by_type=false, max_distance=nothing, ndim=3)

Compute the radial distribution function for a collection of LAMMPS dump
files (`filenames`). 

The width of the radial bins is set to `bin_width=0.05` by default.

If `by_type=false` (default) then the overall RDF is computed. If 
`true`, the RDF between each type is computed.

If given, the `max_distance` dictates the position of the last bin. 
Otherwise, it is dictated by the box size (equivalent to `max_distance=L`
where `L` is the box length).

The number of dimensions is set by `ndim` (default 3).
If `ndim=2`, then the 3rd dimension is not used to computed
periodic images.

## Returns
`bin_edges`: the maximum of each bin
`rdf`: the value of the RDF for each bin

"""
function radial_dist_func(filenames...; bin_width=0.05, by_type=false, max_distance=nothing, ndim=3)
    ndim == 2 || ndim == 3 || throw(ArgumentError("argument ndim must be 2 or 3"))
    frame = LammpsFiles.read_dump(filenames[1])
    natoms = frame.natoms
    
    type_col = findfirst(x->x=="type", frame.properties)
    if by_type
        types = Int.(frame.atoms[type_col, : ])
        unique_types = unique(types)
        ntypes = length(unique_types)
        type_combos = [(i, j) for i in 1:ntypes for j in i:ntypes]
        type_dict = Dict(type=>i+1 for (i, type) in enumerate(type_combos))
    else
        ntypes = 0
    end

    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["x", "y", "z"][1:ndim]]
    box_dims = frame.box[ : , 2] - frame.box[ : , 1]
    if max_distance === nothing
        max_distance = minimum(box_dims[1:ndim])
    end
    density = natoms / prod(box_dims[1:ndim])
    
    bin_edges = range(start=0, stop=max_distance, step=bin_width)
    counts = zeros(Int, length(bin_edges) - 1, ntypes * (ntypes + 1) / 2 + 1)
    max_distance = bin_edges[end]

    for filename in filenames
        frame = LammpsFiles.read_dump(filename)
        coords = frame.atoms[coord_cols, :]
        for i = 1 : natoms-1
            for j = i+1 : natoms
                diff_vec = coords[ : , i] - coords[ : , j]
                for (d, l) in zip(diff_vec, box_dims)
                    (abs(d) > l / 2) && (diff_vec -= l * sign(d))
                end
                dist = hypot(diff_vec)
                (dist >= max_distance) && continue
                bin = fld(dist, bin_width) + 1
                counts[bin, 1] += 1
                if by_type
                    ind = type_dict[minmax(types[i], types[j])]
                    counts[bin, ind] += 1
                end
            end
        end
    end

    if ndim == 3
        rdf = counts ./ (density * 4pi/3 * (bin_edges.^3 - (bin_edges - bin_width).^3) * length(filenames))
    else
        rdf = counts ./ (density * pi * (bin_width^2 - 2bin_width * bin_edges) * length(filenames))
    end

    return (bin_edges, rdf)
end
