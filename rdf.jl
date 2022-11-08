using Statistics

import LammpsFiles

"""
    rdf(filenames...; bin_width=0.05, by_type=false, max_distance=nothing, ndim=3)

    Compute the radial distribution function for a collection of LAMMPS dump
    files (`filenames`). 

    The width of the radial bins is set to `bin_width=0.05` by default.

    If `by_type` is false (default) then the overall RDF is computed. If 
    true, the RDF between each type is computed.

    If given, the `max_distance` dictates the position of the last bin. 
    Otherwise, it is dictated by the box size (equivalent to `max_distance=L`
    where `L` is the box length).

    The number of dimensions is set by `ndim` (default 3).
    If `ndim == 2`, then the 3rd dimension is not used to computed
    periodic images.

"""
function rdf(filenames...; bin_width=0.05, by_type=false, max_distance=nothing, ndim=3)
    ndim == 2 || ndim == 3 || throw(ArgumentError("argument ndim must be 2 or 3"))
    frame = LammpsFiles.read_dump(filenames[1])
    
    type_col = findfirst(x->x=="type", frame.properties)
    if by_type
        types = frame.atoms[type_col, : ]
        unique_types = unique(types)
        ntypes = length(unique_types)
        type_combos = [(i, j) for i in 1:ntypes for j in i:ntypes]
        type_dict = Dict(type_combos[i]=>i for i in eachindex(type_combos))
    else
        ntypes = 0
    end


    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["x", "y", "z"][1:ndim]]
    box_dims = frame.box[ : , 2] - frame.box[ : , 1]
    if max_distance === nothing
        max_distance = min(box_dims[1:ndim])
    end
    coords = frame.atoms[coord_cols, :]
    density = frame.natoms / prod(box_dims[1:ndim])

    # Do it for one frame first, by_type=false
    natoms = frame.natoms
    bin_edges = range(start=bin_width, stop=max_distance, step=bin_width)
    num_bins = length(bin_edges)
    counts = zeros(num_bins, ntypes * (ntypes + 1) / 2 + 1)
    for i = 1 : natoms-1
        for j = i+1 : natoms
            diff_vec = coords[ : , i] - coords[ : , j]
            for (d, l) in zip(diff_vec, box_dims)
                (abs(d) > l / 2) && (diff_vec -= l * sign(d))
            end
            dist = hypot(diff_vec)
            if dist > max_distance
                continue
            end
            bin = fld(dist, bin_width) + 1
            counts[bin, 1] += 1
            if by_type
                ind = type_dict[minmax(types[i], types[j])]
                counts[bin, ind + 1] += 1
            end
        end
    end

    if ndim == 3
        rdf = counts ./ (density * 4pi * bin_width * bin_edges.^2)
    else
        rdf = counts ./ (density * 2pi * bin_width * bin_edges)
    end

    return (bin_edges, rdf)
end
