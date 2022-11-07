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
    frame = LammpsFiles.read_dump(filenames[1])
    id_col = findfirst(x->x=="id", frame.properties)
    type_col = findfirst(x->x=="type", frame.properties)
    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["x", "y", "z"][1:ndim]]
    box_dims = frame.box[ : , 2] - frame.box[ : , 1]
    if max_distance === nothing
        max_distance = min(box_dims[1:ndim])
    end
    coords = frame.atoms[coord_cols, :]
    density = frame.natoms / prod(box_dims[1:ndim])

    # Do it for one frame first, by_type=false
    natoms = frame.natoms
    bin_centers = range(start=0.5*bin_width, stop=max_distance, step=bin_width)
    num_bins = length(bin_centers)
    counts = zeros(num_bins)
    for i = 1 : natoms-1
        for j = i+1 : natoms
            diff_vec = coords[ : , i] - coords[ : , j]
            for (d, l) in zip(diff_vec, box_dims)
                (abs(d) > l / 2) && (diff_vec -= l * sign(d))
            end
            dist = hypot(diff_vec)
            (dist <= max_distance) && (counts[fld(dist, bin_width) + 1] += 1)
        end
    end

    if ndim == 3
        rdf = counts ./ (density * 4pi * bin_width * bin_centers.^2)
    elseif ndim == 2
        rdf = counts ./ (density * 2pi * bin_width * bin_centers)
    end

    return (bin_centers, rdf)

end