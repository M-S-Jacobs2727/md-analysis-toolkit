using Statistics

import LammpsFiles

"""
    msd(filenames...; max_num_frames=100, sorted=false, by_type=false, mol=false, masses=nothing)

    [Implemented]
    Computes the MSD for a collection of LAMMPS dump files. The range of
    increments is from every 1 snapshot to every 'max_num_frames' 
    snapshots (100 by default). 

    If sorted == true, then the atoms are assumed to be in the same order
    in each file. If false (default), they are sorted by the column 'id'.

    If by_type is false (default), the MSD for each particle is averaged. 
    If by_type is true, then each type is independently evaluated, and an 
    overall MSD is given as well. Types are only checked for the first dump
    file, so if atoms change type, that will invalidate the calculation.
    
    [Not implemented]
    If mol is set to true, then the molecule property must be available in
    the dump file. The center of mass of each molecule will be analyzed 
    instead of each atom. This assumes that the mass of each atom type is
    1. If not, then the masses keyword should be set to a vector where each 
    element mass[i] is the mass of atoms of type i.
"""
function msd(filenames...; max_num_frames=100, sorted=false, by_type=false, mol=false, masses=nothing)
    frame = LammpsFiles.read_dump(filenames[1])
    id_col = findfirst(x->x=="id", frame.properties)
    type_col = findfirst(x->x=="type", frame.properties)
    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["xu", "yu", "zu"]]
    max_num_frames = min(max_num_frames, length(filenames) - 1)

    coords = frame.atoms[coord_cols, : ]
    if !sorted
        sort_inds = sortperm(frame.atoms[id_col, : ])
        coords = coords[:, sort_inds]
    end

    all_coords = zeros(3, frame.natoms, length(filenames))
    all_coords[ : , : , 1] = coords

    if by_type
        types = unique(frame.atoms[type_col, : ])
        selections = [findall(x->round(x)==t, frame.atoms[type_col, : ]) for t in types]
    else
        types = []
    end

    for (i, f) in enumerate(filenames[2:end])
        frame = LammpsFiles.read_dump(f)
        coords = frame.atoms[coord_cols, : ]
        if !sorted
            sort_inds = sortperm(frame.atoms[id_col, : ])
            coords = coords[:, sort_inds]
        end
        all_coords[ : , : , i+1] = coords[ : , selection]
    end

    msd = zeros(length(types)+1, max_num_frames)
    for dt=1:max_num_frames
        diff = all_coords[ : , : , dt+1:end] - all_coords[ : , : , 1:end-dt]
        msd[1, dt] = mean(sum(diff .* diff, dims=1))
        if by_type
            for (i, sel) in selections
                msd[i+1, dt] = mean(sum(diff[ : , sel, : ] .* diff[ : , sel, : ], dims=1))
            end
        end
    end
    return msd
end

