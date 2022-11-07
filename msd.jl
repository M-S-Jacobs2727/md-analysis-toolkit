using Statistics

import LammpsFiles

"""
    msd(filenames...; max_num_frames=100, sorted=false, types=(0,), mol=false, masses=[1,])

    [Implemented]
    Computes the MSD for a collection of LAMMPS dump files. The range of
    increments is from every 1 snapshot to every 'max_num_frames' 
    snapshots (100 by default). 
    
    If sorted == true, then the atoms are assumed to be in the same order
    in each file. If false (default), they are sorted by the column 'id'.
    
    [Not implemented]
    If types == (0,) (default), the MSD for each particle is averaged. If
    types is given and is a tuple of ints >= 1, then those types are 
    independently accounted for, and an overall MSD is given as well.
    
    If mol is set to true, then the molecule property must be available in
    the dump file. The center of mass of each molecule will be analyzed 
    instead of each atom. This assumes that the mass of each atom type is
    1. If not, then set the masses keyword to a vector where each element
    mass[i] is the mass of atoms of type i.
"""
function msd(filenames...; max_num_frames=100, sorted=false, types=(0,), mol=false, masses=[1,])
    frame = LammpsFiles.read_dump(filenames[1])
    all_coords = zeros(3, frame.natoms, length(filenames))
    for (i, f) in enumerate(filenames)
        frame = LammpsFiles.read_dump(f)
        inds = [findfirst(x->x==s, frame.properties) for s in ["xu", "yu", "zu"]]
        coords = frame.atoms[inds, : ]
        if !sorted
            sort_inds = sortperm(frame.atoms[findfirst(x->x=="id", frame.properties), : ])
            coords = coords[:, sort_inds]
        end
        all_coords[ : , : , i] = coords
    end

    msd = zeros(max_num_frames)
    for dt=1:max_num_frames
        diff = all_coords[ : , : , dt+1:end] - all_coords[ : , : , 1:end-dt]
        msd[dt] = mean(sum(diff .* diff, dims=1))
    end
    return msd
end

