using Statistics

import LammpsFiles

"""
    mean_square_displacement(filenames...; max_num_frames=100, types=[], bymol=false, masses=nothing)

## Implemented
Computes the MSD for a collection of LAMMPS dump files. The range of
increments is from every 1 snapshot to every `max_num_frames`
snapshots (100 by default).

If `types` is given, then each type is independently evaluated, and an
overall MSD is given as well. Types are only checked for the first dump
file, so if atoms change type, that will invalidate the calculation.
By default, the MSD for each particle is averaged.

## Not implemented
If `by_mol` is set to true, then the molecule property must be available in
the dump file. The center of mass of each molecule will be analyzed
instead of each atom. This assumes that the mass of each atom type is 1.
If not, then the masses keyword should be set to a vector where each
element `mass[i]` is the mass of atoms of type `i`.

Note: cannot analyze by type and by molecule at the same time.
"""
function meanSquareDisplacement(filenames...; max_num_frames=100, types=nothing, bymol=false, masses=nothing)
    # Read basic info
    frame = LammpsFiles.read_dump(filenames[1])
    id_col = findfirst(x->x=="id", frame.properties)
    type_col = findfirst(x->x=="type", frame.properties)
    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["xu", "yu", "zu"]]
    max_num_frames = min(max_num_frames, length(filenames) - 1)

    # Set up sorting
    coords = frame.atoms[coord_cols, : ]
    ids = frame.atoms[id_col, : ]
    sortindices = sortperm(ids)
    sorted = true
    if !issorted(ids)
        sorted = false
        coords = coords[:, sortindices]
    end
    
    # Set up types
    bytype = types !== nothing
    if bytype
        selections = [
            findall(x->round(x)==t,
                    frame.atoms[type_col, : ][sortindices])
            for t in types
        ]
    end

    # Read in coordinates
    trajectories = zeros(Float64, 3, frame.natoms, length(filenames))
    trajectories[ : , : , 1] = coords

    for (i, f) in enumerate(filenames[2:end])
        frame = LammpsFiles.read_dump(f)
        coords = frame.atoms[coord_cols, : ]
        if sorted
            trajectories[ : , : , i+1] = coords
        else
            trajectories[ : , : , i+1] = coords[ : , sortperm(frame.atoms[id_col, : ])]
        end
    end

    # Compute msd for all and for types
    msd = zeros(
        Float64,
        bytype ? length(types) + 1 : 1,
        max_num_frames
    )
    for dt=1:max_num_frames
        diff = trajectories[ : , : , dt+1:end] - trajectories[ : , : , 1:end-dt]
        msd[1, dt] = mean(sum(diff .* diff, dims=1))

        if !bytype continue end

        for (i, sel) in selections
            msd[i+1, dt] = mean(sum(diff[ : , sel, : ] .* diff[ : , sel, : ], dims=1))
        end
    end
    return msd
end

