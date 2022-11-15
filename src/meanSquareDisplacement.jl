using Statistics

import LammpsFiles

"""
    meanSquareDisplacement(filenames::Vector{<:AbstractString}; max_num_frames::Integer=100, types::Vector{<:Integer}=nothing, bymol::Bool=false, masses::Vector{<:Real}=nothing)

Computes the mean-squared displacement of particles or molecules from a
collection of LAMMPS dump files.

## Positional Arguments

`filenames` (Vector of Strings): Paths leading to LAMMPS dump files.

## Keyword Arguments

`max_num_frames` (Integer): The time increments `dt` over which the MSD
is computed ranges from 1 to this value (default 100).

`types` (Vector of Integers): The atom types to independently analyze.
By default, only the global average is computed.

### Not yet implemented

`bymol` (Bool): Default is `false`. If `true`, the MSD of the center of 
mass of each molecule is computed instead. This assumes that the mass of
each atom type is 1. If not, then the masses keyword should be set.

`masses` (Vector of Reals): Only relevent when `bymol` is set to `true`
to compute the center of mass per molecule correctly. Should contain
one value per atom type such that `masses[i]` is the mass of atom type `i`.

## Return Values

`msd` (Matrix of Reals): Values of the mean square displacement. If `types`
is given, the size is `(max_num_frames, length(types)+1)`, where the first
column is the global average and each other column corresponds to the
respective atom type in `types`. Otherwise, the size is `(max_num_frames, 1)`.

If the time between dump frames is `dt`, then row `i` corresponds to a time
increment of `i*dt`.

Note: meanSquareDisplacement cannot analyze by type and by molecule at the
same time.
"""
function meanSquareDisplacement(filenames::Vector{<:AbstractString};
    max_num_frames::Integer=100, types::Vector{<:Integer}=nothing,
    bymol::Bool=false, masses::Vector{<:Real}=nothing
)
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
        max_num_frames,
        bytype ? length(types) + 1 : 1
    )
    for dt=1:max_num_frames
        diff = trajectories[ : , : , dt+1:end] - trajectories[ : , : , 1:end-dt]
        msd[dt, 1] = mean(sum(diff .* diff, dims=1))

        if !bytype continue end

        for (i, sel) in selections
            msd[dt, i+1] = mean(sum(diff[ : , sel, : ] .* diff[ : , sel, : ], dims=1))
        end
    end
    return msd
end

