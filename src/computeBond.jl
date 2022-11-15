using Statistics

using LammpsFiles

"""
    computeBond(datafile::AbstractString, bytype::Bool=false, dumpfiles::AbstractString...=nothing;
    atomstyle::AbstractString="full", return_bondlengths::Bool=false)

Compute the average, standard deviation, mean-squared average, fluctuations,
and distribution of bond lengths in a LAMMPS data file, optionally including
a series of dump files for averaging. By default, computes the global
values independent of bond types. If `bytype` is true, then the values
are computed independently for each bond type in the data file and
returned as a vector, and the optional field `bonds.bondlengths` is returned
as a vector of vectors of floats.

Returns a named tuple `bonds` with the following fields:
- `average` (Float): The standard mean of the bond lengths.
- `stddev` (Float): The standard deviation of the bond lengths
(normalized by the number of bonds).
- `meansq` (Float): The standard mean of the squared bond lengths.
- `fluctuation` (Float): The difference between `meansq` and the
square of `average`.
- `bondlengths` (Vector of Floats, optional): The length of each
bond in each frame read. To be used for distribution analysis.
"""
function computeBond(datafile::AbstractString, bytype::Bool=false, dumpfiles::AbstractString...=nothing;
    atomstyle::AbstractString="full", return_bondlengths::Bool=false)

    data = LammpsFiles.readData(datafile, atomstyle=atomstyle)
    if bytype
        ubondtypes = data.bond_types |> sort |> unique!
        return computeBond(datafile, ubondtypes, dumpfiles, atomstyle=atomstyle, return_bondlengths=return_bondlengths)
    end

    bondlengths = hypot.(
        data.coords[ : , data.idtoindex[data.bonds[1, : ]]] 
        - data.coords[ : , data.idtoindex[data.bonds[2, : ]]]
    )
    
    average = mean(bondlengths)
    stddev = stdm(bondlengths, average)
    meansq = mean(bondlengths.^2)
    fluctuation = meansq-average^2

    return return_bondlengths ? (
        average = average,
        stddev = stddev,
        meansq = meansq,
        fluctuation = fluctuation,
        bondlengths = bondlengths,
    ) : (
        average = average,
        stddev = stddev,
        meansq = meansq,
        fluctuation = fluctuation,
    )
end

"""
    computeBond(datafile::AbstractString, bytype::Vector{Integer}, dumpfiles::AbstractString...=nothing;
    atomstyle::AbstractString="full", return_bondlengths::Bool=false)

Compute the same values, but independently for each bond type given.

Returns the same named tuple `bonds`, but each value is instead a vector of
length `length(bytype)`.

"""
function computeBond(datafile::AbstractString, bytype::Vector{Integer}, dumpfiles::AbstractString...=nothing;
    atomstyle::AbstractString="full", return_bondlengths::Bool=false)

    data = LammpsFiles.readData(datafile, atomstyle=atomstyle)
    bondsbytype = [
        data.bonds[ : , findall(t->(t==bt), data.bond_types)] 
        for bt in bytype
    ]
    bondlengths = [hypot.(
        data.coords[ : , data.idtoindex[bonds[1, : ]]] 
        - data.coords[ : , data.idtoindex[bonds[2, : ]]]
    ) for bonds in bondsbytype]
    
    average = mean.(bondlengths)
    stddev = stdm.(bondlengths, average)
    meansq = mean.([b.^2 for b in bondlengths])
    fluctuation = meansq-average.^2
    
    return return_bondlengths ? (
        average = average,
        stddev = stddev,
        meansq = meansq,
        fluctuation = fluctuation,
        bondlengths = bondlengths,
    ) : (
        average = average,
        stddev = stddev,
        meansq = meansq,
        fluctuation = fluctuation,
    )
end
