using Statistics

using LammpsFiles

"""
    computeBond(datafile::AbstractString, dumpfiles::AbstractString...=nothing; bondtypes=nothing)

Compute the average, standard deviation, mean-squared average, fluctuations,
and distribution of bond lengths in a LAMMPS data file, optionally including
a series of dump files for averaging.

## Positional Arguments

`datafile` (String): Path leading to a LAMMPS data file.

`dumpfiles` (Strings, optional): Paths leading to LAMMPS dump files with all
the same atoms and bonds as the data file.

## Keyword Arguments

`bondtypes` (Vector of Integers): By default, all bonds are averaged together.
By specifying a list of bond types, the return values will only account for 
the types given.

## Return Values

`bonds` (Named Tuple): 
- `bondlengths` (Vector of Floats): The length of each bond in each frame read.
To be used for distribution analysis.
- `average` (Float): The standard mean of the bond lengths.
- `stddev` (Float): The standard deviation of the bond lengths (normalized
by the number of bonds).
- `meansq` (Float): The standard mean of the squared bond lengths.
- `fluctuation` (Float): The difference between `meansq` and the square of
`average`.

## TODO

- Change `bondtypes` to `bytype` so that it accounts for each bond type
separately, like `meanSquareDisplacement`.
- Implement dump file analysis.
"""
function computeBond(datafile::AbstractString, dumpfiles::AbstractString...=nothing; atomstyle="full", bondtypes=nothing)
    data = LammpsFiles.readData(datafile, atomstyle=atomstyle)
    bonds = data.bonds[ : , findall(t->t in bondtypes, data.bond_types)]
    bondlengths = hypot.(
        data.coords[ : , data.idtoindex[bonds[1, : ]]] 
        - data.coords[ : , data.idtoindex[bonds[2, : ]]]
    )
    average = mean(bondlengths)
    stddev = stdm(bondlengths, average)
    meansq = mean(bondlengths.^2)
    fluctuation = meansq-average^2
    return (
        bondlengths = bondlengths,
        average = average,
        stddev = stddev,
        meansq = meansq,
        fluctuation = fluctuation,
    )
end