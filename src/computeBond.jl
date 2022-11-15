using Statistics

using LammpsFiles

"""
    computeBond(datafile::AbstractString, dumpfiles::AbstractString...=nothing; bondtypes=nothing)

Compute the mean, standard deviation, variance, and distribution of bond
lengths in a LAMMPS data file, optionally including a series of dump files
for averaging.
"""
function computeBond(datafile::AbstractString, dumpfiles::AbstractString...=nothing; atomstyle="full", bondtypes=nothing)
    data = LammpsFiles.readData(datafile, atomstyle=atomstyle)
    bonds = data.bonds
    bond_types = data.bond_types
    coords = data.coords
    atom_ids = data.atom_ids
    idtoindex = data.idtoindex
    bondlengths = hypot.(coords[ : , idtoindex[bonds[1, : ]]] - coords[ : , idtoindex[bonds[2, : ]]])
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