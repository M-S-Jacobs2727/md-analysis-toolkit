using Statistics

import LammpsFiles

"""
    rdf(filenames...; bin_width=0.05, by_type=false)

    Compute the radial distribution function for a collection of LAMMPS dump
    files (`filenames`). 

    The width of the radial bins is set to `bin_width=0.05` by default.

    If `by_type` is false (default) then the overall RDF is computed. If 
    true, the RDF between each type is computed.

"""
function rdf(filenames...; bin_width=0.05, by_type=false)
    
end