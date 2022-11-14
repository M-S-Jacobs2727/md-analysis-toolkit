import LammpsFiles

"""
    square_Rg(filenames...; sorted=false, masses=nothing)
    
Compute the mean square radius of gyration of a system of polymer chains for
each snapshot, returning a vector of length equal to the number of files
passed.
"""
function square_Rg(filenames...; sorted=false, masses=nothing)
    # Read basic info
    frame = LammpsFiles.read_dump(filenames[1])
    mol_col = findfirst(x->x=="mol", frame.properties)
    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["xu", "yu", "zu"]]
    
    # Index molecules
    mol_ids = frame.atoms[mol_col, : ]
    molecules = sort!(unique(mol_ids))
    nmolecules = length(molecules)
    natoms_per_mol = [count(x->x==m, mol_ids) for m in molecules]
    mol_indices = [1; (cumsum(natoms_per_mol)[1:end-1] .+ 1)...]
    
    # Get square radius of gyration for each frame
    rg2 = zeros(Float64, length(filenames))
    for (i, filename) in enumerate(filenames)
        frame = LammpsFiles.read_dump(filename)
        molecules = frame.atoms[mol_col, : ]
        coords = frame.atoms[coord_cols, : ]
        if !sorted
            coords = coords[ : , sortperm(molecules)]
        end
        coords = [coords[ : , mol_indices[i]:mol_indices[i+1]] for i=1:nmolecules-1]

        # TODO: accept masses as an argument for correct COM
        centers_of_mass = reshape.(sum.(coords, dims=2), 3)
        for (mol, com, dp) in zip(coords, centers_of_mass, natoms_per_mol)
            mol_offset = mol .- com
            rg2[i] += sum(mol_offset .* mol_offset) / dp
        end
    end

    rg2 ./= nmolecules
end