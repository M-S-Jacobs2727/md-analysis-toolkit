using Statistics

import LammpsFiles

"""
    We assume that each particle represents a monomer with unit
    mass and there are no other atoms in the dump files.
"""
function rg2(filenames...; sorted=false)
    frame = LammpsFiles.read_dump(filenames[1])
    mol_col = findfirst(x->x=="mol", frame.properties)
    coord_cols = [findfirst(x->x==s, frame.properties) for s in ["xu", "yu", "zu"]]
    
    molecules = frame.atoms[mol_col, : ]
    nmolecules = length(unique(molecules))
    natoms_per_mol = [count(x->x==m, molecules) for m in unique(molecules)]
    inds = [1; cumsum(natoms_per_mol)...]
    # same as `frame.natoms / nmolecules` for no solvent and monodisperse
    
    rg2 = 0.0
    for filename in filenames
        frame = LammpsFiles.read_dump(filename)
        molecules = frame.atoms[mol_col, : ]
        coords = frame.atoms[coord_cols, : ]
        if !sorted
            inds = sortperm(molecules)
            coords = coords[ : , inds]
        end
        coords = [coords[ : , inds[i]:inds[i+1]] for i=1:nmolecules-1]

        centers_of_mass = reshape.(sum.(coords, dims=2), 3)
        for (mol, com, dp) in zip(coords, centers_of_mass, natoms_per_mol)
            mol_offset = mol .- com
            rg2 += sum(mol_offset .* mol_offset) / dp
        end
    end

    rg2 /= nmolecules * length(filenames)
end