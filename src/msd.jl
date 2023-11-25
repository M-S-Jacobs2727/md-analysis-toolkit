import MPI

function msd(filenames::Vector{String}, maxlag::Integer=100)
    MPI.Init()
    comm = MPI.COMM_WORLD
    me = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm) 

    max_atoms_per_proc = dim = natoms = ntimesteps = 0
    if me == 0
        dim, natoms, ntimesteps = size(trajectory)
        trajectory = permutedims(trajectory, (1, 3, 2))
        traj_type = eltype(trajectory)
        # permuting the dims allows for the parallelization by atoms,
        # as they are independent of each other

        max_atoms_per_proc = Int(ceil(natoms / nprocs))
    end
    max_atoms_per_proc = MPI.bcast(max_atoms_per_proc, comm)
    dim = MPI.bcast(dim, comm)
    ntimesteps = MPI.bcast(ntimesteps, comm)
    natoms = MPI.bcast(natoms, comm)
    my_trajectory = Array{Float32}(undef, dim, ntimesteps, max_atoms_per_proc)
    q, r = divrem(natoms, nprocs)
    natoms_each = [i <= r ? q+1 : q for i = 1:n]
    my_natoms = natoms_each[me+1]

    if me == 0
        vbuf = MPI.VBuffer(trajectory, dim*ntimesteps*natoms_each)
    else
        vbuf = MPI.VBuffer(nothing)
    end
    my_trajectory = MPI.Scatterv!(vbuf, zeros{eltype(tra)})
    

    msd = zeros(maxlag)
    for a=1:natoms
        for lag=1:maxlag
            msd[lag] += sum(trajectory[:, lag+1:end, a] - trajectory[:, 1:end-lag, a])
        end
    end

    msd /= natoms*(ntimesteps - (1:maxlag))

    return msd
end