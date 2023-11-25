using Statistics
import LammpsFiles


"""
Fd(q, t) = 1/N <sum_{j=1}^N sum_{k=1}^N exp{iq . (r_j(t0+t) - r_k(t0))}>_t0
         = 1/N <sum_{j=1}^N sum_{k=1}^N cos(q . (r_j(t0+t) - r_k(t0))) + 
                                    + i sin(q . (r_j(t0+t) - r_k(t0)))>_t0   
         (above sine part goes to 0 on averaging q over angles for homogenous systems)
         = 1/N <(sum_{j=1}^N cos(q . r_j(t0+t))) (sum_{k=1}^N cos(q . r_k(t0))) -
              - (sum_{j=1}^N sin(q . r_j(t0+t))) (sum_{k=1}^N sin(q . r_k(t0)))>_t0
        
        So, for each frame, compute sum_{j=1}^N cos(q . r_j(t)) and sum_{j=1}^N sin(q . r_j(t)),
        then compute the above correlation for various lag times.
            
Fs(q, t) = 1/N <sum_{j=1}^N exp{iq . (r_j(t0+t) - r_j(t0))}>_t0
         = 1/N <sum_{j=1}^N cos(q . (r_j(t0+t) - r_j(t0))) + 
                        + i sin(q . (r_j(t0+t) - r_j(t0)))>_t0   
         (above sine part goes to 0 on averaging q over angles for homogenous systems?)
         = 1/N <sum_{j=1}^N cos(q . r_j(t0+t))cos(q . r_j(t0)) + 
                          + sin(q . r_j(t0+t))sin(q . r_j(t0))>_t0

The self part (second set of equations above) can be simplified in a similar manner to
the msd function (msd.jl).

Scattering vector q is defined w.r.t. the box length, 
with each dimension i=x,y,z: q_i = 2 n pi / L_i

We assume that filenames refers to a list of LAMMPS dump files, in order, evenly spaced.
"""
function intermediateScatteringFunc(trajectory::Array{<:Real}; maxlag::Integer=100, framedt::Real=0.02)

    dim, natoms, ntimesteps = size(trajectory)
    

end