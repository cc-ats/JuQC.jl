abstract type SCFAlgorithm{T,RealType<:Real} end

mutable struct Roothaan{T,RealType} <: SCFAlgorithm{T,RealType}
    _is_initialized::Bool
    _the_scf::Union{Nothing, SCFSolver{T,RealType}}
end

function Roothaan(;T=Float64)
    return Roothaan{T,real(T)}(false, nothing)
end

function init_scf_algo!(
    scf_algo::SCFAlgorithm{T,RealType},
    the_scf::RestrictedSCFSolver{T,RealType},
    density_matrix::Hermitian{T,Array{T,2}},
    fock_matrix   ::Hermitian{T,Array{T,2}},
    orb_ene       ::Array{T,1},
    orb_coeff     ::Array{T,2},
    ) where {T,RealType}
    if isnothing(scf_algo._the_scf) && not(scf_algo._is_initialized)
        the_scf.density_matrix   = density_matrix
        the_scf.fock_matrix      = fock_matrix
        the_scf.orb_ene          = orb_ene
        the_scf.orb_coeff        = orb_coeff
        scf_algo._the_scf        = the_scf
        scf_algo._is_initialized = true
    else
        error("The SCF algorithm is is initialized.")
    end
end

function next_step!(scf_algo::Roothaan{T,RealType}) where {T,RealType}
    the_scf           = scf_algo._the_scf
    orb_ene           = the_scf.orb_ene
    orb_coeff         = the_scf.orb_coeff
    density_matrix    = build_density_matrix(the_scf,   orb_coeff)
    fock_matrix       = build_fock_matrix(the_scf, density_matrix)

    next_orb_ene, next_orb_coeff = build_orbitals(the_scf, fock_matrix)
    the_scf.orb_ene              = next_orb_ene
    the_scf.orb_coeff            = next_orb_coeff
    return orb_ene, orb_coeff, fock_matrix, density_matrix
end

# struct DIIS{T} <: SCFAlgorithm{T,RealType}
#     _is_initialized::Bool
#     _num_sub_space::Integer
#     _the_scf::Union{Nothing, ObjectiveFunction{T}}
# end

# function DIIS(;num_sub_space=8, T=Float64)
#     return DIIS{T}(false, num_sub_space, nothing)
# end

function kernel!(
    the_scf::RestrictedSCFSolver{T,RealType};
    max_iter::Integer=100, tol::Number=1e-8,
    scf_algo::SCFAlgorithm{T,RealType}=Roothaan(T=T)
    ) where {T,RealType}
    
    if the_scf.is_init_guess
        init_dm              = the_scf.density_matrix
        init_fock            = build_fock_matrix(the_scf,   init_dm)
        init_norm, init_grad = get_grad(the_scf, init_fock, init_dm)
    else
        hcore                = the_scf._fock_builder._hcore
        init_orb_ene, init_orb_coeff  = build_orbitals(the_scf, hcore)
        init_dm              = build_density_matrix(the_scf, init_orb_coeff)
        init_fock            = build_fock_matrix(the_scf,   init_dm)
        init_norm, init_grad = get_grad(the_scf, init_fock, init_dm)
    end

    e_elec = get_e_elec(the_scf, init_fock, init_dm)
    e_tot  = e_elec  + the_scf._e_nuc
    pre_e_elec = zero(T)
    cur_e_elec = e_elec

    if not(scf_algo._is_initialized)
        init_scf_algo!(
        scf_algo, the_scf, 
        init_fock, init_dm,
        init_orb_ene, init_orb_coeff
        )
    end

    iter::Integer      = 1
    is_converged::Bool = false

    while not(is_converged) && iter < max_iter
        orb_ene, orb_coeff, fock_matrix, density_matrix  = next_step!(scf_algo)

        grad_norm, grad     = get_grad(the_scf, fock_matrix, density_matrix)
        e_elec              = get_e_elec(the_scf, fock_matrix, density_matrix)

        pre_e_elec   = cur_e_elec
        cur_e_elec   = e_elec
        ene_err      = abs(cur_e_elec - pre_e_elec)
        grad_err     = grad_norm
        err          = max(grad_norm, ene_err)
        is_converged = err < tol

        e_tot = e_elec  + the_scf._e_nuc
        @printf("%5d%19.10f%14.2e\n", iter, e_tot, err)
        iter += 1
    end
    println("Number of Fock Build = ", the_scf.num_fock_build)
    println("Time   of Fock Build = ", the_scf.time_fock_build)
    
    return RestrictedSCFResult{T,RealType}(
        is_converged, e_tot, e_elec,
        orb_ene, orb_coeff, grad,
        density_matrix, fock_matrix
        )
end