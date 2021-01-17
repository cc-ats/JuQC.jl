abstract type SCFAlgorithm{T,RealType<:Real} end

mutable struct Roothaan{T,RealType} <: SCFAlgorithm{T,RealType}
    _is_initialized::Bool
    _the_scf::Union{Nothing, SCFSolver{T,RealType}}
end

function Roothaan(;T=Float64)
    return Roothaan{T,real(T)}(false, nothing)
end

function init_scf_algo!(
    scf_algo::Roothaan{T,RealType},
    the_scf::RestrictedSCFSolver{T,RealType},
    density_matrix::Hermitian{T,Array{T,2}},
    fock_matrix   ::Hermitian{T,Array{T,2}},
    orb_ene       ::Array{RealType,1},
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

mutable struct DIIS{T,RealType} <: SCFAlgorithm{T,RealType}
    _is_initialized::Bool
    _the_scf::Union{Nothing, SCFSolver{T,RealType}}

    _size_subspace::Integer
    _vec_err::Union{Nothing, Vector{Array{T, 2}},Vector{Tuple{Array{T, 2},Array{T, 2}}}}
    _vec_fock::Union{Nothing, Vector{Hermitian{T,Array{T,2}}},Vector{Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}}}
end

function DIIS(;size_subspace=8, T=Float64)
    return DIIS{T,real(T)}(false, nothing, size_subspace, nothing, nothing)
end

function init_scf_algo!(
    scf_algo::DIIS{T,RealType},
    the_scf::RestrictedSCFSolver{T,RealType},
    density_matrix::Hermitian{T,Array{T,2}},
    fock_matrix   ::Hermitian{T,Array{T,2}},
    orb_ene       ::Array{RealType,1},
    orb_coeff     ::Array{T,2},
    ) where {T,RealType}
    if isnothing(scf_algo._the_scf) && not(scf_algo._is_initialized)
        the_scf.density_matrix   = density_matrix
        the_scf.fock_matrix      = fock_matrix
        the_scf.orb_ene          = orb_ene
        the_scf.orb_coeff        = orb_coeff
        scf_algo._the_scf        = the_scf
        scf_algo._vec_err        = Vector{Array{T,2}}()
        scf_algo._vec_fock       = Vector{Hermitian{T,Array{T,2}}}()

        norm, grad = get_grad(the_scf, fock_matrix, density_matrix)
        push!(scf_algo._vec_err,         grad)
        push!(scf_algo._vec_fock, fock_matrix)
        
        scf_algo._is_initialized = true
    else
        error("The SCF algorithm is is initialized.")
    end
end

function build_diis_fock(vec_err::Vector{Array{T,2}}, vec_fock::Vector{Hermitian{T,Array{T,2}}}) where {T}
    size_vec_err = size(vec_err, 1)
    if size_vec_err > 2
        a_vec        = zeros(T, size_vec_err+1)
        b_mat        = fill(-1.0, (size_vec_err+1, size_vec_err+1))

        a_vec[size_vec_err+1] = -1.0
        for a in 1:size_vec_err
            for b in 1:size_vec_err
                b_mat[a,b] = dot(vec_err[a], vec_err[b])
            end
        end
        
        b_inv = inv(Hermitian(b_mat))
        coeff = b_inv * a_vec
        diis_fock = vec_fock[1] * coeff[1]

        for a in 2:size_vec_err
            diis_fock += vec_fock[a] * coeff[a]
        end
        return Hermitian(diis_fock)
    else
        return vec_fock[size_vec_err]
    end
end

function build_diis_fock(vec_err::Vector{Array{T,2}}, vec_fock::Vector{Hermitian{T,Array{T,2}}}) where {T}
    size_vec_err = size(vec_err, 1)
    if size_vec_err > 2
        a_vec        = zeros(T, size_vec_err+1)
        b_mat        = fill(-1.0, (size_vec_err+1, size_vec_err+1))

        a_vec[size_vec_err+1] = -1.0
        for a in 1:size_vec_err
            for b in 1:size_vec_err
                b_mat[a,b] = dot(vec_err[a], vec_err[b])
            end
        end
        
        b_inv = inv(Hermitian(b_mat))
        coeff = b_inv * a_vec
        diis_fock = vec_fock[1] * coeff[1]

        for a in 2:size_vec_err
            diis_fock += vec_fock[a] * coeff[a]
        end
        return Hermitian(diis_fock)
    else
        return vec_fock[size_vec_err]
    end
end

function next_step!(scf_algo::DIIS{T,RealType}) where {T,RealType}
    the_scf           = scf_algo._the_scf
    orb_ene           = the_scf.orb_ene
    orb_coeff         = the_scf.orb_coeff
    density_matrix    = build_density_matrix(the_scf,   orb_coeff)
    fock_matrix       = build_fock_matrix(the_scf, density_matrix)
    norm, grad        = get_grad(the_scf, fock_matrix, density_matrix)

    vec_fock = scf_algo._vec_fock
    vec_err  = scf_algo._vec_err

    if size(vec_fock, 1) < scf_algo._size_subspace
        push!(vec_fock, fock_matrix)
        push!(vec_err,  grad)
    else
        deleteat!(vec_fock, 1)
        deleteat!(vec_err,  1)
        push!(vec_fock, fock_matrix)
        push!(vec_err,  grad)
    end

    diis_fock = build_diis_fock(vec_err, vec_fock)

    next_orb_ene, next_orb_coeff = build_orbitals(the_scf, diis_fock)
    the_scf.orb_ene              = next_orb_ene
    the_scf.orb_coeff            = next_orb_coeff
    return orb_ene, orb_coeff, fock_matrix, density_matrix
end

function kernel!(
    the_scf::RestrictedSCFSolver{T,RealType};
    max_iter::Integer=100, tol::Number=1e-8,
    scf_algo::SCFAlgorithm{T,RealType}=Roothaan(T=T)
    ) where {T,RealType}
    
    if the_scf.is_init_guess
        init_dm              = the_scf.density_matrix
        init_fock            = build_fock_matrix(the_scf,   init_dm)
    else
        init_fock            = the_scf._fock_builder._hcore
    end

    orb_ene, orb_coeff  = build_orbitals(the_scf,         init_fock)
    density_matrix      = build_density_matrix(the_scf,   orb_coeff)
    fock_matrix         = build_fock_matrix(the_scf, density_matrix)
    norm, grad          = get_grad(the_scf, fock_matrix, density_matrix)

    e_elec = get_e_elec(the_scf, fock_matrix, density_matrix)
    e_tot  = e_elec  + the_scf._e_nuc
    pre_e_elec = zero(T)
    cur_e_elec = e_elec

    init_scf_algo!(scf_algo, the_scf, fock_matrix, density_matrix, orb_ene, orb_coeff)

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