abstract type FockBuilder end

struct RestrictedFockBuilder{T} <: FockBuilder
    _ovlp ::Array{T, 2}
    _hcore::Array{T, 2}
    _eri  ::TwoElectronIntegralAO{T}
end

abstract type SCFSolver{T} <: ObjectiveFunction{T} end 
mutable struct RestrictedSCFSolver{T} <: SCFSolver{T}
    _nelec::Integer
    _nocc ::Integer
    _nao  ::Integer
    _nmo  ::Integer

    _guess_method  ::InitGuessMethod
    _fock_builder  ::RestrictedFockBuilder{T}

    _is_converged  ::Bool
    _orb_ene       ::Array{T, 1}
    _orb_coeff     ::Array{T, 2}
    _orb_coeff_inv ::Array{T, 2}
    _density_matrix::Array{T, 2}
    _fock_matrix   ::Array{T, 2} 
end

function get_nao(the_scf::RestrictedSCFSolver)
    return the_scf._nao::Integer
end

function get_nmo(the_scf::RestrictedSCFSolver)
    return the_scf._nmo::Integer
end

function get_occ_index(the_scf::RestrictedRestrictedSCFSolver)
    return (1:the_scf._nocc)::UnitRange{Integer}
end

function get_vir_index(the_scf::RestrictedRestrictedSCFSolver)
    return ((the_scf._nocc+1):(the_scf._nmo))::UnitRange{Integer}
end

function get_eri_value(the_scf::RestrictedRestrictedSCFSolver, lm::Integer, sgm::Integer, mu::Integer, nu::Integer)
    return get_value(the_scf._fock_builder._eri, lm, sgm, mu, nu)
end

function build_orbitals(the_scf::RestrictedRestrictedSCFSolver{T}, fock_matrix::Array{T, 2}) where {T}
    ovlp_matrix = the_scf._fock_builder._ovlp
    vals, vecs  = eigen(fock_matrix, ovlp_matrix)
    return vals, vecs
end

function build_density_matrix(the_scf::RestrictedRestrictedSCFSolver{T}, orb_coeff::Array{T, 2}) where {T}
    nao       = get_nao(the_scf)
    occ_index = get_occ_index(the_scf)

    co        = orb_coeff[:, occ_index]
    co_t      = transpose(x)

    dm        = Array{T}(undef, nao, nao)
    mul!(dm, co, co_t, 2.0, false)
    return dm
end

function build_fock_matrix(the_scf::RestrictedRestrictedSCFSolver{T}, dm::Array{T, 2}) where {T}
    nao   = get_nao(the_scf)
    hcore = the_scf._fock_builder.hcore
    jmat  = zeros(T, nao, nao)
    kmat  = zeros(T, nao, nao)
    dm_tot            = get_dm_tot(the_scf,  dm)
    dm_alpha, dm_beta = get_dm_spin(the_scf, dm)

    for lm in 1:nao
        for sgm in 1:nao
            for mu in 1:nao
                for nu in 1:nao
                    jmat[lm, sgm] += get_eri_value(the_scf, lm, sgm, mu, nu) * dm_tot(nu, mu)
                    kmat[lm, sgm] += get_eri_value(the_scf, lm, mu, sgm, nu) * dm_alpha(nu, mu)
                end
            end
        end
    end

    return hcore + jmat - kmat
end

function update!(the_scf::RestrictedRestrictedSCFSolver{T}, is_converged::Bool, e_tot::T, orb_ene::Array{T, 1}, orb_coeff::Array{T, 2}, density_matrix::Array{T, 2}, fock_matrix::Array{T, 2}) where {T}
    the_scf._is_converged   = is_converged
    the_scf._e_tot          = e_tot
    the_scf._orb_ene        = orb_ene
    the_scf._orb_coeff      = orb_coeff
    the_scf._density_matrix = density_matrix
    the_scf._fock_matrix    = fock_matrix
end

function kernel!(the_scf::SCFSolver{T}; max_iter::Integer=100, tol::Number=1e-8) where {T}
    
    the_optimizer::OptimizationAlgorithm{T} = algorithm_selection(the_scf, tol)
    reset!(the_optimizer)

    iter::Integer      = 0
    is_converged::Bool = false

    while not(is_converged) && iter < max_iter

        orb_coff, orb_ene = build_orbitals(the_optimizer, the_scf)
        density_matrix    = build_density_matrix(the_scf, orb_coff)
        fock_matrix       = build_fock_matrix(the_scf, density_matrix)
        e_tot             = total_energy(the_scf, density_matrix, fock_matrix)

        is_converged      = next_step!(the_optimizer, e_tot)
        update!(
            the_scf, is_converged, e_tot,
            orb_ene,        orb_coff,
            density_matrix, fock_matrix
        )
        iter += 1
    end
end