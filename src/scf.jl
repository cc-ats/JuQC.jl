mutable struct FockBuilder{T}
    _nao  ::Integer
    _ovlp ::Hermitian{T,Array{T,2}}
    _hcore::Hermitian{T,Array{T,2}}
    _eri  ::TwoElectronIntegralAO{T}
end

abstract type SCFSolver{T}            <: ObjectiveFunction{T} end 
mutable struct RestrictedSCFSolver{T} <: SCFSolver{T}
    # input parameters
    _nao           ::Integer
    _nmo           ::Integer
    _nocc          ::Integer
    _e_nuc         ::Real
    _fock_builder  ::FockBuilder{T}

    # parameters
    is_init_guess ::Bool
    is_converged  ::Bool
    num_fock_build::Integer

    # results
    energy_tot    ::Real
    energy_elec   ::Real

    orb_ene       ::Array{T,1}
    orb_coff      ::Array{T,2}
    grad          ::Array{T,2}

    density_matrix::Hermitian{T,Array{T,2}}
    fock_matrix   ::Hermitian{T,Array{T,2}}
end

mutable struct UnrestrictedSCFSolver{T} <: SCFSolver{T}
    _nao          ::Integer
    _nmo          ::Integer
    _nocc         ::Tuple{Integer,Integer}
    _e_nuc        ::Real
    _fock_builder ::FockBuilder{T}

    is_init_guess ::Bool
    is_converged  ::Bool
    num_fock_build::Integer

    energy_tot    ::Real
    energy_elec   ::Real

    orb_ene       ::Tuple{Array{T,1},Array{T,1}} 
    orb_coff      ::Tuple{Array{T,2},Array{T,2}}
    grad          ::Tuple{Array{T,2},Array{T,2}}

    density_matrix::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
    fock_matrix   ::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
end

function build_scf_solver(
    nao::Integer, e_nuc::Real,
    nelec::Tuple{Integer,Integer},
    ovlp::Hermitian{T,Array{T,2}},
    hcore::Hermitian{T,Array{T,2}},
    eri::TwoElectronIntegralAO{T};
    is_restricted::Bool=true
    ) where {T}
    
    if is_restricted
        if not(nelec[1]==nelec[2])
            error("RHF must be closed shell")
        end
        fock_builder  = FockBuilder{T}(nao, ovlp, hcore, eri)
        the_scf       = RestrictedSCFSolver{T}(
            nao, nao, nelec[1], e_nuc, fock_builder,
            false, false, 0, 0.0, 0.0,
            zeros(T, nao), zeros(T, nao, nao), zeros(T, nao, nao),
            Hermitian(zeros(T, nao, nao)), Hermitian(zeros(T, nao, nao))
            )
        return the_scf
    else
        error("Not implemented!")
    end
end

function get_nao(the_scf::RestrictedSCFSolver)
    return the_scf._nao::Integer
end

function get_nmo(the_scf::RestrictedSCFSolver)
    return the_scf._nmo::Integer
end

function get_nocc(the_scf::RestrictedSCFSolver)
    return the_scf._nocc
end

function get_nvir(the_scf::RestrictedSCFSolver)
    return the_scf._nmo - the_scf._nocc
end

function get_occ_index(the_scf::RestrictedSCFSolver)
    return (1:the_scf._nocc)
end

function get_vir_index(the_scf::RestrictedSCFSolver)
    return ((the_scf._nocc+1):(the_scf._nmo))
end

function get_density_matrix_tot(the_scf::RestrictedSCFSolver{T}, p::Hermitian{T,Array{T,2}}) where {T}
    return p::Hermitian{T,Array{T,2}}
end

function get_density_matrix_spin(the_scf::RestrictedSCFSolver{T}, p::Hermitian{T,Array{T,2}}) where {T}
    return ((p/2, p/2))::Tuple{Hermitian{T,Array{T,2}}, Hermitian{T,Array{T,2}}}
end

function get_e_elec(the_scf::RestrictedSCFSolver{T}, f::Hermitian{T,Array{T,2}}, p::Hermitian{T,Array{T,2}}) where {T}
    h      = the_scf._fock_builder._hcore
    tmp    = f + h
    return real(dot(p, tmp)/2)::Real
end

function get_grad(the_scf::RestrictedSCFSolver{T}, f::Hermitian{T,Array{T,2}}, p::Hermitian{T,Array{T,2}}) where {T}
    nao  = get_nao(the_scf)
    s    = the_scf._fock_builder._ovlp
    fps  = f * p * s
    g    = fps - transpose(fps)
    n    = norm(g)/nao
    return n::Real, g::Array{T,2}
end

function get_value(the_scf::RestrictedSCFSolver{T}) where {T}
    return the_scf.energy_elec::Real
end

function get_grad(the_scf::RestrictedSCFSolver{T}) where {T}
    return the_scf.grad::Array{T, 1}
end

function build_orbitals(the_scf::RestrictedSCFSolver{T}, fock_matrix::Hermitian{T,Array{T,2}}) where {T}
    ovlp_matrix = the_scf._fock_builder._ovlp
    vals, vecs  = eigen(fock_matrix::Hermitian{T,Array{T,2}}, ovlp_matrix::Hermitian{T,Array{T,2}})
    return vals::Array{T,1}, vecs::Array{T,2}
end

function build_density_matrix(the_scf::RestrictedSCFSolver{T}, orb_coeff::Array{T, 2}) where {T}
    nao       = get_nao(the_scf)
    occ_index = get_occ_index(the_scf)

    co        = orb_coeff[:, occ_index]
    co_t      = transpose(co)

    dm        = Array{T}(undef, nao, nao)
    mul!(dm, co, co_t, 2.0, 0.0)

    return Hermitian(dm::Array{T,2})
end

function build_fock_matrix(
    the_scf::RestrictedSCFSolver{T},
    dm::Hermitian{T,Array{T,2}}
    ) where {T}
    nao    = the_scf._nao
    hcore  = the_scf._fock_builder._hcore

    dm_tot            = get_density_matrix_tot(the_scf::RestrictedSCFSolver{T},  dm)
    dm_alpha, dm_beta = get_density_matrix_spin(the_scf::RestrictedSCFSolver{T}, dm)

    jmat   = zeros(T, nao, nao)
    kmat   = zeros(T, nao, nao)

    for lm in 1:nao
        for sgm in 1:nao
            for mu in 1:nao
                for nu in 1:nao
                    jmat[lm, sgm] += get_value(
                        the_scf._fock_builder._eri, lm, sgm, mu, nu
                        ) * dm_tot[nu, mu]

                    kmat[lm, sgm] += get_value(
                        the_scf._fock_builder._eri, lm, mu, sgm, nu
                        ) * dm_alpha[nu, mu]
                end
            end
        end
    end

    fock = hcore + jmat - kmat
    the_scf.num_fock_build += 1

    return Hermitian(fock::Array{T,2})
end

function set_init_guess!(
    the_scf::RestrictedSCFSolver{T}, 
    is_orb_coeff::Bool, tmp::Array{T,2}
    ) where {T}
    nao = the_scf._nao

    if not(sizeof(tmp) == (nao, nao))
        error("Wrong init guess dimension!!!")
    end

    if is_orb_coeff
        the_scf.is_init_guess  = true
        the_scf.density_matrix = build_density_matrix(the_scf, tmp_array)
    else
        the_scf.is_init_guess  = true
        the_scf.density_matrix = Hermitian(tmp_array)
    end
end

function update!(
    the_scf       ::RestrictedSCFSolver{T},
    is_converged  ::Bool,
    energy_tot    ::Real,
    energy_elec   ::Real,
    orb_ene       ::Array{T,1},
    orb_coff      ::Array{T,2},
    grad          ::Array{T,2},
    density_matrix::Hermitian{T,Array{T,2}},
    fock_matrix   ::Hermitian{T,Array{T,2}},
    ) where {T}

    the_scf.is_converged   = is_converged
    the_scf.energy_tot     = energy_tot
    the_scf.energy_elec    = energy_elec
    the_scf.orb_ene        = orb_ene
    the_scf.orb_coff       = orb_coff
    the_scf.grad           = grad
    the_scf.density_matrix = density_matrix
    the_scf.fock_matrix    = fock_matrix
end

function kernel!(the_scf::RestrictedSCFSolver{T};
    max_iter::Integer=100, tol::Number=1e-8,
    the_opt::OptimizationAlgorithm{T}=RoothaanOptimizer{T}(false, nothing)
    ) where {T}
    
    if the_scf.is_init_guess
        init_dm           = the_scf.density_matrix
        fock_matrix       = build_fock_matrix(the_scf,  init_dm)
        orb_ene, orb_coff = build_orbitals(the_scf, fock_matrix)
        density_matrix    = build_density_matrix(the_scf, orb_coff)
    else
        hmat               = the_scf._fock_builder._hcore
        orb_ene, orb_coff  = build_orbitals(the_scf, hmat)
        density_matrix     = build_density_matrix(the_scf, orb_coff)
        fock_matrix        = build_fock_matrix(the_scf, density_matrix)
    end

    e_elec  = get_e_elec(the_scf, fock_matrix, density_matrix)
    pre_e_elec = zero(T)
    cur_e_elec = e_elec

    if not(the_opt._is_initialized)
        init_opt_algo!(the_opt, the_scf)
    end

    iter::Integer      = 1
    is_converged::Bool = false

    while not(is_converged) && iter < max_iter
        is_opt_converged  = next_step!(the_opt)
        if (iter > 1)
            density_matrix    = get_density_matrix(the_opt)
        end
        fock_matrix       = build_fock_matrix(the_scf, density_matrix)
        orb_ene, orb_coff = build_orbitals(the_scf, fock_matrix)
        grad_norm, grad   = get_grad(the_scf,   fock_matrix, density_matrix)

        density_matrix    = build_density_matrix(the_scf, orb_coff)
        e_elec            = get_e_elec(the_scf, fock_matrix, density_matrix)

        pre_e_elec   = cur_e_elec
        cur_e_elec   = e_elec
        ene_err      = abs(cur_e_elec - pre_e_elec)
        grad_err     = grad_norm
        err          = max(grad_norm, ene_err)
        is_converged = err < tol

        e_tot = e_elec  + the_scf._e_nuc
        update!(the_scf, is_converged, e_tot, e_elec, orb_ene, orb_coff, grad, density_matrix, fock_matrix)

        @printf("%5d%19.10f%14.2e\n", iter, e_tot, err)
        iter += 1
    end
    println("NFock = ", the_scf.num_fock_build)
end