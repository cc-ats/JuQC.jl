mutable struct FockBuilder{T}
    _nao  ::Integer
    _ovlp ::Hermitian{T,Array{T,2}}
    _hcore::Hermitian{T,Array{T,2}}
    _eri  ::TwoElectronIntegralAO{T}
end

abstract type  SCFSolver{T} end 
mutable struct RestrictedSCFSolver{T}  <: SCFSolver{T}
    # input parameters
    _nao           ::Integer
    _nmo           ::Integer
    _nocc          ::Integer
    _e_nuc         ::real(T)
    _fock_builder  ::FockBuilder{T}

    # parameters
    is_init_guess ::Bool
    is_converged  ::Bool

    # results
    num_fock_build ::Integer
    time_fock_build::real(T)

    energy_tot    ::real(T)
    energy_elec   ::real(T)

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
    _e_nuc        ::real(T)
    _fock_builder ::FockBuilder{T}

    is_init_guess ::Bool
    is_converged  ::Bool

    num_fock_build ::Integer
    time_fock_build::real(T)

    energy_tot    ::real(T)
    energy_elec   ::real(T)

    orb_ene       ::Tuple{Array{T,1},Array{T,1}} 
    orb_coff      ::Tuple{Array{T,2},Array{T,2}}
    grad          ::Tuple{Array{T,2},Array{T,2}}

    density_matrix::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
    fock_matrix   ::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
end

abstract type  SCFResult{T} end 
mutable struct RestrictedSCFResult{T} <: SCFResult{T}
    is_converged  ::Bool

    # results
    energy_tot    ::real(T)
    energy_elec   ::real(T)

    orb_ene       ::Array{T,1}
    orb_coff      ::Array{T,2}
    grad          ::Array{T,2}

    density_matrix::Hermitian{T,Array{T,2}}
    fock_matrix   ::Hermitian{T,Array{T,2}}
end

struct UnrestrictedSCFResult{T}       <: SCFResult{T}
    is_converged  ::Bool

    # results
    energy_tot    ::real(T)
    energy_elec   ::real(T)

    orb_ene       ::Tuple{Array{T,1},Array{T,1}} 
    orb_coff      ::Tuple{Array{T,2},Array{T,2}}
    grad          ::Tuple{Array{T,2},Array{T,2}}

    density_matrix::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
    fock_matrix   ::Tuple{Hermitian{T,Array{T,2}},Hermitian{T,Array{T,2}}}
end

function build_scf_solver(
    nao::Integer, e_nuc::real(T),
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
            false, false, 0, 0.0, 0.0, 0.0,
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
    return real(dot(p, f + h)/2)::real(T)
end

function get_grad(the_scf::RestrictedSCFSolver{T}, f::Hermitian{T,Array{T,2}}, p::Hermitian{T,Array{T,2}}) where {T}
    s    = the_scf._fock_builder._ovlp
    fps  = f * p * s
    g    = fps - transpose(fps)
    n    = norm(g)
    return n::real(T), g::Array{T,2}
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

    t0 = time()
    jmat   = zeros(T, nao, nao)
    kmat   = zeros(T, nao, nao)

    for lm in 1:nao
        for sgm in 1:nao
            for mu in 1:nao
                for nu in 1:nao
                    @inbounds jmat[lm, sgm] += get_value(
                        the_scf._fock_builder._eri, lm, sgm, mu, nu
                        ) * dm_tot[nu, mu]

                    @inbounds kmat[lm, sgm] += get_value(
                        the_scf._fock_builder._eri, lm, mu, sgm, nu
                        ) * dm_alpha[nu, mu]
                end
            end
        end
    end

    fock = hcore + jmat - kmat
    the_scf.num_fock_build  += 1
    the_scf.time_fock_build += time() - t0

    return Hermitian(fock::Array{T,2})
end

function set_init_guess!(
    the_scf::RestrictedSCFSolver{T}, 
    tmp::Array{T,2};
    is_orb_coeff::Bool=true, is_density_matrix::Bool=false
    ) where {T}
    nao = get_nao(the_scf)
    nmo = get_nmo(the_scf)

    if is_orb_coeff && not(is_density_matrix)
        if not(sizeof(tmp) == (nao, nmo))
            error("Wrong init guess dimension!!!")
        end
        the_scf.is_init_guess  = true
        the_scf.density_matrix = build_density_matrix(the_scf, tmp)
    else if not(is_orb_coeff) && is_density_matrix
        if not(sizeof(tmp) == (nao, nao))
            error("Wrong init guess dimension!!!")
        end
        the_scf.is_init_guess  = true
        the_scf.density_matrix = Hermitian(tmp)
    else
        error("Unknown init guess type!!!")
    end
end

function update!(
    the_scf       ::RestrictedSCFSolver{T},
    is_converged  ::Bool,
    energy_tot    ::real(T),
    energy_elec   ::real(T),
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

function kernel!(
    the_scf::RestrictedSCFSolver{T};
    max_iter::Integer=100, tol::Number=1e-8,
    scf_algo::SCFAlgorithm=Roothaan(T=T)
    ) where {T}
    
    if the_scf.is_init_guess
        init_dm              = the_scf.density_matrix
        init_fock            = build_fock_matrix(the_scf,   init_dm)
        init_norm, init_grad = get_grad(the_scf, init_fock, init_dm)
    else
        hcore                = the_scf._fock_builder._hcore
        init_orb_ene, init_orb_coff  = build_orbitals(the_scf, hcore)
        init_dm              = build_density_matrix(the_scf, init_orb_coff)
        init_fock            = build_fock_matrix(the_scf,   init_dm)
        init_norm, init_grad = get_grad(the_scf, init_fock, init_dm)
    end

    e_elec     = get_e_elec(the_scf, fock_matrix, density_matrix)
    pre_e_elec = zero(T)
    cur_e_elec = e_elec

    if not(scf_algo._is_initialized)
        init_scf_algo!(scf_algo, the_scf, init_fock, init_dm)
    end

    iter::Integer      = 1
    is_converged::Bool = false

    while not(is_converged) && iter < max_iter
        is_scf_algo_converged  = next_step!(scf_algo)
        if (iter > 1)
            density_matrix     = get_density_matrix(scf_algo)
        end
        fock_matrix       = build_fock_matrix(the_scf, density_matrix)
        orb_ene, orb_coff = build_orbitals(the_scf, fock_matrix)
        grad_norm, grad   = get_grad(the_scf,   fock_matrix, density_matrix)

        new_density_matrix  = build_density_matrix(the_scf, orb_coff)
        e_elec              = get_e_elec(the_scf, fock_matrix, density_matrix)

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
    println("Number of Fock Build = ", the_scf._fock_builder.num_fock_build)
    println("Time   of Fock Build = ", the_scf._fock_builder.time_fock_build)
    
    return RestrictedSCFResult{T}(is_converged, e_tot, e_elec, orb_ene, orb_coff, grad, density_matrix, fock_matrix, )
end