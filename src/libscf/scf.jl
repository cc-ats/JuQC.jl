mutable struct SelfConsistentField
    # this part are the inputs
    sys::QuantumChemistrySystem
    conv_tol::Float64
    max_cycle::Integer
    init_density_guess::Array{Float64,2}
    orbtype::String
    iprint::Integer

    # OneEMtrx oneE; //!< OneEMtrx class used for obtaining H_core and NOrb
    # BasisSet s1; //!< BasisSet class (used for AOInts calls)
    # ShlPrs s2; //!< ShlPrs class (used for AOInts calls)
    # this part are the intermediate results
    nuclear_energy::Float64
    nbas::Integer
    nao::Integer
    nmo::Integer
    overlp::Array{Float64,2}
    hcore::Array{Float64,2}
    eri::Array{Float64,4}
    
    # this part are the final results
    is_converged::Bool
    scf_energy::Float64
    mo_occ::Array{Integer,1}
    mo_energy::Array{Float64,1}
    mo_coeff::Array{Float64,2}
    function SelfConsistentField(
        sys::QuantumChemistrySystem,
        conv_tol::Float64, max_cycle::Integer,
        init_density_guess::Array{Float64,2};
        orbtype::String="RHF",
        iprint::Integer=3
        )
        
        nuclear_energy = get_nuclear_energy(sys)
        nao, nmo   = get_dimension(sys)
        nbas       = nao
        overlp = get_overlap_matrix(sys)
        hcore  = get_hcore_matrix(sys)
        eri    = get_eri_array(sys)
        
        scf_energy    = 0.0
        mo_occ        = Array{Integer,1}(undef, nmo)
        mo_energy     = Array{Float64,1}(undef, nmo)
        mo_coeff      = Array{Float64,2}(undef, nao, nmo)
        
        if iprint >= 1
            println("")
            println("the SCF is initialized with the ", typeof(sys), " instance.")
            @printf("max_cycle = % 4d, the conv_tol = % 4.1e\n", iprint, conv_tol)
            @printf("nao       = % 4d, nmo          = % 4d\n", nao, nmo)
            @printf("orbtype = %s\n", orbtype)
            @printf("iprint  = %d\n", iprint)
        end

        if iprint >= 2
            @printf("nuclear_energy = % 12.4f\n", nuclear_energy)
            println("")
        end

        new(
            sys, conv_tol, max_cycle, init_density_guess,  # this part are the inputs
            orbtype, iprint,
            nuclear_energy, nbas, nao, nmo,                # this part are the intermediate results
            overlp, hcore, eri,
            false, scf_energy, mo_occ, mo_energy, mo_coeff # this part are the final results
            )
    end
end

function get_jk(the_scf::SelfConsistentField, rdm1::Array{Float64,2})
    g = Array{Float64,2}(undef, the_scf.nao, the_scf.nao)
    if the_scf.orbtype == "RHF"
        for lm=1:the_scf.nao, sg=1:the_scf.nao, mu=1:the_scf.nao, nu=1:the_scf.nao
            g[lm,sg] = g[sg,lm] += (the_scf.eri[lm,mu,sg,nu]-the_scf.eri[lm,mu,sg,nu]/2) * rdm1[mu,nu]
        end
    elseif the_scf.orbtype == "UHF"
        error("UHF not implemented yet!")
    end
    return (g+the_scf.hcore)::Array{Float64,2}
end

function get_fock(the_scf::SelfConsistentField, rdm1::Array{Float64,2})
    g = get_jk(the_scf, rdm1)
    if the_scf.orbtype == "RHF"
        g = get_jk(the_scf, rdm1)
        return (g+the_scf.hcore)::Array{Float64,2}
    elseif the_scf.orbtype == "UHF"
        error("UHF not implemented yet!")
    end
end

function get_mo_occ(the_scf::SelfConsistentField)
    if the_scf.orbtype == "RHF"
        return get_mo_occ(the_scf.sys)
    elseif the_scf.orbtype == "UHF"
        error("UHF not implemented yet!")
    end
end

function get_rdm1(the_scf::SelfConsistentField, mo_coeff::Array{Float64,2}, mo_occ::Array{Integer,1})
    dm = zeros(Float64, the_scf.nao, the_scf.nao)
    for lm=1:the_scf.nao, sg=1:the_scf.nao
        for p = 1:the_scf.nmo
            dm[lm, sg] = dm[sg, lm] = mo_coeff[p,lm] * mo_coeff[p,sg] * float(mo_occ[p])
        end
    end
    return dm
end

function kernel!(the_scf::SelfConsistentField)
    get_diis_err = (f::Array{Float64,2}, p::Array{Float64,2}) -> norm(reduce(*, [f,p,the_scf.overlp])-reduce(*, [the_scf.overlp,p,f]))

    the_scf.mo_occ = get_mo_occ(the_scf)::Array{Integer,1}
    fock           = get_fock(the_scf, the_scf.init_density_guess)::Array{Float64,2}
    
    iter       = 0
    diis_err   = get_diis_err(fock, the_scf.init_density_guess)

    while !the_scf.is_converged && iter <= the_scf.max_cycle
        eigen_obj = eigen(fock, the_scf.overlp)
        the_scf.mo_energy = eigen_obj.values
        the_scf.mo_coeff  = eigen_obj.vectors

        rdm1              = get_rdm1(the_scf, the_scf.mo_coeff, the_scf.mo_occ)::Array{Float64,2}
        fock              = get_fock(the_scf, rdm1)::Array{Float64,2}
        diis_err          = get_diis_err(fock, rdm1)
        the_scf.is_converged = (diis_err < the_scf.conv_tol)
        println("rdm1 = ",  rdm1)
        iter = iter + 1
    end
end