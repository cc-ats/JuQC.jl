struct HubbardModel <: QuantumChemistrySystem
    dim::Integer
    para_t::Float64
    para_u::Float64
    is_pbc::Bool
    function HubbardModel(dim::Integer, para_t::Float64, para_u::Float64, is_pbc::Bool=true)
        println("Initializing a Hubbard moldel.")
        @printf("dim    = % 4d , is_pbc = % 4s\n", dim, is_pbc)
        @printf("para_t = % 4.2f, para_u = % 4.2f\n", para_t, para_u)
        println("")
        new(dim, para_t, para_u, is_pbc)
    end
end

function get_nuclear_energy(hm::HubbardModel)
    return 0.0
end

function get_mo_occ(hm::HubbardModel; orbtype::String="RHF")
    if orbtype == "RHF"
        occ = ones(Integer, hm.dim)
        return occ
    elseif orbtype == "UHF"
        error("UHF not implemented yet!")
    end
end

function get_dimension(hm::HubbardModel)
    return (hm.dim, hm.dim)
end

function get_hcore_matrix(hm::HubbardModel)
    h1 = zeros(Float64, hm.dim, hm.dim)
    for i = 1:(hm.dim-1)
        h1[i,i+1] = h1[i+1,i] = -hm.para_t
    end
    if hm.is_pbc
        h1[hm.dim,1]   = h1[hm.dim,1]   = -hm.para_t
    end
    return h1
end

function get_overlap_matrix(hm::HubbardModel)
    return Matrix{Float64}(1.0*I, hm.dim, hm.dim)
end

function get_eri_array(hm::HubbardModel)
    hubbard_eri = zeros(Float64, hm.dim, hm.dim, hm.dim, hm.dim)
    for i in 1:hm.dim
        hubbard_eri[i,i,i,i] = hm.para_u
    end
    return hubbard_eri
end