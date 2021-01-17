abstract type AbstractTensor end
#
struct TwoElectronIntegralAO{T<:Number} <: AbstractTensor
    _index::Array{Integer, 2}
    _data ::Array{T}
end

function get_eri_index(lm::Integer, sgm::Integer, mu::Integer, nu::Integer)::Vector{Integer}
    if lm >= sgm
        a = lm
        b = sgm
    else
        a = sgm
        b = lm
    end

    if mu >= nu
        c = mu
        d = nu
    else
        c = nu
        d = mu
    end

    if a > c || (a == c && b > d)
        return Vector{Integer}([a, b, c, d])
    else
        return Vector{Integer}([c, d, a, b])
    end
end

function get_value(two_elec_int::TwoElectronIntegralAO{T}, lm::Integer, sgm::Integer, mu::Integer, nu::Integer)::T where {T}
    eri_index = get_eri_index(lm, sgm, mu, nu)::Vector{Integer}
    tmp_val   = zero(T)

    is_index_in_two_elec_int = false
    for i in 1:size(two_elec_int._index,1)
        @inbounds is_index_in_two_elec_int = (two_elec_int._index[i,:] == eri_index)
        if is_index_in_two_elec_int
            @inbounds tmp_val = two_elec_int._data[i]::T
            break
        end
    end

    return tmp_val
end

function read_data(int_path::AbstractString; FloatType=Float64, IntType=Int64)
# Reads integrals and places them in array format for future use in 
# quantum chemical programs
    RealType      = real(FloatType)
    e_nuc_path    = join(["./", int_path, "/nuclear_repulsion.dat"])
    num_elec_path = join(["./", int_path, "/number_electrons.dat"])
    nbas_path     = join(["./", int_path, "/number_basis_functions.dat"])
    
    e_nuc::RealType             = readdlm(e_nuc_path,    ',', RealType)[1]
    num_elec::Array{IntType, 2} = readdlm(num_elec_path, ',', IntType)
    nbas::IntType               = readdlm(nbas_path,     ',', IntType)[1]

    return e_nuc::RealType, (num_elec[1], num_elec[2]), nbas
end

function read_ovlp(int_path::AbstractString; FloatType=Float64, IntType=Int64)
    RealType = real(FloatType)
    s_int_path                   = join(["./", int_path, "/overlap.dat"])
    s_matrix::Array{FloatType,2} = readdlm(s_int_path, ',', RealType)
    return Hermitian(s_matrix)
end

function read_hcore(int_path::AbstractString; FloatType=Float64, IntType=Int64)
    RealType = real(FloatType)
    t_int_path    = join(["./", int_path, "/kinetic_energy.dat"])
    v_int_path    = join(["./", int_path, "/potential_energy.dat"])
    t_matrix::Array{FloatType,2} = readdlm(t_int_path, ',', RealType)
    v_matrix::Array{FloatType,2} = readdlm(v_int_path, ',', RealType)
    h_matrix = t_matrix-v_matrix
    return Hermitian(h_matrix)
end

function read_eri(int_path::AbstractString; FloatType=Float64, IntType=Int64)
    RealType = real(FloatType)
    eri_index_raw_path  = join(["./", int_path, "/eri_index.dat"])
    eri_val_raw_path    = join(["./", int_path, "/eri_val.dat"])
    eri_index::Array{IntType,2} = readdlm(eri_index_raw_path, ',', IntType)
    eri_val::Array{FloatType}   = readdlm(eri_val_raw_path, ',', RealType)
    return TwoElectronIntegralAO{FloatType}(eri_index, eri_val)
end