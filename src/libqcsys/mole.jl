struct Mole <: QuantumChemistrySystem
    overlp::Array{Float64, 2}
    hcore::Array{Float64, 2}
    eri::Array{Float64, 4}
end

function get_nuclear_energy(mol::Mole)
end

function get_mo_occ(mol::Mole)
end

function get_dimension(mol::Mole)
end

function get_hcore_matrix(mol::Mole)
end

function get_overlap_matrix(mol::Mole)
end

function get_eri_array(mol::Mole)
end
