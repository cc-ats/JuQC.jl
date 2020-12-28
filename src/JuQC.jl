"""
Julia Quantum Chemistry Package
"""
module JuQC

using LinearAlgebra, DelimitedFiles
using Printf

include("./read_integrals.jl")
export read_data, read_ovlp, read_hcore, read_eri
export get_value

end # module