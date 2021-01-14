"""
Julia Quantum Chemistry Package
"""
module JuQC

using LinearAlgebra, DelimitedFiles
using Printf

include("./utils.jl")
include("./read_integrals.jl")
include("./opt.jl")
include("./scf.jl")

export read_data, read_ovlp, read_hcore, read_eri
export get_value, build_scf_solver, kernel

end # module