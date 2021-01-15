"""
Julia Quantum Chemistry Package
"""
module JuQC

using LinearAlgebra, DelimitedFiles
using Printf

include("./utils.jl")
include("./read_integrals.jl")
include("./scf.jl")
include("./scf_algo.jl")


export SCFSolver, RestrictedSCFSolver, UnrestrictedSCFSolver
export SCFAlgorithm, Roothaan
export read_data, read_ovlp, read_hcore, read_eri
export get_value, build_scf_solver, init_scf_algo!, kernel!

end # module