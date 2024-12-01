using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using Oscar:partition,ZZRingElem,Partition

# Include the SymmetricPolynomials module
include(joinpath(@__DIR__, "..", "src", "SymmetricPolynomials.jl"))
using .SymmetricPolynomials

include("test_qmult.jl")
include("test_remove_rim_hooks.jl")
include("test_qschur_term.jl")

# test_basic_rim_hook()
test_qschur_term()
test_qmult()
