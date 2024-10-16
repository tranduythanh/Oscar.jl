module SymmetricPolynomials

using Oscar

# Export the main functions
export schur, grothendieck, dual_grothendieck
export p, e, h
export schur_expansion, grothendieck_expansion, dual_grothendieck_expansion
export pieri_rule, dual_pieri_rule, murnaghan_nakayama_rule
export G_pieri_rule, dual_G_pieri_rule, G_murnaghan_nakayama_rule
export g_pieri_rule, dual_g_pieri_rule, g_murnaghan_nakayama_rule
export littlewood_richardson_rule, G_littlewood_richardson_rule, g_littlewood_richardson_rule
export Partition, draw, is_in_range, remove_rim_hooks, skew_partition_height, is_non_increasing, mult, qmult, qmult_schur

include("schur.jl")
include("symmetric_functions.jl")
include("grothendieck_polynomials.jl")
include("dual_grothendieck_polynomials.jl")
include("expansion_functions.jl")
include("rule_functions.jl")
include("remove_rim_hooks.jl")

end
