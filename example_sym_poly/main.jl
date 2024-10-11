using Oscar

# path to SymmetricPolynomials module
push!(LOAD_PATH, joinpath(@__DIR__, "..", "experimental", "SymmetricPolynomials", "src"))

# include the module file
include(joinpath(@__DIR__, "..", "experimental", "SymmetricPolynomials", "src", "SymmetricPolynomials.jl"))

# Import the module
using .SymmetricPolynomials

# Define a partition
lambda = [2, 1]

# Calculate the Schur polynomial for this partition with 3 variables
result = schur(lambda, 3)

# Print the result
println("Schur polynomial for partition [2, 1] with 3 variables:")
println(result)

# Let's try another function, like multiplying two Schur polynomials using the Littlewood-Richardson rule
lambda1 = [2, 1]
lambda2 = [1, 1]
mult_result = littlewood_richardson_rule(lambda1, lambda2)

println("\nMultiplication of Schur polynomials [2, 1] and [1, 1]:")
for (coeff, part) in mult_result
    println("$coeff * s$part")
end

# Let's also demonstrate a Grothendieck polynomial calculation
g_result = grothendieck(lambda, 3)

println("\nGrothendieck polynomial for partition [2, 1] with 3 variables:")
println(g_result)

# And a dual Grothendieck polynomial calculation
dg_result = dual_grothendieck(lambda, 3)

println("\nDual Grothendieck polynomial for partition [2, 1] with 3 variables:")
println(dg_result)

# Demonstrate the use of symmetric functions
println("\np(2, 3) = ", p(2, 3))
println("e(2, 3) = ", e(2, 3))
println("h(2, 3) = ", h(2, 3))

println(dual_pieri_rule(3,[2,1]))
println(dual_G_pieri_rule(3,[2,1]))
println(dual_g_pieri_rule(3,[2,1]))
println(g_murnaghan_nakayama_rule(3,[2,1]))
println(g_littlewood_richardson_rule([2,1],[2,1]))
println(G_littlewood_richardson_rule([2,1],[2,1]))

println(littlewood_richardson_rule([2,1],[2,1]))
println(murnaghan_nakayama_rule(3,[2,1]))
println(G_murnaghan_nakayama_rule(3,[2,1]))
