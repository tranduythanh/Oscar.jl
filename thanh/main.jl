# Import the library from the specified path
include("../experimental/SymmetricIntersections/src/symmetric_polynomials.jl")

# Now you can use any function from the library
# For example, let's calculate a Schur polynomial

# Define a partition
lambda = [2, 1]

# Calculate the Schur polynomial for this partition with 3 variables
result = Schur_polynomial(lambda, 3)

# Print the result
println("Schur polynomial for partition [2, 1] with 3 variables:")
println(result)

# Let's try another function, like multiplying two Schur polynomials
lambda1 = [2, 1]
lambda2 = [1, 1]
mult_result = mult(lambda1, lambda2)

println("\nMultiplication of Schur polynomials [2, 1] and [1, 1]:")
for (coeff, partition) in mult_result
    println("$coeff * s$partition")
end
