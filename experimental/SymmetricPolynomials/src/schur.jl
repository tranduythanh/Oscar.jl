"""
    schur(lambda::Vector{Int}, n::Int = length(lambda))

Compute the Schur polynomial for the partition given by `lambda` with `n` variables.
"""
function schur(lambda::Vector{Int}, n::Int = length(lambda))
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return schur(R, lambda, n)
end

"""
    schur(lambda::Partition, n::Int = length(lambda))

Compute the Schur polynomial for the given partition `lambda` with `n` variables.
"""
function schur(lambda::Partition, n::Int = length(lambda))
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return schur(R, lambda, n)
end

"""
    schur(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))

Compute the Schur polynomial in the ring `R` for the partition given by `lambda` with `n` variables.
"""
function schur(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    return schur_polynomial(R, partition(lambda), n)
end

"""
    schur(R::ZZMPolyRing, lambda::Partition, n::Int = length(lambda))

Compute the Schur polynomial in the ring `R` for the given partition `lambda` with `n` variables.
"""
function schur(R::ZZMPolyRing, lambda::Partition, n::Int = length(lambda))
    return schur_polynomial(R, lambda, n)
end
