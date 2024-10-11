# Dual Grothendieck polynomial functions
function dual_grothendieck(lambda::Vector{Int}, n::Int = length(lambda))
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return dual_grothendieck(R, lambda, n)
end

function dual_grothendieck(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    if n == 0 || n < length(lambda)
        if isempty(lambda)
            return one(R)
        else
            return zero(R)
        end
    end
    @req n <= nvars(R) "n <= nvars(R) required"
    return dual_grothendieck_polynomial_bf(R, lambda, n)
end

function dual_grothendieck_polynomial_bf(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= length(lambda) "number of variables must be at least the length of the partition"
    while length(lambda) < n
        push!(lambda, 0)
    end
    x = gens(R)[1:n]
    A = zero_matrix(R,n,n)
    for j = 1:n
        A[1,j] = x[j]^(lambda[1]+n-1)
        for i = 2:n
            A[i,j] = sum(binomial(i+k-2,k)*x[j]^(lambda[i] + n - i - k) for k = 0:(lambda[i] + n - i))
        end
    end
    sp = det(A)
    for i = 1:n - 1
        for j = i + 1:n
            sp = divexact(sp, x[i] - x[j])
        end
    end
    return sp
end
