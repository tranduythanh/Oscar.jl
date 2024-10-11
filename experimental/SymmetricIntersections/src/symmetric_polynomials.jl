using Oscar

# Schur polynomial functions
function schur(lambda::Vector{Int}, n::Int = length(lambda))
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return schur(R, lambda, n)
end

function schur(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    return schur_polynomial(R, partition(lambda), n)
end

# Power sum, elementary symmetric, and complete homogeneous symmetric functions
function p(k::Int, n::Int)
    @req n >= 0 "n >= 0 required"
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return p(R, k, n)
end

function p(R::ZZMPolyRing, k::Int, n::Int)
    x = gens(R)[1:n]
    return sum(x[i]^k for i in 1:n)
end

function e(k::Int, n::Int)
    @req n >= 0 "n >= 0 required"
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return e(R, k, n)
end

function e(R::ZZMPolyRing, k::Int, n::Int)
    return schur_polynomial(R, partitions(k)[length(partitions(k))], n)
end

function h(k::Int, n::Int)
    @req n >= 0 "n >= 0 required"
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return h(R, k, n)
end

function h(R::ZZMPolyRing, k::Int, n::Int)
    return schur_polynomial(R, partition(k), n)
end

# Grothendieck polynomial functions
function grothendieck(lambda::Vector{Int}, n::Int = length(lambda))
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return grothendieck(R, lambda, n)
end

function grothendieck(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    if n == 0 || n < length(lambda)
        if isempty(lambda)
            return one(R)
        else
            return zero(R)
        end
    end
    @req n <= nvars(R) "n <= nvars(R) required"
    return grothendieck_polynomial_bf(R, lambda, n)
end

function grothendieck_polynomial_bf(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= length(lambda) "number of variables must be at least the length of the partition"
    while length(lambda) < n
        push!(lambda, 0)
    end
    x = gens(R)[1:n]
    A = zero_matrix(R, n, n)
    for i = 1:n
        for j = 1:n
            A[i,j] = x[i]^(lambda[j]+n-j)*(1-x[i])^(j-1)
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

# Expansion functions
function schur_expansion(f)
    R = parent(f)
    n = length(vars(f))
    result = []
    while f != zero(f)
        l = filter!(x -> x != 0, collect(exponents(f))[end])
        c = collect(coefficients(f))[end]
        push!(result, (c, reverse(l)))
        f -= c * schur(R, reverse(l), n)
    end
    return result
end

function grothendieck_expansion(f)
    R = parent(f)
    n = length(vars(f))
    result = []
    while f != zero(f)
        l = filter!(x -> x != 0, collect(exponents(f))[end])
        c = collect(coefficients(f))[end]
        push!(result, (c, reverse(l)))
        f -= c * grothendieck(R, reverse(l), n)
    end
    return result
end

function dual_grothendieck_expansion(f)
    R = parent(f)
    n = length(vars(f))
    result = []
    while f != zero(f)
        l = filter!(x -> x != 0, reverse(collect(exponents(f))[1]))
        c = collect(coefficients(f))[1]
        push!(result, (c, reverse(l)))
        f -= c * dual_grothendieck(R, reverse(l), n)
    end
    return result
end

# Rule functions
function pieri_rule(k::Int, lambda::Vector{Int})
    n = 1 + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return schur_expansion(h(R, k, n) * schur(R, lambda, n))
end

function dual_pieri_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return schur_expansion(e(R, k, n) * schur(R, lambda, n))
end

function murnaghan_nakayama_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return schur_expansion(p(R, k, n) * schur(R, lambda, n))
end

function G_pieri_rule(k::Int, lambda::Vector{Int})
    n = 1 + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return grothendieck_expansion(h(R, k, n) * grothendieck(R, lambda, n))
end

function dual_G_pieri_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return grothendieck_expansion(e(R, k, n) * grothendieck(R, lambda, n))
end

function G_murnaghan_nakayama_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return grothendieck_expansion(p(R, k, n) * grothendieck(R, lambda, n))
end

function g_pieri_rule(k::Int, lambda::Vector{Int})
    n = 1 + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return dual_grothendieck_expansion(h(R, k, n) * dual_grothendieck(R, lambda, n))
end

function dual_g_pieri_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return dual_grothendieck_expansion(e(R, k, n) * dual_grothendieck(R, lambda, n))
end

function g_murnaghan_nakayama_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return dual_grothendieck_expansion(p(R, k, n) * dual_grothendieck(R, lambda, n))
end

# Littlewood-Richardson rule functions
function littlewood_richardson_rule(lambda::Vector{Int}, mu::Vector{Int}, n::Int=length(lambda)+length(mu))
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return schur_expansion(schur(R, lambda, n) * schur(R, mu, n))
end

function G_littlewood_richardson_rule(lambda::Vector{Int}, mu::Vector{Int}, n::Int=length(lambda)+length(mu))
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return grothendieck_expansion(grothendieck(R, lambda, n) * grothendieck(R, mu, n))
end

function g_littlewood_richardson_rule(lambda::Vector{Int}, mu::Vector{Int}, n::Int=length(lambda)+length(mu))
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return dual_grothendieck_expansion(dual_grothendieck(R, lambda, n) * dual_grothendieck(R, mu, n))
end
