using Oscar

# Schur polynomial functions
function Schur_polynomial(lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= 0 "n >= 0 required"
    while length(lambda) < n
        push!(lambda, 0)  # Append zeros to the partition if needed
    end
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return Schur_polynomial(R, lambda, n)
end

function Schur_polynomial(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= 0 "n >= 0 required"
    if n == 0 || n < length(lambda)
        if isempty(lambda)
            return one(R)
        else
            return zero(R)
        end
    end
    @req n <= nvars(R) "n <= nvars(R) required"
    return Schur_polynomial_bf(R, lambda, n) #bialternant formula
end

function Schur_polynomial_bf(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= 0 "n >= 0 required"
    while length(lambda) < n
        push!(lambda, 0)  # Append zeros to the partition if needed
    end
    x = gens(R)[1:n]
    A = zero_matrix(R,n,n)
    for i = 1:n
        for j = 1:n
            A[i,j] = x[i]^(lambda[j]+n-j)
        end
    end
    sp = det(A)
    # divide by the product
    for i = 1:n - 1
        for j = i + 1:n
            sp = divexact(sp, x[i] - x[j])
        end
    end
    return sp
end

# Power sum, elementary symmetric, and complete homogeneous symmetric functions
function p(k::Int, n::Int)
    @req n >= 0 "n >= 0 required"
    R, _ = PolynomialRing(ZZ, n, cached = false)
    return p(R, k, n)
end

function p(R::ZZMPolyRing, k::Int, n::Int)
    x = gens(R)[1:n]
    return sum(x[i]^k for i in 1:n)
end

function e(k::Int,n::Int)
    @req n >= 0 "n >= 0 required"
    R, _ = PolynomialRing(ZZ, n, cached = false)
    return e(R,k,n)
end

function e(R::ZZMPolyRing, k::Int, n::Int)
    p = partitions(k)[length(partitions(k))]
    return Schur_polynomial(R, p, n)
end

function h(k::Int, n::Int)
    @req n >= 0 "n >= 0 required"
    R, _ = PolynomialRing(ZZ, n, cached = false)
    return h(R, k, n)
end

function h(R::ZZMPolyRing, k::Int, n::Int)
    p = partition(k)
    return Schur_polynomial(R, p, n)
end

# Grothendieck polynomial functions
function Grothendieck_polynomial(lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= 0 "n >= 0 required"
    while length(lambda) < n
        push!(lambda, 0)  # Append zeros to the partition if needed
    end
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return Grothendieck_polynomial(R, lambda, n)
end

function Grothendieck_polynomial(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= 0 "n >= 0 required"
    if n == 0 || n < length(lambda)
        if isempty(lambda)
            return one(R)
        else
            return zero(R)
        end
    end
    @req n <= nvars(R) "n <= nvars(R) required"
    return Grothendieck_polynomial_bf(R, lambda, n) #bialternant formula
end

function Grothendieck_polynomial_bf(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    @req n > 0 "number of variables must be > 0"
    @req n >= length(lambda) "number of variables must be at least the length of the partition"
    while length(lambda) < n
        push!(lambda, 0)  # Append zeros to the partition if needed
    end
    x = gens(R)[1:n]
    A = zero_matrix(R,n,n)
    for i = 1:n
        for j = 1:n
            A[i,j] = x[i]^(lambda[j]+n-j)*(1-x[i])^(j-1)
        end
    end
    sp = det(A)
    # divide by the product
    for i = 1:n - 1
        for j = i + 1:n
            sp = divexact(sp, x[i] - x[j])
        end
    end
    return sp
end

# Dual Grothendieck polynomial functions
function dual_Grothendieck_polynomial(lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= 0 "n >= 0 required"
    while length(lambda) < n
        push!(lambda, 0)  # Append zeros to the partition if needed
    end
    R, _ = polynomial_ring(ZZ, n, cached = false)
    return dual_Grothendieck_polynomial(R, lambda, n)
end

function dual_Grothendieck_polynomial(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    @req n >= 0 "n >= 0 required"
    if n == 0 || n < length(lambda)
        if isempty(lambda)
            return one(R)
        else
            return zero(R)
        end
    end
    @req n <= nvars(R) "n <= nvars(R) required"
    return dual_Grothendieck_polynomial_bf(R, lambda, n) #bialternant formula
end

function dual_Grothendieck_polynomial_bf(R::ZZMPolyRing, lambda::Vector{Int}, n::Int = length(lambda))
    @req n > 0 "number of variables must be > 0"
    @req n >= length(lambda) "number of variables must be at least the length of the partition"
    while length(lambda) < n
        push!(lambda, 0)  # Append zeros to the partition if needed
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
    # divide by the product
    for i = 1:n - 1
        for j = i + 1:n
            sp = divexact(sp, x[i] - x[j])
        end
    end
    return sp
end

# Modified Schur_expansion function
function Schur_expansion(f)
  R = parent(f)
  n = length(vars(f))  # Number of variables
  result = []  # List to store the results
  while f != zero(f)
      l = collect(exponents(f))
      c = collect(coefficients(f))
      push!(result, (c[1], trim_zeros(l[1])))
      f -= c[1]*Schur_polynomial(R, l[1], n)
  end
  return result
end

# Modified Grothendieck_expansion function
function Grothendieck_expansion(f)
  R = parent(f)
  n = length(vars(f))
  result = []
  while f != zero(f)
      l = collect(exponents(f))
      c = collect(coefficients(f))
      t = reverse(l[end])
      push!(result, (c[end], trim_zeros(t)))
      f -= c[end]*Grothendieck_polynomial(R, t, n)
  end
  return result
end

# Modified dual_Grothendieck_expansion function
function dual_Grothendieck_expansion(f)
  R = parent(f)
  n = length(vars(f))
  result = []
  while f != zero(f)
      l = collect(exponents(f))
      c = collect(coefficients(f))
      push!(result, (c[1], trim_zeros(l[1])))
      f -= c[1]*dual_Grothendieck_polynomial(R, l[1], n)
  end
  return result
end

# Rule functions
function pieri_rule(k::Int, lambda::Vector{Int})
    n = 1 + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = h(R,k,n)*Schur_polynomial(R,lambda,n)
    return Schur_expansion(f)
end

function dual_pieri_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = e(R,k,n)*Schur_polynomial(R,lambda,n)
    return Schur_expansion(f)
end

function murnaghan_nakayama_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = p(R,k,n)*Schur_polynomial(R,lambda,n)
    return Schur_expansion(f)
end

function G_pieri_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = h(R,k,n)*Grothendieck_polynomial(R,lambda,n)
    return Grothendieck_expansion(f)
end

function dual_G_pieri_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = e(R,k,n)*Grothendieck_polynomial(R,lambda,n)
    return Grothendieck_expansion(f)
end

function G_murnaghan_nakayama_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = p(R,k,n)*Grothendieck_polynomial(R,lambda,n)
    return Grothendieck_expansion(f)
end

function g_pieri_rule(k::Int, lambda::Vector{Int})
    n = 1 + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = h(R,k,n)*dual_Grothendieck_polynomial(R,lambda,n)
    return dual_Grothendieck_expansion(f)
end

function dual_g_pieri_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = e(R,k,n)*dual_Grothendieck_polynomial(R,lambda,n)
    return dual_Grothendieck_expansion(f)
end

function g_murnaghan_nakayama_rule(k::Int, lambda::Vector{Int})
    n = k + length(lambda)
    R, _ = polynomial_ring(ZZ, n, cached = false)
    f = p(R,k,n)*dual_Grothendieck_polynomial(R,lambda,n)
    return dual_Grothendieck_expansion(f)
end

# Multiplication functions
# Helper function to trim trailing zeros
function trim_zeros(arr::Vector{Int})
  i = findlast(x -> x != 0, arr)
  return isnothing(i) ? Int[] : arr[1:i]
end

# Modified mult function (already trimmed, but included for completeness)
function mult(lambda::Vector{Int}, mu::Vector{Int}, n::Int=length(lambda)+length(mu))
  R, _ = polynomial_ring(ZZ, n, cached = false)
  f = Schur_polynomial(R,lambda,n)*Schur_polynomial(R,mu,n)
  expansion = Schur_expansion(f)
  return expansion  # Schur_expansion now returns trimmed results
end

# Modified G_mult function
function G_mult(lambda::Vector{Int}, mu::Vector{Int}, n::Int=length(lambda)+length(mu))
  R, _ = polynomial_ring(ZZ, n, cached = false)
  f = Grothendieck_polynomial(R,lambda,n)*Grothendieck_polynomial(R,mu,n)
  return Grothendieck_expansion(f)  # Grothendieck_expansion now returns trimmed results
end

# Modified g_mult function
function g_mult(lambda::Vector{Int}, mu::Vector{Int}, n::Int=length(lambda)+length(mu))
  R, _ = polynomial_ring(ZZ, n, cached = false)
  f = dual_Grothendieck_polynomial(R,lambda,n)*dual_Grothendieck_polynomial(R,mu,n)
  return dual_Grothendieck_expansion(f)  # dual_Grothendieck_expansion now returns trimmed results
end

# Elementary symmetric polynomial
function es(m, j)
    R, x = PolynomialRing(ZZ, m, :x)
    return elementary_symmetric_polynomial(R, x, j)
end
