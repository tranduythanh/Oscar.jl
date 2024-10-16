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
  parts = collect(partitions(k))
  return schur(R, parts[end], n)
end

function h(k::Int, n::Int)
  @req n >= 0 "n >= 0 required"
  R, _ = polynomial_ring(ZZ, n, cached = false)
  return h(R, k, n)
end

function h(R::ZZMPolyRing, k::Int, n::Int)
  return schur(R, [k], n)
end
