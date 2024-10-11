function pieri_rule(k::Int, lambda::Vector{Int})
  n = 1 + length(lambda)
  R, _ = polynomial_ring(ZZ, n, cached = false)
  return schur_expansion(h(R, k, n) * schur(R, lambda, n))
end

function pieri_rule(k::Int, lambda::Partition)
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
