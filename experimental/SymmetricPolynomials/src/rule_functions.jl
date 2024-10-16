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

function mult(lambda::Partition, mu::Partition, n::Int=length(lambda)+length(mu))
  return littlewood_richardson_rule(collect(lambda), collect(mu), n)
end

# Define a custom type to represent the result
struct QMultTerm
    coeff::ZZRingElem
    q_power::Int
end

function Base.show(io::IO, term::QMultTerm)
    if term.q_power == 0
        print(io, term.coeff)
    else
        print(io, "$(term.coeff) * q^$(term.q_power)")
    end
end

function qmult(
  lambda::Partition, 
  mu::Partition,
  rim_size::Int,
  acceptable_grid::Tuple{Int,Int}
)
  # First, compute the regular Littlewood-Richardson expansion using mult
  expansion = mult(lambda, mu)

  function __rm_rim_hook(part::Vector{Int})
      p = Partition(part)
      if is_in_range(p, acceptable_grid[1], acceptable_grid[2])
          return (p, 0, 0)  # No rim hook removal needed
      end
      return remove_rim_hooks(p, rim_size, acceptable_grid)
  end

  # Apply rim hook removal to each term in the expansion
  new_expansion = Dict{Vector{Int}, QMultTerm}()
  for (coeff, part) in expansion
      new_p, q_num, height = __rm_rim_hook(part)
      
      if !isempty(new_p)
          sign = (-1)^(height - length(lambda))
          new_coeff = coeff * sign
          
          current_term = get(new_expansion, collect(new_p), QMultTerm(ZZRingElem(0), 0))
          new_expansion[collect(new_p)] = QMultTerm(
              current_term.coeff + new_coeff,
              current_term.q_power + q_num
          )
      end
  end

  # Convert the result to the desired format
  result = Vector{Tuple{ZZRingElem, Int, Vector{Int}}}()
  for (part, term) in new_expansion
      push!(result, (term.coeff, term.q_power, part))
  end

  return result
end

function express_as_schur_sum(qmult_result::Vector{Tuple{ZZRingElem, Int, Vector{Int}}})
  terms = String[]
  for (coeff, q_power, part) in qmult_result
      term = ""
      
      # Handle the sign
      if coeff < 0
          term = " - "
          coeff = abs(coeff)
      elseif length(terms) > 0
          term = " + "
      end
      
      # Handle the coefficient
      if coeff != 1
          term *= "$(coeff)*"
      end
      
      # Handle the q factor
      if q_power > 0
          if q_power == 1
              term *= "q*"
          else
              term *= "q^$q_power*"
          end
      end
      
      # Add the Schur function notation
      term *= "S$(part)"
      
      push!(terms, term)
  end
  
  # Join all terms
  result = join(terms, "")
  
  return result
end

function qmult_schur(lambda::Partition, mu::Partition, rim_size::Int, acceptable_grid::Tuple{Int,Int})
  result = qmult(lambda, mu, rim_size, acceptable_grid)
  return express_as_schur_sum(result)
end
