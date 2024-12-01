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
