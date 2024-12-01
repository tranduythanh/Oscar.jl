function get_k_n(Gr::AbstractVariety)::Tuple{Int,Int}
    k = rank(tautological_bundles(Gr)[1])
    n = dim(Gr)
    return (k, n)
end

function extend_kn(k::Int, n::Int)::Tuple{Int,Int}
    return (k, 2 * n - k)
end

struct QSchurTerm
    coeff::ZZRingElem
    q_power::Int
    partition::Partition
end

# Convert QSchurTerm to string
function Base.string(term::QSchurTerm)::String
    # Use IOBuffer to capture the output of show
    io = IOBuffer()
    show(io, term)
    return String(take!(io))
end

function Base.show(io::IO, term::QSchurTerm)
    if term.coeff == 0
        print(io, "0")
        return nothing
    end

    if term.coeff != 1
        if term.coeff == -1
            print(io, "-")
        else
            print(io, term.coeff)
        end

        if term.q_power > 0 || !isempty(term.partition.p)
            if term.coeff != 1 && term.coeff != -1
                print(io, "*")
            end
        end
    end

    if term.q_power > 0
        print(io, "q")
        if term.q_power > 1
            print(io, "^", term.q_power)
        end
        print(io, "*")
    end

    if isempty(term.partition.p)
        if term.coeff != 0
            print(io, "S[]")
        end
    else
        print(io, "S[", join(term.partition.p, ","), "]")
    end
end

function Base.string(terms::Vector{QSchurTerm})::String
    if isempty(terms)
        return "0"
    end
    if length(terms) == 1
        return string(terms[1])  # Use the existing string method for single QSchurTerm
    end
    return join(string.(terms), " + ")  # For multiple terms
end

function qmult(
    Gr::AbstractVariety,
    lambda::Partition,
    mu::Partition,
)::Vector{QSchurTerm}

    if length(mu) == 0
        return [QSchurTerm(ZZRingElem(1), 0, lambda)]
    end

    expansion = mult(lambda, mu)

    (k, n) = get_k_n(Gr)
    (k, n) = extend_kn(k, n)

    rim_size = n
    acceptable_grid = (k, n - k)

    function __rm_rim_hook(part::Vector{Int})
        p = Partition(part)
        if is_in_range(p, acceptable_grid[1], acceptable_grid[2])
            return (p, 0, 0)
        end
        return remove_rim_hooks(p, rim_size, acceptable_grid)
    end

    new_expansion = Dict{Vector{Int},QSchurTerm}()
    for (coeff, part) in expansion
        new_p, q_num, height = __rm_rim_hook(part)

        if !isempty(new_p)
            sign = (-1)^(height - length(lambda))
            new_coeff = coeff * sign

            part_vec = collect(new_p)

            current_term = get(
                new_expansion, part_vec, QSchurTerm(ZZRingElem(0), 0, Partition(Int[]))
            )

            new_expansion[part_vec] = QSchurTerm(
                current_term.coeff + new_coeff,
                current_term.q_power + q_num,
                Partition(part_vec),
            )
        end
    end

    println(lambda)
    println(mu)
    println(new_expansion)

    return collect(values(new_expansion))
end

function qmult_schur(Gr::AbstractVariety, lambda::Partition{Int}, mu::Partition{Int})
    result = qmult(Gr, lambda, mu)
    return express_as_schur_sum(result)
end
