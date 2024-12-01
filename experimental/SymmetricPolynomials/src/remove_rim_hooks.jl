# Extend Base.show for our specific use case if not already defined
if !hasmethod(Base.show, Tuple{IO, Partition})
  Base.show(io::IO, p::Partition) = print(io, collect(p))
end

function part_clip(lambda::Vector{Int})::Vector{Int}
  # Trims or removes trailing zeros from the vector lambda.
  i = findlast(!iszero, lambda)
  isnothing(i) ? Int[] : lambda[1:i]
end

function is_non_increasing(part::Vector{Int})::Bool
  # Check if the vector is in non-increasing order
  isempty(part) && return true
  all(i -> part[i] >= part[i+1], 1:length(part)-1)
end

function draw(p::Partition)
  isempty(p) && return println("0")
  foreach(row -> println("[]"^row), p)
end

function is_in_range(p::Partition, nrow::Int, ncol::Int)::Bool
  (nrow >= 0 && ncol >= 0) || throw(ArgumentError("nrow and ncol must be non-negative"))
  isempty(p) && return true
  length(p) <= nrow || return false
  p[1] <= ncol || return false
  return true
end

function remove_rim_hooks(p::Partition, rim_size::Int, acceptable_grid::Tuple{Int,Int})::Tuple{Partition,Int,Int}
  isempty(p) && return Partition(Int[]), 0, 0
  rim_size <= 0 && return Partition(collect(p)), 0, 0
  rim_size > sum(p) + length(p) - 1 && return Partition(Int[]), 0, 0

  nrow, ncol = acceptable_grid
  current_partition = collect(p)
  total_rim_hooks_removed = 0
  
  while true
      _partition, rim_hook_removed = _remove_single_rim_hook(current_partition, rim_size)
      
      if rim_hook_removed
          if !is_non_increasing(_partition)
              return Partition(Int[]), total_rim_hooks_removed, skew_partition_height(Partition(Int[]), p)
          end

          new_partition = Partition(_partition)
          total_rim_hooks_removed += 1

          if is_in_range(new_partition, nrow, ncol)
              return new_partition, total_rim_hooks_removed, skew_partition_height(new_partition, p)
          end

          if collect(new_partition) == current_partition
              return Partition(Int[]), total_rim_hooks_removed, skew_partition_height(Partition(Int[]), p)
          end

          current_partition = collect(new_partition)
      else
          return Partition(Int[]), 0, 0
      end
  end
end

function skew_partition_height(lambda::Partition, mu::Partition)::Int
  isempty(lambda) && return length(mu)

  _validate_partitions(lambda, mu)

  count(i -> mu[i] > (i <= length(lambda) ? lambda[i] : 0), 1:length(mu))
end

# Helper functions

function _remove_single_rim_hook(partition::Vector{Int}, rim_size::Int)::Tuple{Vector{Int},Bool}
  _partition = vcat(partition, 0)
  _rim_size = rim_size
  
  for i in 1:length(_partition)-1
      delta = _partition[i] - _partition[i + 1] + (_partition[i + 1] > 0)
      if delta >= _rim_size
          _partition[i] -= _rim_size
          return part_clip(_partition), true
      end
      _rim_size -= delta
      _partition[i] -= delta
      if _rim_size <= 0 && _partition[i] >= _partition[i+1]
          return part_clip(_partition), true
      end
  end
  
  return partition, false
end

function _validate_partitions(lambda::Partition, mu::Partition)
  for (i, lambda_i) in enumerate(lambda)
      i > length(mu) && throw(ArgumentError("Partition lambda is not contained within partition mu (lambda has more rows)."))
      lambda_i > mu[i] && throw(ArgumentError("Partition lambda is not contained within partition mu (lambda[$i] > mu[$i])."))
  end
end
