function test_qmult()
  @testset "Quantum Multiplication" begin
      @testset "Test case 1" begin
          Gr = abstract_grassmannian(1, 2)
          res = qmult(Gr, Partition([1]), Partition(Int[]))
          println(res)
          @test res[1].partition.p == [1]
          @test res[1].q_power == 0
          @test res[1].coeff == 1
      end
  end
end


# Helper functions to parse test cases
function parse_variety(str::AbstractString)
  m = match(r"Gr\((\d+),(\d+)\)", str)
  if m === nothing
      error("Invalid Grassmannian format: $str")
  end
  k, n = parse.(Int, m.captures)
  return abstract_grassmannian(k, n)
end

function parse_partition(str::AbstractString)
  if str == "S[]"
      return Partition(Int[])
  end
  
  m = match(r"S\[([\d,]+)\]", str)
  if m === nothing
      error("Invalid partition format: $str")
  end
  
  nums = split(m.captures[1], ",")
  return Partition([parse(Int, n) for n in nums if !isempty(n)])
end

function parse_expected_result(str::AbstractString)
  # Handle empty result
  if str == "S[]"
      return (ZZRingElem(1), 0, Partition(Int[]))
  end
  
  # Parse q coefficient
  q_power = 0
  coeff = ZZRingElem(1)
  
  if startswith(str, "q")
      q_power = 1  # For Grassmannian cases, q power is always 1
      str = replace(str, "q*" => "")
  end
  
  # Parse partition
  partition = parse_partition(str)
  
  return (coeff, q_power, partition)
end

function test_qmult()
  @testset "Grassmannian Quantum Multiplication" begin
      test_cases = readlines(joinpath(@__DIR__, "..", "test", "testcases.csv"))
      
      for (i, line) in enumerate(test_cases)
          parts = split(line, ";")
          Gr = parse_variety(parts[1])
          k, n = get_k_n(Gr)
          
          @testset "Gr($k,$n)" begin
              lambda = parse_partition(parts[3])
              mu = parse_partition(parts[4])
              
              # Run test
              result = qmult(Gr, lambda, mu)
              result_str = string(result)

              # Debug information
              println("Testcase: $line")
              println("Got:      $result_str")
              
              # Verify results
              @test result_str == parts[5]
          end
      end
  end
end
