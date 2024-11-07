function test_qschur_term()
  @testset "QSchurTerm String Formatting" begin
      test_cases = [
          # Basic partitions
          (QSchurTerm(ZZRingElem(1), 0, Partition([1])), "S[1]"),
          (QSchurTerm(ZZRingElem(1), 0, Partition([2,1])), "S[2,1]"),
          (QSchurTerm(ZZRingElem(1), 0, Partition(Int[])), "S[]"),
          
          # With coefficients
          (QSchurTerm(ZZRingElem(2), 0, Partition([1])), "2*S[1]"),
          (QSchurTerm(ZZRingElem(-1), 0, Partition([1])), "-S[1]"),
          (QSchurTerm(ZZRingElem(-2), 0, Partition([1])), "-2*S[1]"),
          
          # With q-powers
          (QSchurTerm(ZZRingElem(1), 1, Partition([1])), "q*S[1]"),
          (QSchurTerm(ZZRingElem(1), 2, Partition([1])), "q^2*S[1]"),
          (QSchurTerm(ZZRingElem(1), 3, Partition([1])), "q^3*S[1]"),
          
          # With both coefficients and q-powers
          (QSchurTerm(ZZRingElem(2), 1, Partition([1])), "2*q*S[1]"),
          (QSchurTerm(ZZRingElem(-1), 1, Partition([1])), "-q*S[1]"),
          (QSchurTerm(ZZRingElem(-2), 2, Partition([1])), "-2*q^2*S[1]"),
          
          # Special cases
          (QSchurTerm(ZZRingElem(0), 0, Partition(Int[])), "0"),
          (QSchurTerm(ZZRingElem(1), 1, Partition(Int[])), "q*S[]"),
          (QSchurTerm(ZZRingElem(-1), 1, Partition(Int[])), "-q*S[]"),
          
          # Complex cases
          (QSchurTerm(ZZRingElem(-2), 3, Partition([2,1])), "-2*q^3*S[2,1]"),
          (QSchurTerm(ZZRingElem(1), 2, Partition([3,2,1])), "q^2*S[3,2,1]"),
          (QSchurTerm(ZZRingElem(-1), 0, Partition([1,1])), "-S[1,1]"),
      ]
      
      for (i, (term, expected)) in enumerate(test_cases)
          @testset "Case $i" begin
              result = string(term)
              @test result == expected
              println("Test $i:")
              println("Input:    ", "coeff=$(term.coeff), q^$(term.q_power), partition=$(term.partition)")
              println("Expected: ", expected)
              println("Got:      ", result)
              println("---")
          end
      end
  end
  
  @testset "Vector{QSchurTerm} String Formatting" begin
      vector_test_cases = [
          # Single term
          ([QSchurTerm(ZZRingElem(1), 0, Partition([1]))], "S[1]"),
          
          # Multiple terms
          ([
              QSchurTerm(ZZRingElem(1), 0, Partition([1])),
              QSchurTerm(ZZRingElem(2), 1, Partition([2]))
          ], "S[1] + 2*q*S[2]"),
          
          # With negative terms
          ([
              QSchurTerm(ZZRingElem(-1), 1, Partition([1])),
              QSchurTerm(ZZRingElem(2), 0, Partition([2]))
          ], "-q*S[1] + 2*S[2]"),
          
          # Empty vector
          (QSchurTerm[], "0"),
          
          # Complex example
          ([
              QSchurTerm(ZZRingElem(-2), 3, Partition([2,1])),
              QSchurTerm(ZZRingElem(1), 0, Partition([1])),
              QSchurTerm(ZZRingElem(-1), 1, Partition(Int[]))
          ], "-2*q^3*S[2,1] + S[1] + -q*S[]")
      ]
      
      for (i, (terms, expected)) in enumerate(vector_test_cases)
          @testset "Vector Case $i" begin
              result = string(terms)
              @test result == expected
              println("Vector Test $i:")
              println("Expected: ", expected)
              println("Got:      ", result)
              println("---")
          end
      end
  end
end
