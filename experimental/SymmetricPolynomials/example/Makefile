JULIA_VERSION := 1.10.6

precompile:
	julia precompile.jl

run:
	julia -i --sysimage Oscar.so main.jl

test:
	julia --sysimage Oscar.so ../test/test.jl
