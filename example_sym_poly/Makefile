.PHONY: precompile run

JULIA_VERSION := 1.10.4

precompile:
	julia precompile.jl

run:
	julia -J Oscar_$(JULIA_VERSION).so main.jl
