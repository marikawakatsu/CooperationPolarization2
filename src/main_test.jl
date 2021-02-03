# main.jl
include("mk_test.jl")

# using .mk_test

for y=1.0:-0.05:-1.0
    for x=-1.75:0.0315:0.5
        abs(mandelbrot(complex(x, y))) < 2 ? print("*") : print(" ")
    end
    println()
end

if "test" != "mari"
    @warn "test and mari are different"
end

hello()

# Taken from: https://rosettacode.org/wiki/Mandelbrot_set#Julia