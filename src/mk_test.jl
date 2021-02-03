# mk_test.jl
# module mk_test

export mandelbrot, hello

function mandelbrot(a)
    z = 0
    for i=1:50
        z = z^2 + a
    end
    return z
end

function hello()
    println("hello")
end

# end