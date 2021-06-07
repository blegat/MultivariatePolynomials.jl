@testset "Differentiation" begin
    Mod.@polyvar x y z
    @test differentiate(3, y) == 0
    @test differentiate.([x, y], y) == [0, 1]
    @test differentiate([x, y], (x, y)) == Matrix(1.0I, 2, 2) # TODO: this can be just `I` on v0.7 and above
    @test differentiate(true*x+true*x^2, y) == 0
    @inferred differentiate(true*x+true*x^2, y)
    @test differentiate(x*y + 3y^2 , [x, y]) == [y, x+6y]
    @inferred differentiate(x*y + 3y^2 , x)
    @inferred differentiate(x*y + 3y^2 , [x, y])
    @test differentiate(1 / x , [x, y]) == [-1/x^2, 0]
    @test differentiate((x - y) / (x * y) , [x, y]) == [y^2 / (x * y)^2, -x^2 / (x * y)^2]
    p = x^2-y+x*y
    @test differentiate(p, x, 2) == 2
    @test differentiate(p, x, 0) === p
    @test_throws DomainError differentiate(p, x, -1)
    @test differentiate(x, x, 1) == 1
    #@inferred differentiate(x, x, 0) # FIXME failing at the moment
    @inferred differentiate(x, x, 1)
    @test differentiate(4x*y^3, y) == 12x*y^2
    @inferred differentiate(4x*y^3, y)
    @test differentiate(differentiate(x*y^4, y), y) == 12x*y^2
    @inferred differentiate(differentiate(x*y^4, y), y)
    @test differentiate(x*y^4, y, 3) == 24x*y
    @inferred differentiate(x*y^4, y, 3)
    @test differentiate(x*y^4, y, 0) == x*y^4
    @test differentiate(2x^2, x, 2) == 4
    @inferred differentiate(2x^2, x, 2)
    @inferred differentiate(2x^2 + x*y + y^2, [x, y])
    p = differentiate(2x^2 + 3x*y + y^2, [x, y], 2)
    @test isa(p, Matrix{<:AbstractPolynomial{Int}})
    @test p == [4 3; 3 2]
    p = differentiate(2x^2 + 3x*y^2 + 4y^3 + 2.0, (x, y), 2)
    @test isa(p, Matrix{<:AbstractPolynomial{Float64}})
    @test p == [4.0 6.0y; 6.0y 6.0x+24.0y]

    f = [x^2+y, z^2+4x]
    @test differentiate(f, [x, y, z]) == [2x 1 0; 4 0 2z]
    @test differentiate(f, (x, y, z)) == [2x 1 0; 4 0 2z]

    @testset "Differentiate empty polynomial" begin
        p = x^0 - 1
        @test iszero(@inferred differentiate(p, x))
        @test all(iszero, @inferred differentiate(p, [x, y]))
    end

    @testset "differentiation with Val{}" begin
        @test @inferred(differentiate(x, x, Val{0}())) == x
        @test @inferred(differentiate(x, x, Val{1}())) == 1
        @test @inferred(differentiate(x^2, x, Val{1}())) == 2x
        @test @inferred(differentiate(x^2, y, Val{1}())) == 0
        @test @inferred(differentiate(2x^2 + 3y, x, Val{1}())) == 4x
        @test @inferred(differentiate(2x^2 + 3y, x, Val{2}())) == 4
        p = differentiate(2x^2 + 3x*y + y^2, [x, y], Val{2}())
        @test isa(p, Matrix{<:AbstractPolynomial{Int}})
        @test p == [4 3; 3 2]
    end
end
