@testset "Measure" begin
    @polyvar x y
    @test_throws ErrorException Measure([1, 2], [x, x*y, y])
    @test_throws ErrorException Measure([1, 2, 3, 4], MonomialVector([x, x*y, y]))
    m = Measure([1, 0, 2, 3], [x^2*y^2, y*x^2, x*y*x^2, x*y^2])
    @test m.a == [2, 1, 0, 3]
end

# [HL05] Henrion, D. & Lasserre, J-B.
# Detecting Global Optimality and Extracting Solutions of GloptiPoly 2005

@testset "[HL05] Section 2.3" begin
    @polyvar x y
    ν = AtomicMeasure([x, y], [0.4132, 0.3391, 0.2477], [[1, 2], [2, 2], [2, 3]])
    μ = Measure(ν, [x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1])
    μ = MatMeasure(μ, [1, x, y, x^2, x*y, y^2])
    atoms = extractatoms(μ, 1e-4, 1e-14)
    @test isapprox(atoms, ν)
end

#   @testset "[HL05] Section 4" begin
#       @polyvar x y
#       μ = Measure([1/9,     0,     1/9,     0, 1/9,   0,     0,     0,   0, 1/3,   0, 1/3, 0, 0, 1],
#                   [x^4, x^3*y, x^2*y^2, x*y^3, y^4, x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1])
#       μ = MatMeasure(μ, [1, x, y, x^2, x*y, y^2])
#       @show extractatoms(μ, 1e-16, 1e-16)
#   end
