using Test, SparseArrays, LARgenerators, Base.Threads, SparseArrays, IntervalTrees, LinearAlgebra
import LinearAlgebraicRepresentation as Lar

V, (VV, EV, FV) = Lar.cuboidGrid([1, 1], true)
square = V, EV

assembly = Lar.Struct([
    Lar.Struct([Lar.t(0, 0), Lar.r(0), square])
    Lar.Struct([Lar.t(0, 0.1), Lar.r(0.1), square])
    Lar.Struct([Lar.t(0, 0.2), Lar.r(0.2), square])
    Lar.Struct([Lar.t(0, 0.3), Lar.r(0.3), square])
    Lar.Struct([Lar.t(0, 0.4), Lar.r(0.4), square])
])

V, EV = Lar.struct2lar(assembly)
cop_EW = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells))
Z = convert(Lar.Points, V')
W, copEV, copFE = Lar.planar_arrangement(Z::Lar.Points, cop_EW::Lar.ChainOp)
innerpoints = LARgenerators.internalpoints2d(W, copEV, copFE)
listOfModels = Lar.evalStruct(assembly)
boolmatrix = BitArray(undef, copFE.m, length(listOfModels))
containmenttest = LARgenerators.testinternalpoint2d(listOfModels)(innerpoints[1])
W, copEV, copFE, boolmatrix = LARgenerators.bool2d(assembly)

@testset "@code_warntype for bool2d" begin
    err = nothing
    try
        @inferred LARgenerators.internalpoints2d(W, copEV, copFE)
    catch err
    end
    @test isa(err, ErrorException) == false

    err = nothing
    try
        @inferred LARgenerators.testinternalpoint2d(listOfModels)(innerpoints[1])
    catch err
    end
    @test isa(err, ErrorException) == false

    err = nothing
    try
        @inferred LARgenerators.bool2d(assembly)
    catch err
    end
    @test isa(err, ErrorException) == false

end

@testset "Generic result of bool2d" begin
    @test boolmatrix.dims[1] > 0 && boolmatrix.dims[2] > 0 # dimensione array > ( 0, 0 )
    @test length(filter(x -> x > 0, boolmatrix)) > 0
    @test boolmatrix.dims[1] * boolmatrix.dims[2] >= length(filter(x -> x > 0, boolmatrix)) > 0
end

@testset "Specific result of bool2d" begin
    @test boolmatrix.dims[1] == 30 && boolmatrix.dims[2] == 5
end