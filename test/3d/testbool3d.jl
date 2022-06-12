using LARgenerators, SparseArrays, Base.Threads
import LinearAlgebraicRepresentation as Lar
import ViewerGL as GL
using Test

n, m, p = 1, 1, 1
V, (VV, EV, FV, CV) = Lar.cuboidGrid([n, m, p], true)
cube = V, FV, EV

assembly = Lar.Struct([
    Lar.Struct([Lar.t(0, 0, 0), Lar.r(0, 0, 0), cube])
    Lar.Struct([Lar.t(0.2, 0, 0.2), Lar.r(0.1, 0, 0), cube])
    Lar.Struct([Lar.t(0, 0.1, 0.2), Lar.r(0, 0.2, 0), cube])
])

V, FV, EV = Lar.struct2lar(assembly)
cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells))
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells)
W = convert(Lar.Points, V')

Z, copEV, copFE, copCF = Lar.space_arrangement(W, cop_EV, cop_FE)
W = convert(Lar.Points, Z')

innerpoints = LARgenerators.internalpoints(W, copEV, copFE, copCF[2:end, :])

listOfModels = Lar.evalStruct(assembly)
inputfacenumbers = [length(listOfModels[k][2]) for k = 1:length(listOfModels)]
cumulative = cumsum([0; inputfacenumbers]) .+ 1
fspans = collect(zip(cumulative[1:end-1], cumulative[2:end] .- 1))
span(h) = [j for j = 1:length(fspans) if fspans[j][1] <= h <= fspans[j][2]]

W, copEV, copFE, copCF, boolmatrix = LARgenerators.bool3d(assembly)

@testset "@code_warntype for bool3d" begin
    err = nothing
    try
        @inferred LARgenerators.internalpoints(W, copEV, copFE, copCF[2:end, :])
    catch err
    end
    @test isa(err, ErrorException) == false

    err = nothing
    try
        @inferred LARgenerators.testinternalpoint(V, EV, FV, x)(innerpoints[1])
    catch err
    end
    @test isa(err, ErrorException) == false

    err = nothing
    try
        @inferred LARgenerators.bool3d(assembly)
    catch err
    end
    @test isa(err, ErrorException) == false

end

@testset "Generic result of bool3d" begin
    @test boolmatrix.dims[1] > 0 && boolmatrix.dims[2] > 0 # dimensione array > ( 0, 0 )
    @test length(filter(x -> x > 0, boolmatrix)) > 0
    @test boolmatrix.dims[1] * boolmatrix.dims[2] >= length(filter(x -> x > 0, boolmatrix)) > 0
end

@testset "Specific result of bool3d" begin
    @test boolmatrix.dims[1] == 8 && boolmatrix.dims[2] == 4
end