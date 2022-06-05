using SparseArrays, LARgenerators, Base.Threads, SparseArrays, IntervalTrees, LinearAlgebra
import ViewerGL as GL
import LinearAlgebraicRepresentation as Lar
using Statistics
using BenchmarkTools

function chainbasis2polygonsEx(V, copEV, copFE)
    FE = [findnz(copFE[k, :])[1] for k = 1:copFE.m]
    EV = [findnz(copEV[k, :])[1] for k = 1:copEV.m]

    FEs = []
    EVs = []
    FVs = []

    for f = 1:copFE.m
        push!(FEs, collect(Set(cat([e for e in FE[f]]; dims=1))))
        push!(EVs, [EV[e] for e in FE[f]])
        push!(FVs, collect(Set(cat([EV[e] for e in FE[f]]; dims=1))))
    end
    polygons = collect(zip(EVs, FVs, FEs))
    W = convert(Lar.Points, V')
    return W, polygons, FE
end

function internalpoints2dEx(W, copEV, copFE)
    U, pols, FE = chainbasis2polygonsEx(W, copEV, copFE)
    internalpoints = []
    for f = 1:length(pols)
        (EV, FV, FE) = pols[f]
        internalpoint = Lar.getinternalpoint2d(W, EV, FV, f, copEV, copFE)
        push!(internalpoints, internalpoint)
    end
    return internalpoints
end

function bool2dEx(assembly)

    V, EV = Lar.struct2lar(assembly)
    cop_EW = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells))
    W = convert(Lar.Points, V')

    W, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

    innerpoints = internalpoints2dEx(W, copEV, copFE[1:end, :])

    listOfModels = Lar.evalStruct(assembly)
    inputfacenumbers = [length(listOfModels[k][2]) for k = 1:length(listOfModels)]

    boolmatrix = BitArray(undef, length(innerpoints), length(listOfModels))
    containmenttest = Lar.testinternalpoint2d(listOfModels)
    for (k, point) in enumerate(innerpoints)
        cells = containmenttest(point)
        for l in cells
            boolmatrix[k, l] = 1
        end
    end
    return W, copEV, copFE, boolmatrix
end

T = 50

V, (VV, EV, FV) = Lar.cuboidGrid([1, 1], true)
square = V, EV

assembly = Lar.Struct([
    Lar.Struct([Lar.t(rand(0:0.05:0.5), rand(0:0.05:0.5)), Lar.r(rand(0:0.05:0.5)), square]) for i in 1:T
])

V, EV = Lar.struct2lar(assembly)
GL.VIEW([GL.GLGrid(V, EV, GL.COLORS[1], 1), GL.GLFrame2]);

W, copEV, copFE, boolmatrix = LARgenerators.bool2d(assembly)

function meanElapsedTime(T) #0.040
    a = []
    for _ in 1:T
        push!(a, @elapsed bool2dEx(assembly))
    end
    print(mean(a))
end

meanElapsedTime(100)

@btime LARgenerators.bool2d($assembly) #7.790 s
@btime bool2dEx($assembly) #8.803 s

V, EV = Lar.struct2lar(assembly)
cop_EW = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells)) #da -1 a 1 CPD9.characteristicMatrix non ha segno
Z = convert(Lar.Points, V')
W, copEV, copFE = Lar.planar_arrangement(Z::Lar.Points, cop_EW::Lar.ChainOp)
@btime Lar.planar_arrangement($Z, $cop_EW) #6.131 s

innerpoints = LARgenerators.internalpoints2d(W, copEV, copFE)
@btime LARgenerators.internalpoints2d($W, $copEV, $copFE) #84.684 ms
listOfModels = Lar.evalStruct(assembly)
@btime Lar.evalStruct(assembly) #140.400 μs

boolmatrix = BitArray(undef, copFE.m, length(listOfModels))

x = Vector{Any}()

for _ in 1:Threads.nthreads()
    push!(x, Any[])
end

containmenttest = LARgenerators.testinternalpoint2d(listOfModels, x)

function createBoolMatrix()
    for (k, point) in collect(enumerate(innerpoints))
        cells = containmenttest(point)
        for l in cells
            boolmatrix[k, l] = 1
        end
    end
end

createBoolMatrix()
@btime createBoolMatrix()#163 ms
