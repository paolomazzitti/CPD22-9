import LinearAlgebraicRepresentation as Lar
using Base.Threads

"""
    setTileBool(box)(point)

Returns the tile code of an edge.

# Arguments
- `box::Tuple{Float64, Float64, Float64, Float64}`
- `point::Vector{Float64}`

# Return
- `tileCodeBool::Int64`
"""
@inline function setTileBool(box)
    b1, b2, b3, b4 = box
    function tileCodeBool(point)
        x, y = point
        code = 0
        if y > b1
            code = code | 1
        end
        if y < b2
            code = code | 2
        end
        if x > b3
            code = code | 4
        end
        if x < b4
            code = code | 8
        end
        return code
    end
    return tileCodeBool
end

"""
    pointInPolygonClassificationBool(V, EV)(pnt)

Determines whether a point is inside or outside a polygon. It is a specific version used for the bool2d function.

# Arguments
- `V::Adjoint{Float64, Matrix{Float64}}`
- `EV::Vector{Vector{Int64}}`
- `pnt::Vector{Float64}`

# Return
- `p_in::String` or `p_out::String`
"""
function pointInPolygonClassificationBool(V, EV)
    function pointInPolygonClassification0Bool(pnt)
        x, y = pnt
        xmin, xmax, ymin, ymax = x, x, y, y
        tilecodeBool = setTileBool([ymax, ymin, xmax, xmin])
        count = 0
        for edge in EV
            p1, p2 = V[:, edge[1]], V[:, edge[2]]
            c1, c2 = tilecodeBool(p1), tilecodeBool(p2)
            (x1, y1), (x2, y2) = p1, p2
            c_edge, c_int = c1 ⊻ c2, c1 & c2

            if c_edge == 3 && c_int == 4
                count += 1
            elseif c_edge == 15
                x_int = ((y - y2) * (x1 - x2) / (y1 - y2)) + x2
                if x_int > x
                    count += 1
                end
            end
        end
        if (round(count) % 2) == 1
            return "p_in"
        else
            return "p_out"
        end
    end
    return pointInPolygonClassification0Bool
end

"""
    settestpoints2d(W, f, copEV, copFE)

Given a polygon, it returns two points, one inside and one outside it.

# Arguments
- `W::Matrix{Float64}`
- `f::Int64`
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`

# Return
- `(ptest1, ptest2)::Tuple{Vector{Float64}, Vector{Float64}}`
"""
@inline function settestpoints2d(W, f, copEV, copFE)
    e = findnz(copFE[f, :])[1][1]
    v1, v2 = findnz(copEV[e, :])[1]
    t = W[v2, :] - W[v1, :]
    n = [-t[2], t[1]]
    p0 = (W[v1, :] + W[v2, :]) ./ 2
    ϵ = 1.0e-4
    ptest1 = p0 + ϵ * n
    ptest2 = p0 - ϵ * n
    return ptest1, ptest2
end

"""
    getinternalpoint2d(W, f, copEV, copFE)

Returns an interior point in the polygon. 

# Arguments

- `W::Matrix{Float64}`
- `f::Int64`
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`

# Return

- `ptest1::Vector{Float64}` or `ptest2::Vector{Float64}`
"""
function getinternalpoint2d(W, f, copEV, copFE)
    ptest1, ptest2 = settestpoints2d(W, f, copEV, copFE)
    edges = [findnz(copEV[e, :])[1] for e in findnz(copFE[f, :])[1]]
    V = W'
    classify = pointInPolygonClassificationBool(V, edges)
    if classify(ptest1) == "p_in"
        return ptest1
    else
        return ptest2
    end
end

"""
    internalpoints2d(W, copEV, copFE)

Returns, for each atom of the decomposition, an internal point.

# Arguments
- `W::Matrix{Float64}`
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`

# Return
- `internalpoints::Vector{Vector{Float64}}`
"""
function internalpoints2d(W, copEV, copFE)
    internalpoints = Vector{Vector{Float64}}(undef, copFE.m)
    for f = 1:copFE.m
        internalpoint = getinternalpoint2d(W, f, copEV, copFE)
        @inbounds internalpoints[f] = internalpoint
    end
    return internalpoints
end

"""
	testinternalpoint2d(listOfModels)(testpoint)

Given a point, it returns the list of models containing it.
    
# Arguments
- `listOfModels::Vector{Any}`
- `testpoint::Vector{Float64}`

# Return
- `intersectedFaces: Vector{Int64}`
"""
function testinternalpoint2d(listOfModels)
    function testinternalpoint0(testpoint)

        intersectedFaces = Vector{Int64}()

        for (k, model) in enumerate(listOfModels)
            verts, edges = model
            classify = pointInPolygonClassificationBool(verts, edges)
            inOut = classify(testpoint)
            if inOut == "p_in"
                push!(intersectedFaces, k)
            end
        end

        return intersectedFaces
    end
    return testinternalpoint0
end

"""
    bool2d(assembly)

# Arguments
- `assembly::LinearAlgebraicRepresentation.Struct`
    
# Return
- `W::Matrix{Float64}`
- `copEV::SparseMatrixCSC{Int8,Int64}`
- `copFE::SparseMatrixCSC{Int8,Int64}`
- `boolmatrix::BitMatrix`

# Usage
```julia
using SparseArrays, LARgenerators, Base.Threads, SparseArrays, IntervalTrees, LinearAlgebra
import LinearAlgebraicRepresentation as Lar

V, (VV, EV, FV) = Lar.cuboidGrid([1, 1], true)
square = V, EV

assembly = Lar.Struct([
    Lar.Struct([Lar.t(0,0), Lar.r(0), square])
    Lar.Struct([Lar.t(0,0.1), Lar.r(0.1), square])
    Lar.Struct([Lar.t(0,0.2), Lar.r(0.2), square])
    Lar.Struct([Lar.t(0,0.3), Lar.r(0.3), square])
    Lar.Struct([Lar.t(0,0.4), Lar.r(0.4), square])
])

W, copEV, copFE, boolmatrix = LARgenerators.bool2d(assembly)
```
"""
function bool2d(assembly)

    V::Matrix{Float64}, EV::Vector{Vector{Int64}} = Lar.struct2lar(assembly)
    cop_EW = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells))
    Z = convert(Lar.Points, V')
    W::Matrix{Float64}, copEV::SparseMatrixCSC{Int8,Int64}, copFE::SparseMatrixCSC{Int8,Int64} = Lar.planar_arrangement(Z::Lar.Points, cop_EW::Lar.ChainOp)
    innerpoints = internalpoints2d(W, copEV, copFE)
    listOfModels = Lar.evalStruct(assembly)
    boolmatrix = BitArray(undef, copFE.m, length(listOfModels))

    containmenttest = testinternalpoint2d(listOfModels)

    for (k, point) in collect(enumerate(innerpoints))
        cells = containmenttest(point)
        for l in cells
            @inbounds boolmatrix[k, l] = 1
        end
    end
    return W, copEV, copFE, boolmatrix
end
