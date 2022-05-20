using Base.Threads
using CPD9

function setTileBool(box)
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

function settestpoints2d(W, f, copEV, copFE)
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

function getinternalpoint2d(W, f, copEV, copFE)
    ptest1, ptest2 = settestpoints2d(W, f, copEV, copFE)
    edges = [findnz(copEV[e, :])[1] for e in findnz(copFE[f, :])[1]]
    V = W'
    classify = CPD9.pointInPolygonClassificationBool(V, edges)
    if classify(ptest1) == "p_in"
        return ptest1
    else
        return ptest2
    end
end

function internalpoints2d(W, copEV, copFE)
    internalpoints = Vector{Vector{Float64}}(undef, copFE.m)
    @threads for f = 1:copFE.m
        internalpoint = getinternalpoint2d(W, f, copEV, copFE)
        @inbounds internalpoints[f] = internalpoint
    end
    return internalpoints
end

"""
	testinternalpoint2d(V::CPD9.Points, EV::CPD9.Cells, FV::CPD9.Cells)

"""
function testinternalpoint2d(listOfModels)
    function testinternalpoint0(testpoint)
        intersectedfaces = []
        for (k, model) in enumerate(listOfModels)
            verts, edges = model
            classify = CPD9.pointInPolygonClassificationBool(verts, edges)
            inOut = classify(testpoint)
            if inOut == "p_in"
                push!(intersectedfaces, k)
            end
        end
        return intersectedfaces
    end
    return testinternalpoint0
end

function bool2d(assembly)

    V, EV = CPD9.struct2lar(assembly)
    cop_EW = convert(CPD9.ChainOp, CPD9.coboundary_0(EV::CPD9.Cells)) #da -1 a 1 CPD9.characteristicMatrix non ha segno
    Z = convert(CPD9.Points, V')
    W, copEV, copFE = CPD9.Arrangement.planar_arrangement(Z::CPD9.Points, cop_EW::CPD9.ChainOp)
    innerpoints = internalpoints2d(W, copEV, copFE)
    listOfModels = CPD9.evalStruct(assembly)
    boolmatrix = BitArray(undef, copFE.m, length(listOfModels))
    containmenttest = testinternalpoint2d(listOfModels)

    Threads.@threads for (k, point) in collect(enumerate(innerpoints))
        cells = containmenttest(point)
        for l in cells
            boolmatrix[k, l] = 1
        end
    end
    return W, copEV, copFE, boolmatrix
end

"""
	boolops( assembly::CPD9.Struct, op::Symbol )

User interface to 2d Boolean ops, where `op` symbol ``in`` {`:|`, `:&`, `:-`}.
Return a pair `(V,EV)`.
# Example
```julia
assembly = CPD9.Struct([wall, openings])
V,EV = boolops(assembly, :-)
```
"""
function boolops(assembly::CPD9.Struct, op::Symbol)
    W, copEV, copFE, boolmatrix = bool2d(assembly)
    boolvars = [boolmatrix[:, k] for k = 1:size(boolmatrix, 2)]
    if eval(op) == -
        solution = boolvars[1] .& .!(.|(boolvars[2:end]...))
    else
        operator = eval(op)
        solution = operator.(boolvars...)
    end
    chain1d = sparse(copFE') * Int.(solution)
    EV = CPD9.cop2lar(copEV)
    EVop = [ev for (k, ev) in enumerate(EV) if abs(chain1d[k]) == 1]
    V = convert(CPD9.Points, W')
    return V, EVop
end

