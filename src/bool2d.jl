import LinearAlgebraicRepresentation as Lar
using Base.Threads

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
    classify = pointInPolygonClassificationBool(V, edges)
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
        internalpoints[f] = internalpoint
    end
    return internalpoints
end

"""
	testinternalpoint2d(V::CPD9.Points, EV::CPD9.Cells, FV::CPD9.Cells)
    ciao
"""
function testinternalpoint2d(listOfModels)
    function testinternalpoint0(testpoint)

        intersectedFaces = zeros(Int64, length(listOfModels))

        for (k, model) in enumerate(listOfModels)
            verts, edges = model
            classify = pointInPolygonClassificationBool(verts, edges)
            inOut = classify(testpoint)
            if inOut == "p_in"
                intersectedFaces[k] = k
            end
        end

        return intersectedFaces
    end
    return testinternalpoint0
end
"""
    bool2d(assembly)
    some docstring

"""
function bool2d(assembly)

    main(Vector{Any}(undef, 2))

    V, EV = Lar.struct2lar(assembly)
    cop_EW = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells)) #da -1 a 1 CPD9.characteristicMatrix non ha segno
    Z = convert(Lar.Points, V')
    W, copEV, copFE = Lar.planar_arrangement(Z::Lar.Points, cop_EW::Lar.ChainOp)
    innerpoints = internalpoints2d(W, copEV, copFE)
    listOfModels = Lar.evalStruct(assembly)
    boolmatrix = BitArray(undef, copFE.m, length(listOfModels))

    containmenttest = testinternalpoint2d(listOfModels)

    @threads for (k, point) in collect(enumerate(innerpoints))
        cells = containmenttest(point)
        for l in cells
            if l != 0
                boolmatrix[k, l] = 1
            end
        end
    end
    return W, copEV, copFE, boolmatrix
end
