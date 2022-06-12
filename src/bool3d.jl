using SparseArrays, LARgenerators, Base.Threads, SparseArrays, IntervalTrees, LinearAlgebra
import ViewerGL as GL
import LinearAlgebraicRepresentation as Lar

#
#	Method to compute an internal point to a polyhedron.
#	----------------------------------------------------
#
#	1. Take two of points close to the opposite sides of any face of a polyhedron, e.g., the first face.
#	2. For each of the two points compute the intersections of a (vertical) ray with the planes (of the faces) intersected by the ray (positive direction of the half-line).
#	3. Transform each such plane and face (and intersection point) to 2D space.
#	4. Test for point-in-polygon intersection.
#	5. Compute the parity of the intersection points for each ray.
#	6. Invariant:  if one is even; the other is odd.
#	7. The initial point with odd number of intersection points is interior to the polyhedron. The other is exterior.
#

"""
	spaceindex(point3d)(model)

Compute the set of face boxes of possible intersection with a point-ray.
Work in 3D, where the ray direction is parallel to the z-axis.
Return an array of indices of face.

# Usage
```jldoctest
julia> V,(VV,EV,FV,CV) = Lar.cuboidGrid([1,1,1],true)

julia> spaceindex([.5,.5,.5])((V,FV))
3-element Array{Int64,1}:
 5
 6
```
"""
function spaceindex(point3d::Array{Float64,1})::Function
    function spaceindex0(model::Lar.LAR)::Array{Int,1}
        V, CV = copy(model[1]), copy(model[2])
        V = [V point3d]
        dim, idx = size(V)
        push!(CV, [idx, idx, idx])
        cellpoints = [V[:, CV[k]]::Lar.Points for k = 1:length(CV)]
        #----------------------------------------------------------
        bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
        xboxdict = Lar.coordintervals(1, bboxes)
        yboxdict = Lar.coordintervals(2, bboxes)
        # xs,ys are IntervalTree type
        xs = IntervalTrees.IntervalMap{Float64,Array}()
        for (key, boxset) in xboxdict
            xs[tuple(key...)] = boxset
        end
        ys = IntervalTrees.IntervalMap{Float64,Array}()
        for (key, boxset) in yboxdict
            ys[tuple(key...)] = boxset
        end
        xcovers = Lar.boxcovering(bboxes, 1, xs)
        ycovers = Lar.boxcovering(bboxes, 2, ys)
        covers = [intersect(pair...) for pair in zip(xcovers, ycovers)]

        # remove each cell from its cover
        pointcover = setdiff(covers[end], [idx + 1])
        return pointcover[1:end-1]
    end
    return spaceindex0
end

"""
	rayintersection(point3d::Array{Float64})(V,FV,face::Int)

# Explanation
Compute the intersection point of the vertical line through `point3d` w `face`.
If the face is parallel to `z axis` return `false`.
# Usage
```jldoctest
julia> V,(VV,EV,FV,CV) = Lar.simplex(3,true);

julia> V
3×4 Array{Float64,2}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0

julia> FV
4-element Array{Array{Int64,1},1}:
 [1, 2, 3]
 [1, 2, 4]
 [1, 3, 4]
 [2, 3, 4]

 julia> Lar.rayintersection([.333,.333,0])(V,FV,4)
 3-element Array{Float64,1}:
  0.333
  0.333
  0.3340000000000001
```
"""
function rayintersection(point3d)
    function rayintersection0(V, FV, face::Int)
        l0, l = point3d, [0, 0, 1.0]
        ps = V[:, FV[face]]  # face points
        p0 = ps[:, 1]
        v1, v2 = ps[:, 2] - p0, ps[:, 3] - p0
        n = LinearAlgebra.normalize(cross(v1, v2))

        denom = dot(n, l)
        if (abs(denom) > 1e-8) #1e-6
            p0l0 = p0 - l0
            t = dot(p0l0, n) / denom
            if t > 0
                return l0 + t * l
            end
        else
            #error("ray and face are parallel")
            return false
        end
    end
    return rayintersection0
end


"""
	planemap(V,copEV,copFE,face)(point)

# Arguments

- `V::Matrix{Float64}`
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`
- `face::Int64`

# Explanation

Tranform the 3D face and the 3D point in their homologous 2D, in order to test for containment.
"""
function planemap(V, copEV, copFE, face)
    fv, edges = Lar.vcycle(copEV, copFE, face)
    function planemap0(point)
        vs = V[:, fv]
        point = point .- vs[:, 1]
        vs = vs .- vs[:, 1]
        u, v = edges[1]
        z, w = [[z, w] for (z, w) in edges if z == v][1]
        v1 = vs[:, u] - vs[:, v]
        v2 = vs[:, w] - vs[:, v]
        v3 = cross(v2, v1)
        M = [v1 v2 v3]
        vs = inv(M) * [vs point]
        outvs = vs[1:2, 1:end-1]
        outpoint = vs[1:2, end]
        return outvs, edges, outpoint
    end
    return planemap0
end

"""
    settestpoints(V, FV, Fs, copEV, copFE)

# Arguments
- `V::Matrix{Float64}`
- `FV::Vector{Vector{Int64}}`
- `Fs::Vector{Int64}`
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`

# Return
- `(ptest1, ptest2)::Tuple{Vector{Float64}, Vector{Float64}}`
"""
function settestpoints(V, FV, Fs, copEV, copFE)
    fdict = Dict(zip(Fs, 1:length(Fs)))
    countKey = 1
    f = Fs[1]
    e = findnz(copFE[f, :])[1][1] 
    f1 = findnz(copFE[:, e])[1][countKey]
    while !haskey(fdict, f1)
        countKey = countKey + 1
        f1 = findnz(copFE[:, e])[1][countKey]
    end
    f2 = findnz(copFE[:, e])[1][countKey+1]
    while !haskey(fdict, f2)
        countKey = countKey + 1
        f2 = findnz(copFE[:, e])[1][countKey]
    end
    v1, v2 = findnz(copEV[e, :])[1] 
    V1 = FV[fdict[f1]]
    V2 = FV[fdict[f2]]
    v1, v2 = intersect(V1, V2)
    t1 = V[:, v1], V[:, v2], V[:, [v for v in V1 if v ≠ v1 && v ≠ v2][1]]
    t2 = V[:, v2], V[:, v1], V[:, [v for v in V2 if v ≠ v1 && v ≠ v2][1]]
    n1 = LinearAlgebra.normalize(cross(t1[2] - t1[1], t1[3] - t1[1]))
    n2 = LinearAlgebra.normalize(cross(t2[2] - t2[1], t2[3] - t2[1]))
    p0 = (V[:, v1] + V[:, v2]) ./ 2
    n = n1 + n2
    ϵ = 1.0e-4
    ptest1 = p0 + ϵ * n
    ptest2 = p0 - ϵ * n
    return ptest1, ptest2
end



"""
	testinternalpoint(V, EV, FV)

# Arguments
- `V::Matrix{Float64}`
- `EV::Vector{Vector{Int64}}`
- `FV::Vector{Vector{Int64}}`

# Return
- `intersectedfaces::Vector{Int64}`
"""
function testinternalpoint(V, EV, FV)
    copEV = Lar.lar2cop(EV)
    copFV = Lar.lar2cop(FV)
    copFE = copFV * copEV'
    I, J, Val = findnz(copFE)
    triple = zip([(i, j, 1) for (i, j, v) in zip(I, J, Val) if v == 2]...)
    I, J, Val = map(collect, triple)
    Val = convert(Array{Int8,1}, Val)
    copFE = sparse(I, J, Val)
    function testinternalpoint0(testpoint)

        intersectedfaces = []

        faces = spaceindex(testpoint)((V, FV))
        depot = []

        for face in faces
            value = LARgenerators.rayintersection(testpoint)(V, FV, face)
            if typeof(value) == Array{Float64,1}
                push!(depot, (face, value))
            end
        end

        for (face, point3d) in depot
            vs, edges, point2d = LARgenerators.planemap(V, copEV, copFE, face)(point3d)
            classify = LARgenerators.pointInPolygonClassificationBool(vs, edges)
            inOut = classify(point2d)
            if inOut != "p_out"
                push!(intersectedfaces, face)
            end
        end
        return intersectedfaces
    end
    return testinternalpoint0
end



"""
    getinternalpoint(V,EV,FV,Fs, copEV,copFE,z,ptest1,ptest2)

# Arguments
- `V::Matrix{Float64}`
- `EV::Vector{Vector{Int64}}`
- `FV::Vector{Vector{Int64}}`
- `Fs::Vector{Int64}`
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`
- `z::Vector{Any}`
- `ptest1::Vector{Int64}`
- `ptest2::Vector{Int64}`

# Return
- `ptest1::Vector{Int64}` or `ptest2::Vector{Int64}`

# Explanation
- For each test point compute the face planes intersected by vertical ray
- Transform each plane in 2D and look whether the intersection point is internal
- Return the test point with odd numeber of ray intersections
"""
function getinternalpoint(V, FV, Fs, copEV, copFE, z, ptest1, ptest2)

    for (f, face) in collect(enumerate(Fs))
        ret1 = LARgenerators.rayintersection(ptest1)(V, FV, f)
        if typeof(ret1) == Array{Float64,1}
            push!(z[threadid()], (face, ret1))
        end
    end

    k1 = 0
    for (face, point3d) in z[threadid()]
        vs, edges, point2d = planemap(V, copEV, copFE, face)(point3d)
        classify = LARgenerators.pointInPolygonClassificationBool(vs, edges)
        inOut = classify(point2d)
        if inOut == "p_in"
            k1 += 1
        end
    end
    if k1 % 2 == 1
        return ptest1
    else
        return ptest2
    end
end

"""
    chainbasis2solids(copEV, copFE, copCF)

# Arguments
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`
- `copCF::SparseMatrixCSC{Int8, Int64}`

# Return
- (`FVs, CF)::Tuple{Vector{Vector{Vector{Int64}}}, Vector{Vector{Int64}}}`
"""
function chainbasis2solids(copEV, copFE, copCF)
    CF = [findnz(copCF[k, :])[1] for k = 1:copCF.m]
    FE = [findnz(copFE[k, :])[1] for k = 1:copFE.m]
    EV = [findnz(copEV[k, :])[1] for k = 1:copEV.m]

    FVs = Vector{Vector{Vector{Int64}}}(undef, copCF.m)

    @threads for k = 1:copCF.m    #= @threads =#
        FVs[k] = [collect(Set(GL.Cat([EV[e] for e in FE[f]]))) for f in CF[k]]
    end

    return FVs, CF
end

"""
    internalpoints(V, copEV, copFE, copCF)
# Arguments
- `V::Matrix{Float64}`
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`
- `copCF::SparseMatrixCSC{Int8, Int64}`

# Return
- `innerpoints::Vector{Any}`

# Explanation
- Transform each 3-cell in a solid
- Compute, for each 3-cell, one *internal point*
"""
function internalpoints(V, copEV, copFE, copCF)

    FVs, CF = LARgenerators.chainbasis2solids(copEV, copFE, copCF)

    innerpoints = Array{Any}(undef, length(FVs))
    ptestArray = Array{Any}(undef, length(FVs))

    z = Vector{Any}()

    for _ in 1:Threads.nthreads()
        push!(z, Any[])
    end

    for l in 1:length(FVs)
        ptestArray[l] = LARgenerators.settestpoints(V, FVs[l], CF[l], copEV, copFE)
    end

    for k in 1:length(FVs)    #= @threads =#
        points = LARgenerators.getinternalpoint(V, FVs[k], CF[k], copEV, copFE, z, ptestArray[k][1], ptestArray[k][2])
        innerpoints[k] = points
        z[threadid()] = []
    end

    return innerpoints
end

"""
    bool3d(assembly)

# Arguments
- `assembly::LinearAlgebraicRepresentation.Struct`

# Returns
- `W::Matrix{Float64}`
- `copEV::SparseMatrixCSC{Int8, Int64}`
- `copFE::SparseMatrixCSC{Int8, Int64}`
- `copCF::SparseMatrixCSC{Int8, Int64}`
- `boolmatrix::BitMatrix`

# Usage

```julia
using LARgenerators, SparseArrays, Base.Threads
import LinearAlgebraicRepresentation as Lar

n, m, p = 1, 1, 1
V, (VV, EV, FV, CV) = Lar.cuboidGrid([n, m, p], true)
cube = V, FV, EV

assembly = Lar.Struct([ cube,
    Lar.t(.3,.4,.25), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), cube,
    Lar.t(-.2,.4,-.2), Lar.r(0,pi/5,0), Lar.r(0,pi/12,0), cube
])

W, copEV, copFE, copCF, boolmatrix = LARgenerators.bool3d(assembly)
```

# Explanation 
1. After the arrangement, extract all the d-cells from (d-1)-coboundary as isolated polyhedra.
2. Then compute a single interior point for each of them.
3. Then compare each such point against all input boundaries, in order to compute those which it was interior to. Extend this point membership as 3-cell containment within the relative input solids.
4. The point membership with a boundary consists in the parity count of the intersection points of a vertical ray starting at the test point, with the boundary surface.
"""
function bool3d(assembly)

    V, FV, EV = Lar.struct2lar(assembly)
    cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells))
    cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells)
    W = convert(Lar.Points, V')

    Z::Matrix{Float64}, copEV::SparseMatrixCSC{Int8,Int64}, copFE::SparseMatrixCSC{Int8,Int64}, copCF::SparseMatrixCSC{Int8,Int64} = Lar.space_arrangement(W, cop_EV, cop_FE)
    W = convert(Lar.Points, Z')

    innerpoints = internalpoints(W, copEV, copFE, copCF[2:end, :])

    listOfModels = Lar.evalStruct(assembly)
    inputfacenumbers = [length(listOfModels[k][2]) for k = 1:length(listOfModels)]
    cumulative = cumsum([0; inputfacenumbers]) .+ 1
    fspans = collect(zip(cumulative[1:end-1], cumulative[2:end] .- 1))
    span(h) = [j for j = 1:length(fspans) if fspans[j][1] <= h <= fspans[j][2]]

    containmenttest = testinternalpoint(V, EV, FV)
    boolmatrix = BitArray(undef, length(innerpoints) + 1, length(fspans) + 1)
    boolmatrix[1, 1] = 1

    for (k, point) in collect(enumerate(innerpoints))
        cells = containmenttest(point)
        rows = [span(h) for h in cells]
        for l in GL.Cat(rows)
            boolmatrix[k+1, l+1] = 1
        end
    end
    return W, copEV, copFE, copCF, boolmatrix
end