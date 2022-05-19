module Arrangement
    using CPD9
    using IntervalTrees
    using NearestNeighbors
   #using Triangle
    using Triangulate
	using SparseArrays
	using LinearAlgebra    
	using Distributed    
	Lar = CPD9

    include("./minimal_cycles.jl")
    include("./dimension_travel.jl")
    include("./planar_arrangement.jl")
    include("./spatial_arrangement.jl")

    export planar_arrangement, spatial_arrangement

end