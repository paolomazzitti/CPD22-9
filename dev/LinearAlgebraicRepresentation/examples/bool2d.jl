using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL

T=3
n,m = 1,1
V,(VV,EV,FV) = Lar.cuboidGrid([n,m],true)
square = V,EV
q = sqrt(2)
assembly = Lar.Struct([
    Lar.Struct([ Lar.t( rand(0:0.1:0.5), rand(0:0.1:0.5) ), Lar.r(rand(0:0.1:2)), square ]) for i in 1:T
])

V,EV = Lar.struct2lar(assembly)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);

W, copEV, copFE, boolmatrix = Lar.bool2d(assembly)

A = boolmatrix[:,1]
B = boolmatrix[:,2]
C = boolmatrix[:,3]

AorB = A .| B .| C
AandB = A .& B .& C
AxorB = A .⊻ B .⊻ C

union = Matrix(copFE)' * Int.(AorB)
intersection = Matrix(copFE)' * Int.(AandB)
xor = Matrix(copFE)' * AxorB

V = convert(Lar.Points,W')
EV = Lar.cop2lar(copEV)
EVor = [ev for (k,ev) in enumerate(EV) if abs(union[k])==1 ]
EVand = [ev for (k,ev) in enumerate(EV) if abs(intersection[k])==1 ]
EVxor = [ev for (k,ev) in enumerate(EV) if abs(xor[k])==1 ]
GL.VIEW([ GL.GLGrid(V,EVor, GL.COLORS[1],1), GL.GLFrame2 ]);
GL.VIEW([ GL.GLGrid(V,EVand, GL.COLORS[1],1), GL.GLFrame2 ]);
GL.VIEW([ GL.GLGrid(V,EVxor, GL.COLORS[1],1), GL.GLFrame2 ]);

model = (V,[VV, EVxor])
GL.VIEW( GL.numbering(.5)( model,GL.COLORS[1],0.1 ) );
