using CPD9, ViewerGL, SparseArrays
GL = ViewerGL

T=3
n,m = 1,1
V,(VV,EV,FV) = CPD9.cuboidGrid([n,m],true)
square = V,EV
q = sqrt(2)
assembly = CPD9.Struct([
    CPD9.Struct([ CPD9.t( rand(0:0.1:0.5), rand(0:0.1:0.5) ), CPD9.r(rand(0:0.1:2)), square ]) for i in 1:T
])

V,EV = CPD9.struct2lar(assembly)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);

W, copEV, copFE, boolmatrix = CPD9.bool2d(assembly)

A = boolmatrix[:,1]
B = boolmatrix[:,2]
C = boolmatrix[:,3]

AorB = A .| B .| C
AandB = A .& B .& C
AxorB = A .⊻ B .⊻ C

union = Matrix(copFE)' * Int.(AorB)
intersection = Matrix(copFE)' * Int.(AandB)
xor = Matrix(copFE)' * AxorB

V = convert(CPD9.Points,W')
EV = CPD9.cop2lar(copEV)
EVor = [ev for (k,ev) in enumerate(EV) if abs(union[k])==1 ]
EVand = [ev for (k,ev) in enumerate(EV) if abs(intersection[k])==1 ]
EVxor = [ev for (k,ev) in enumerate(EV) if abs(xor[k])==1 ]
GL.VIEW([ GL.GLGrid(V,EVor, GL.COLORS[1],1), GL.GLFrame2 ]);
GL.VIEW([ GL.GLGrid(V,EVand, GL.COLORS[1],1), GL.GLFrame2 ]);
GL.VIEW([ GL.GLGrid(V,EVxor, GL.COLORS[1],1), GL.GLFrame2 ]);

model = (V,[VV, EVxor])
GL.VIEW( GL.numbering(.5)( model,GL.COLORS[1],0.1 ) );
