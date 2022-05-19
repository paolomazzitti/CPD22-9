using CPD9, ViewerGL, SparseArrays
GL = ViewerGL

wall = CPD9.svg2lar("dev/CPD9/test/svg/assembly/wall.svg", flag=false)
openings = CPD9.svg2lar("dev/CPD9/test/svg/assembly/openings.svg", flag=false)
rectangle = CPD9.svg2lar("dev/CPD9/test/svg/assembly/rectangle.svg", flag=false)
box = CPD9.svg2lar("dev/CPD9/test/svg/assembly/box.svg", flag=false)

assembly = CPD9.Struct([wall, openings, rectangle, box])
V,EV = CPD9.struct2lar(assembly)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);

assembly = CPD9.Struct([wall, openings, rectangle])
V,EV = CPD9.boolops(assembly, :|)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
V,FVs,EVs = CPD9.arrange2D(V,EV)
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));

assembly = CPD9.Struct([wall, openings, rectangle])
V,EV = CPD9.boolops(assembly, :-)
diff = CPD9.Struct([ (V,EV) ])
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
V,FVs,EVs = CPD9.arrange2D(V,EV)
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));

assembly = CPD9.Struct([wall, openings, box])
V,EV = CPD9.boolops(assembly, :&)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
V,FVs,EVs = CPD9.arrange2D(V,EV)
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));

inner = CPD9.Struct([ box,diff ])
V,EV = CPD9.boolops( inner, :-)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
V,FVs,EVs = CPD9.arrange2D(V,EV)
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));

W, copEV, copFE, boolmatrix = CPD9.bool2d(inner)
FVs = CPD9.triangulate2D(W, [copEV, copFE])
V = convert(CPD9.Points, W')
GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
