gap> START_TEST( "LargeAlgebras.tst");

#
gap> QG:=GroupRing(Rationals,SmallGroup(1920,236000));;
gap> W:=WedderburnDecompositionInfo(QG);;
gap> A:=W[Length(W)];
[ 1, NF(5,[ 1, 4 ]), 120, [ [ 2, 11, 0 ], [ 2, 19, 60 ], [ 2, 29, 0 ], [ 2, 31, 0 ] ], [ [ 0, 0, 0 ], [ 60, 0 ], [ 0 ] ] ]
gap> LocalIndicesOfCyclotomicAlgebra(A);
[ [ infinity, 2 ] ]

#
gap> QG:=GroupRing(Rationals,SmallGroup(1920,236001));;
gap> W:=WedderburnDecompositionInfo(QG);;
gap> A:=W[Length(W)];
[ 1, NF(5,[ 1, 4 ]), 120, [ [ 2, 11, 0 ], [ 2, 19, 0 ], [ 2, 29, 60 ], [ 2, 31, 0 ] ], [ [ 0, 0, 0 ], [ 60, 0 ], [ 0 ] ] ]
gap> LocalIndicesOfCyclotomicAlgebra(A);
[  ]

#
gap> QG:=GroupRing(Rationals,SmallGroup(1920,236002));;
gap> W:=WedderburnDecompositionInfo(QG);;
gap> A:=W[Length(W)];
[ 1, NF(5,[ 1, 4 ]), 120, [ [ 2, 11, 0 ], [ 2, 19, 60 ], [ 2, 29, 0 ], [ 2, 31, 0 ] ], [ [ 60, 0, 60 ], [ 0, 60 ], [ 0 ] ] 
 ]
gap> LocalIndicesOfCyclotomicAlgebra(A);
[  ]

#
gap> QG:=GroupRing(Rationals,SmallGroup(1920,236003));;
gap> W:=WedderburnDecompositionInfo(QG);;
gap> A:=W[Length(W)];
[ 1, NF(5,[ 1, 4 ]), 120, [ [ 2, 11, 0 ], [ 2, 19, 0 ], [ 2, 29, 60 ], [ 2, 31, 0 ] ], [ [ 0, 60, 0 ], [ 0, 0 ], [ 0 ] ] ]
gap> LocalIndicesOfCyclotomicAlgebra(A);
[ [ 5, 2 ] ]

#
gap> QG:=GroupRing(Rationals,SmallGroup(1920,236007));;
gap> W:=WedderburnDecompositionInfo(QG);;
gap> A:=W[Length(W)];
[ 1, NF(5,[ 1, 4 ]), 120, [ [ 2, 11, 0 ], [ 2, 19, 0 ], [ 2, 29, 60 ], [ 2, 31, 0 ] ], [ [ 0, 0, 0 ], [ 60, 0 ], [ 60 ] ] ]
gap> LocalIndicesOfCyclotomicAlgebra(A);
[ [ 5, 2 ] ]

#
gap> QG:=GroupRing(Rationals,SmallGroup(1920,236008));;
gap> W:=WedderburnDecompositionInfo(QG);;
gap> A:=W[Length(W)];
[ 1, NF(5,[ 1, 4 ]), 120, [ [ 2, 11, 0 ], [ 2, 19, 0 ], [ 2, 29, 0 ], [ 2, 31, 0 ] ], [ [ 0, 0, 0 ], [ 0, 60 ], [ 0 ] ] ]
gap> LocalIndicesOfCyclotomicAlgebra(A);
[ [ infinity, 2 ] ]

#
gap> QG:=GroupRing(Rationals,SmallGroup(1920,236009));;
gap> W:=WedderburnDecompositionInfo(QG);;
gap> A:=W[Length(W)];
[ 1, NF(5,[ 1, 4 ]), 120, [ [ 2, 11, 0 ], [ 2, 19, 0 ], [ 2, 29, 60 ], [ 2, 31, 0 ] ], [ [ 0, 60, 0 ], [ 0, 0 ], [ 60 ] ] ]
gap> LocalIndicesOfCyclotomicAlgebra(A);
[  ]

#
gap> STOP_TEST( "LargeAlgebras.tst", 1 );