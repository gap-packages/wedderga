gap> START_TEST( "div-alg.tst");

#
gap> A := [1,Rationals,120,
>         [[2,31,0],[2,61,0],[2,41,0],[4,97,0]],[[60,0,0],[0,0],[60]]];;
gap> ReducingCyclotomicAlgebra(A);
[ 4, Rationals, 30, [ [ 2, 11, 0 ], [ 4, 7, 0 ] ], [ [ 15 ] ] ]

#
gap> A := [1,Rationals,120,
>         [[2,31,0],[2,61,0],[4,73,0],[2,41,0]],[[60,0,0],[0,0],[60]]];;
gap> ReducingCyclotomicAlgebra(A);
[ 4, Rationals, 30, [ [ 4, 13, 0 ], [ 2, 11, 0 ] ], [ [ 15 ] ] ]

#
gap> A1:=[ 1, Rationals, 8, [ [ 2, 7, 0 ], [ 2, 5, 0 ] ], [ [ 4 ] ] ];
[ 1, Rationals, 8, [ [ 2, 7, 0 ], [ 2, 5, 0 ] ], [ [ 4 ] ] ]
gap> ReducingCyclotomicAlgebra(A1);
fail
gap> SchurIndex(A1);
1

#
gap> A2:=[ 1, Rationals, 30, [ [ 2, 11, 0 ], [ 4, 7, 0 ] ], [ [ 15 ] ] ];
[ 1, Rationals, 30, [ [ 2, 11, 0 ], [ 4, 7, 0 ] ], [ [ 15 ] ] ]
gap> ReducingCyclotomicAlgebra(A2);
fail
gap> SchurIndex(A2);
2

# Example from PR #64
gap> A:=[1,Rationals,120,[[2,31,0],[2,61,0],[2,41,0],[4,97,0]],[[60,0,0],[0,0],[60]]];;
gap> LocalIndicesOfCyclotomicAlgebra(A);
[ [ 3, 2 ], [ 5, 2 ] ]
gap> B:=ReducingCyclotomicAlgebra(A);
[ 4, Rationals, 30, [ [ 2, 11, 0 ], [ 4, 7, 0 ] ], [ [ 15 ] ] ]
gap> LocalIndicesOfCyclotomicAlgebra(B);
[ [ 3, 2 ], [ 5, 2 ] ]
gap> SchurIndex(A);
2

# Example from PR #78
gap> QG:=GroupRing(Rationals,SmallGroup(672,622));
<algebra-with-one over Rationals, with 7 generators>
gap> wd:=WedderburnDecompositionInfo(QG);;
gap> A:=wd[28];;
gap> A = KillingCocycle(A);
true
gap> LocalIndicesOfCyclotomicAlgebra(A);
[ [ 7, 2 ] ]

# Length 4 example from PR #80
gap> DecomposeCyclotomicAlgebra([1,NF(7,[1,2,4]),6,[3,2,3]]);
[ NF(7,[ 1, 2, 4 ]),NF(21,[ 1, 4, 16 ]),[ -1 ] ]

# Length 5 example from PR #80 with all zeroes in A[5]
gap> A:=[1,NF(7,[1,6]),84,[[2,13,42],[2,29,63],[2,43,0]],[[0,0],[0]]];
[ 1, NF(7,[ 1, 6 ]), 84, [ [ 2, 13, 42 ], [ 2, 29, 63 ], [ 2, 43, 0 ] ], [ [ 0, 0 ], [ 0 ] ] ]
gap> DecomposeCyclotomicAlgebra(A);
[ [ NF(7,[ 1, 6 ]), CF(7), [ -1 ] ],
[ NF(7,[ 1, 6 ]), NF(21,[ 1, 13 ]), [ -E(4) ] ],
[ NF(7,[ 1, 6 ]), NF(28,[ 1, 13 ]), [ 1 ] ] ]

# Some special cases in DecomposeCyclotomicAlgebra (PR #81)
gap> QG:=GroupRing(Rationals,SmallGroup(240,96));;
gap> W:=WedderburnDecompositionInfo(QG);;
gap> Length(W);
18
gap> A:=W[Length(W)];
[ 1, Rationals, 30, [ [ 2, 11, 0 ], [ 4, 7, 0 ] ], [ [ 15 ] ] ]
gap> DecomposeCyclotomicAlgebra(A);
[ [ Rationals, CF(3), [ 1 ] ], [ Rationals, CF(5), [ 9 ] ] ]

# LocalIndicesOfRationalSymbolAlgebra (PR #83)
gap> List([-2..3],a->LocalIndicesOfRationalSymbolAlgebra(a,-2));
[ fail, fail, fail, fail, fail, fail ]
gap> List([-2..3],a->LocalIndicesOfRationalSymbolAlgebra(a,-1));
[ fail, [ [ infinity, 2 ], [ 2, 2 ] ], fail, fail, [  ], [ [ 2, 2 ], [ 3, 2 ] ] ]
gap> List([-2..3],a->LocalIndicesOfRationalSymbolAlgebra(a,2));
[ fail, [  ], fail, fail, [  ], [ [ 2, 2 ], [ 3, 2 ] ] ]
gap> List([-2..3],a->LocalIndicesOfRationalSymbolAlgebra(a,3));
[ fail, [ [ 2, 2 ], [ 3, 2 ] ], fail, fail, [ [ 2, 2 ], [ 3, 2 ] ], [ [ 2, 2 ], [ 3, 2 ] ] ]
gap> List([-2..3],a->LocalIndicesOfRationalSymbolAlgebra(a,5));
[ fail, [  ], fail, fail, [ [ 2, 2 ], [ 5, 2 ] ], [ [ 3, 2 ], [ 5, 2 ] ] ]
gap> List([-2..3],a->LocalIndicesOfRationalSymbolAlgebra(a,7));
[ fail, [ [ 2, 2 ], [ 7, 2 ] ], fail, fail, [  ], [ [ 2, 2 ], [ 7, 2 ] ] ]
gap> List([-2..3],a->LocalIndicesOfRationalSymbolAlgebra(a,11));
[ fail, [ [ 2, 2 ], [ 11, 2 ] ], fail, fail, [ [ 2, 2 ], [ 11, 2 ] ], [ [ 2, 2 ], [ 3, 2 ] ] ]

# IsDyadicSchurGroup (PR #104)
gap> IsDyadicSchurGroup(SmallGroup(8,4));
true
gap> IsDyadicSchurGroup(SmallGroup(160,208));
true
gap> IsDyadicSchurGroup(SmallGroup(160,84));
false

#
gap> STOP_TEST( "div-alg.tst", 1 );
