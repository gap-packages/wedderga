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

#
gap> STOP_TEST( "div-alg.tst", 1 );
