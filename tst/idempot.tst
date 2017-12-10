gap> START_TEST( "idempot.tst");
gap> G:=SymmetricGroup(4);
Sym( [ 1 .. 4 ] )
gap> S:=StrongShodaPairs(G);
[ [ Sym( [ 1 .. 4 ] ), Sym( [ 1 .. 4 ] ) ], 
  [ Sym( [ 1 .. 4 ] ), Group([ (2,4,3), (1,4)(2,3), (1,3)(2,4) ]) ], 
  [ Group([ (2,4,3), (1,4)(2,3), (1,3)(2,4) ]), Group([ (1,4)(2,3), (1,3)
      (2,4) ]) ], [ Group([ (1,3,2,4), (3,4) ]), Group([ (3,4), (1,2)(3,4) ]) 
     ], [ Group([ (3,4), (1,3,2,4) ]), Group([ (1,3,2,4), (1,2)(3,4) ]) ] ]
gap> K:=S[5][1];;
gap> H:=S[5][2];;
gap> QG:=GroupRing(Rationals,G);
<algebra-with-one over Rationals, with 2 generators>
gap> CentralElementBySubgroups(QG,K,H);
(3/8)*()+(-1/8)*(3,4)+(-1/8)*(2,3)+(-1/8)*(2,4)+(-1/8)*(1,2)+(-1/8)*(1,2)
(3,4)+(1/8)*(1,2,3,4)+(1/8)*(1,2,4,3)+(1/8)*(1,3,4,2)+(-1/8)*(1,3)+(-1/8)*
(1,3)(2,4)+(1/8)*(1,3,2,4)+(1/8)*(1,4,3,2)+(-1/8)*(1,4)+(1/8)*(1,4,2,3)+(-1/
8)*(1,4)(2,3)
gap> STOP_TEST( "idempot.tst", 1 );
