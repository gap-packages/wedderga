gap> START_TEST( "others.tst");
gap> G:=SmallGroup(1000,86);                 
<pc group of size 1000 with 6 generators>
gap>  QG:=GroupRing(Rationals,G);
<algebra-with-one over Rationals, with 6 generators>
gap> ShodaPairsAndIdempotents(QG).ShodaPairs;
[ [ Group([ f1, f2, f3, f4, f5, f6 ]), Group([ f6, f4, f5, f3, f2, f1 ]) ], 
  [ Group([ f1, f2, f3, f4, f5, f6 ]), Group([ f6, f4, f5, f3, f2 ]) ], 
  [ Group([ f1, f2, f3, f4, f5, f6 ]), Group([ f6, f4, f5, f3 ]) ], 
  [ Group([ f1, f2, f3, f4, f5, f6 ]), Group([ f6, f4, f5 ]) ], 
  [ Group([ f5, f6, f4 ]), Group([ f6, f5 ]) ], 
  [ Group([ f4^4, f6^4, f5 ]), Group([ f6, f4 ]) ], 
  [ Group([ f4^4*f5*f6^4, f6^4, f5 ]), Group([ f6, f4*f5^4 ]) ], 
  [ Group([ f3, f4^4*f5*f6^2, f6 ]), Group([ f3, f4*f5^4*f6^2 ]) ], 
  [ Group([ f4^4*f5*f6^4, f3*f5^4, f3*f5^4*f6 ]), Group([ f4*f5^4 ]) ] ]
gap> STOP_TEST( "idempot.tst", 1 );
