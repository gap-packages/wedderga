# wedderga, chapter 3
gap> START_TEST( "wedderga03.tst");

# wedderga/doc/SSP.xml:24-52

gap> StrongShodaPairs( SymmetricGroup(4) );
[ [ Sym( [ 1 .. 4 ] ), Sym( [ 1 .. 4 ] ) ], 
  [ Sym( [ 1 .. 4 ] ), Group([ (2,4,3), (1,4)(2,3), (1,3)(2,4)\
 ]) ], 
  [ Group([ (2,4,3), (1,4)(2,3), (1,3)(2,4) ]), Group([ (1,4)(\
2,3), (1,3)
      (2,4) ]) ], [ Group([ (1,3,2,4), (3,4) ]), Group([ (3,4)\
, (1,2)(3,4) ]) 
     ], [ Group([ (3,4), (1,3,2,4) ]), Group([ (1,3,2,4), (1,2\
)(3,4) ]) ] ]

gap> StrongShodaPairs( DihedralGroup(64) );
[ [ <pc group of size 64 with 6 generators>, 
      <pc group of size 64 with 6 generators> ], 
  [ <pc group of size 64 with 6 generators>, 
      Group([ f1*f2*f3*f4*f5*f6, f3, f4, f5, f6 ]) ], 
  [ <pc group of size 64 with 6 generators>, Group([ f2, f3, f\
4, f5, f6 ]) ], 
  [ <pc group of size 64 with 6 generators>, Group([ f1, f3, f\
4, f5, f6 ]) ], 
  [ Group([ f1*f2*f3*f4*f5*f6, f3, f4, f5, f6 ]), 
      Group([ f1*f2*f4*f5*f6, f4, f5, f6 ]) ], 
  [ Group([ f2, f3, f4, f5, f6 ]), Group([ f5, f6 ]) ], 
  [ Group([ f2, f3, f4, f5, f6 ]), Group([ f6 ]) ], 
  [ Group([ f2, f3, f4, f5, f6 ]), Group([  ]) ] ]

# wedderga/doc/SSP.xml:75-87

gap> G:=SymmetricGroup(3);; K:=Group([(1,2,3)]);; H:=Group( () );;
gap> IsStrongShodaPair( G, K, H );
true
gap> IsStrongShodaPair( G, G, H );
false
gap> IsStrongShodaPair( G, K, K );
false
gap> IsStrongShodaPair( G, G, K );
true

# wedderga/doc/SSP.xml:107-117

gap> G:=AlternatingGroup(5);;
gap> K:=AlternatingGroup(4);;
gap> H := Group( (1,2)(3,4), (1,3)(2,4) );;
gap> IsStrongShodaPair( G, K, H );
false
gap> IsShodaPair( G, K, H );
true

# wedderga/doc/SSP.xml:134-150

gap> S4:=SymmetricGroup(4);;
gap> IsStronglyMonomial(S4);
true
gap> G:=SmallGroup(24,3);;
gap> IsStronglyMonomial(G);
false
gap> IsMonomial(G);
false
gap> G:=SmallGroup(1000,86);;
gap> IsMonomial(G);
true
gap> IsStronglyMonomial(G);
false

gap> STOP_TEST("wedderga03.tst", 1 );
