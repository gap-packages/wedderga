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

# wedderga/doc/SSP.xml:75-89

gap> G:=SymmetricGroup(4);; K:=Group( (1,3,2,4), (3,4) );;
gap> H1:=Group( (2,4,3), (1,4)(2,3), (1,3)(2,4) );;
gap> H2:=Group( (3,4), (1,2)(3,4) );;
gap> IsExtremelyStrongShodaPair( G, G, H1 );
true
gap> IsExtremelyStrongShodaPair( G, K, H2 );
false
gap> IsExtremelyStrongShodaPair( G, G, H2 );
false
gap> IsExtremelyStrongShodaPair( G, G, K );
false

# wedderga/doc/SSP.xml:107-121

gap> G:=SymmetricGroup(4);; K:=Group( (1,3,2,4), (3,4) );;
gap> H1:=Group( (2,4,3), (1,4)(2,3), (1,3)(2,4) );;
gap> H2:=Group( (3,4), (1,2)(3,4) );;
gap> IsStrongShodaPair( G, G, H1 );
true
gap> IsExtremelyStrongShodaPair( G, K, H2 );
false
gap> IsStrongShodaPair( G, K, H2 );
true
gap> IsStrongShodaPair( G, G, K );
false

# wedderga/doc/SSP.xml:141-151

gap> G:=AlternatingGroup(5);;
gap> K:=AlternatingGroup(4);;
gap> H := Group( (1,2)(3,4), (1,3)(2,4) );;
gap> IsStrongShodaPair( G, K, H );
false
gap> IsShodaPair( G, K, H );
true

# wedderga/doc/SSP.xml:168-184

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
