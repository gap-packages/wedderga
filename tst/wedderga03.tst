# wedderga, chapter 3
gap> START_TEST( "wedderga03.tst");

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/SSP.xml", 24, 45 ]

gap> StrongShodaPairs( SymmetricGroup(4) );
[ [ Sym( [ 1 .. 4 ] ), Group([ (1,4)(2,3), (1,3)(2,4), (2,4,3), (3,4) ]) ],
  [ Sym( [ 1 .. 4 ] ), Group([ (1,4)(2,3), (1,3)(2,4), (2,4,3) ]) ],
  [ Group([ (3,4), (1,3,2,4) ]), Group([ (1,3,2,4), (1,2)(3,4) ]) ],
  [ Group([ (1,3,2,4), (3,4) ]), Group([ (3,4), (1,2)(3,4) ]) ],
  [ Group([ (2,4,3), (1,4)(2,3) ]), Group([ (1,4)(2,3), (1,3)(2,4) ]) ] ]
gap> StrongShodaPairs( DihedralGroup(64) );
[ [ <pc group of size 64 with 6 generators>,
      Group([ f6, f5, f4, f3, f1, f2 ]) ],
  [ <pc group of size 64 with 6 generators>, Group([ f6, f5, f4, f3, f1*f2 ])
     ],
  [ <pc group of size 64 with 6 generators>, Group([ f6, f5, f4, f3, f2 ]) ],
  [ <pc group of size 64 with 6 generators>, Group([ f6, f5, f4, f3, f1 ]) ],
  [ Group([ f1*f2, f4*f5*f6, f5*f6, f6, f3, f3 ]),
      Group([ f6, f5, f4, f1*f2 ]) ],
  [ Group([ f6, f5, f2, f3, f4 ]), Group([ f6, f5 ]) ],
  [ Group([ f6, f2, f3, f4, f5 ]), Group([ f6 ]) ],
  [ Group([ f2, f3, f4, f5, f6 ]), Group([  ]) ] ]

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/SSP.xml", 68, 80 ]

gap> G:=SymmetricGroup(3);; K:=Group([(1,2,3)]);; H:=Group( () );;
gap> IsStrongShodaPair( G, K, H );
true
gap> IsStrongShodaPair( G, G, H );
false
gap> IsStrongShodaPair( G, K, K );
false
gap> IsStrongShodaPair( G, G, K );
true

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/SSP.xml", 100, 110 ]

gap> G:=AlternatingGroup(5);;
gap> K:=AlternatingGroup(4);;
gap> H := Group( (1,2)(3,4), (1,3)(2,4) );;
gap> IsStrongShodaPair( G, K, H );
false
gap> IsShodaPair( G, K, H );
true

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/SSP.xml", 127, 143 ]

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
