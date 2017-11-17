# wedderga, chapter 1
gap> START_TEST( "wedderga01.tst");

# wedderga/doc/intro.xml:135-158

gap> QG := GroupRing( Rationals, SymmetricGroup(4) );
<algebra-with-one over Rationals, with 2 generators>
gap> WedderburnDecomposition(QG);
[ Rationals, Rationals, <crossed product with center Rationals over CF(
    3) of a group of size 2>, ( Rationals^[ 3, 3 ] ), ( Rationals^[ 3, 3 ] ) ]
gap> FG := GroupRing( CF(5), SymmetricGroup(4) );
<algebra-with-one over CF(5), with 2 generators>
gap> WedderburnDecomposition( FG );
[ CF(5), CF(5), <crossed product with center CF(5) over AsField( CF(5), CF(
    15) ) of a group of size 2>, ( CF(5)^[ 3, 3 ] ), ( CF(5)^[ 3, 3 ] ) ]
gap> FG := GroupRing( GF(5), SymmetricGroup(4) ); 
<algebra-with-one over GF(5), with 2 generators>
gap> WedderburnDecomposition( FG );
[ ( GF(5)^[ 1, 1 ] ), ( GF(5)^[ 1, 1 ] ), ( GF(5)^[ 2, 2 ] ), 
  ( GF(5)^[ 3, 3 ] ), ( GF(5)^[ 3, 3 ] ) ]
gap> FG := GroupRing( GF(5), SmallGroup(24,3) );
<algebra-with-one over GF(5), with 4 generators>
gap> WedderburnDecomposition( FG );
[ ( GF(5)^[ 1, 1 ] ), ( GF(5^2)^[ 1, 1 ] ), ( GF(5)^[ 2, 2 ] ), 
  ( GF(5^2)^[ 2, 2 ] ), ( GF(5)^[ 3, 3 ] ) ]

gap> STOP_TEST("wedderga01.tst", 1 );
