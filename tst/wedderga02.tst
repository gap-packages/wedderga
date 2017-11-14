# wedderga, chapter 2
gap> START_TEST( "wedderga02.tst");

# wedderga/doc/decomp.xml:31-45

gap> WedderburnDecomposition( GroupRing( GF(5), DihedralGroup(16) ) );
[ ( GF(5)^[ 1, 1 ] ), ( GF(5)^[ 1, 1 ] ), ( GF(5)^[ 1, 1 ] ),
  ( GF(5)^[ 1, 1 ] ), ( GF(5)^[ 2, 2 ] ), ( GF(5^2)^[ 2, 2 ] ) ]
gap> WedderburnDecomposition( GroupRing( Rationals, DihedralGroup(16) ) );
[ Rationals, Rationals, Rationals, Rationals, ( Rationals^[ 2, 2 ] ),
  <crossed product with center NF(8,[ 1, 7 ]) over AsField( NF(8,
    [ 1, 7 ]), CF(8) ) of a group of size 2> ]
gap> WedderburnDecomposition( GroupRing( CF(5), DihedralGroup(16) ) );
[ CF(5), CF(5), CF(5), CF(5), ( CF(5)^[ 2, 2 ] ),
  <crossed product with center NF(40,[ 1, 31 ]) over AsField( NF(40,
    [ 1, 31 ]), CF(40) ) of a group of size 2> ]

# wedderga/doc/decomp.xml:77-94

gap> WedderburnDecomposition( GroupRing( Rationals, SmallGroup(48,15) ) );
[ Rationals, Rationals, Rationals, Rationals, ( Rationals^[ 2, 2 ] ),
  <crossed product with center Rationals over CF(3) of a group of size 2>,
  ( CF(3)^[ 2, 2 ] ), <crossed product with center Rationals over CF(
    3) of a group of size 2>, <crossed product with center NF(8,
    [ 1, 7 ]) over AsField( NF(8,[ 1, 7 ]), CF(8) ) of a group of size 2>,
  <crossed product with center Rationals over CF(12) of a group of size 4> ]
gap> WedderburnDecomposition( GroupRing( CF(3), SmallGroup(48,15) ) );
[ CF(3), CF(3), CF(3), CF(3), ( CF(3)^[ 2, 2 ] ), ( CF(3)^[ 2, 2 ] ),
  ( CF(3)^[ 2, 2 ] ), ( CF(3)^[ 2, 2 ] ), ( CF(3)^[ 2, 2 ] ),
  <crossed product with center NF(24,[ 1, 7 ]) over AsField( NF(24,
    [ 1, 7 ]), CF(24) ) of a group of size 2>,
  ( <crossed product with center CF(3) over AsField( CF(3), CF(
    12) ) of a group of size 2>^[ 2, 2 ] ) ]

# wedderga/doc/decomp.xml:106-122

gap> QG:=GroupRing(Rationals,SmallGroup(240,89));
<algebra-with-one over Rationals, with 2 generators>
gap> WedderburnDecomposition(QG);
Wedderga: Warning!!!
Some of the Wedderburn components displayed are FRACTIONAL MATRIX ALGEBRAS!!!

[ Rationals, Rationals, <crossed product with center Rationals over CF(
    5) of a group of size 4>, ( Rationals^[ 4, 4 ] ), ( Rationals^[ 4, 4 ] ),
  ( Rationals^[ 5, 5 ] ), ( Rationals^[ 5, 5 ] ), ( Rationals^[ 6, 6 ] ),
  <crossed product with center NF(12,[ 1, 11 ]) over AsField( NF(12,
    [ 1, 11 ]), NF(60,[ 1, 11 ]) ) of a group of size 4>,
  [ 3/2, <crossed product with center NF(8,[ 1, 7 ]) over AsField( NF(8,
        [ 1, 7 ]), NF(40,[ 1, 31 ]) ) of a group of size 4> ] ]  

# wedderga/doc/decomp.xml:191-200

gap> WedderburnDecompositionInfo( GroupRing( Rationals, DihedralGroup(16) ) );
[ [ 1, Rationals ], [ 1, Rationals ], [ 1, Rationals ], [ 1, Rationals ],
  [ 2, Rationals ], [ 1, NF(8,[ 1, 7 ]), 8, [ 2, 7, 0 ] ] ]
gap> WedderburnDecompositionInfo( GroupRing( CF(5), DihedralGroup(16) ) );
[ [ 1, CF(5) ], [ 1, CF(5) ], [ 1, CF(5) ], [ 1, CF(5) ], [ 2, CF(5) ],
  [ 1, NF(40,[ 1, 31 ]), 8, [ 2, 7, 0 ] ] ]

# wedderga/doc/decomp.xml:218-237

gap> F:=FreeGroup("a","b");;a:=F.1;;b:=F.2;;rel:=[a^8,a^4*b^2,b^-1*a*b*a];;
gap> Q16:=F/rel;; QQ16:=GroupRing( Rationals, Q16 );;
gap> QS4:=GroupRing( Rationals, SymmetricGroup(4) );;
gap> WedderburnDecomposition(QQ16);
[ Rationals, Rationals, Rationals, Rationals, ( Rationals^[ 2, 2 ] ),
  <crossed product with center NF(8,[ 1, 7 ]) over AsField( NF(8,
    [ 1, 7 ]), CF(8) ) of a group of size 2> ]
gap> WedderburnDecomposition( QS4 );
[ Rationals, Rationals, ( Rationals^[ 3, 3 ] ), ( Rationals^[ 3, 3 ] ),
  <crossed product with center Rationals over CF(3) of a group of size 2> ]
gap> WedderburnDecompositionInfo(QQ16);
[ [ 1, Rationals ], [ 1, Rationals ], [ 1, Rationals ], [ 1, Rationals ], 
  [ 2, Rationals ], [ 1, NF(8,[ 1, 7 ]), 8, [ 2, 7, 4 ] ] ]
gap> WedderburnDecompositionInfo(QS4);  
[ [ 1, Rationals ], [ 1, Rationals ], [ 3, Rationals ], [ 3, Rationals ], 
  [ 1, Rationals, 3, [ 2, 2, 0 ] ] ]

# wedderga/doc/decomp.xml:291-304

gap> WedderburnDecompositionInfo( GroupRing( Rationals, SmallGroup(48,15) ) );
[ [ 1, Rationals ], [ 1, Rationals ], [ 1, Rationals ], [ 1, Rationals ], 
  [ 2, Rationals ], [ 1, Rationals, 3, [ 2, 2, 0 ] ], [ 2, CF(3) ], 
  [ 1, Rationals, 6, [ 2, 5, 0 ] ], [ 1, NF(8,[ 1, 7 ]), 8, [ 2, 7, 0 ] ], 
  [ 1, Rationals, 12, [ [ 2, 5, 3 ], [ 2, 7, 0 ] ], [ [ 3 ] ] ] ]
gap> WedderburnDecompositionInfo( GroupRing( CF(3), SmallGroup(48,15) ) );
[ [ 1, CF(3) ], [ 1, CF(3) ], [ 1, CF(3) ], [ 1, CF(3) ], [ 2, CF(3) ],
  [ 2, CF(3), 3, [ 1, 1, 0 ] ], [ 2, CF(3) ], [ 2, CF(3) ],
  [ 2, CF(3), 6, [ 1, 1, 0 ] ], [ 1, NF(24,[ 1, 7 ]), 8, [ 2, 7, 0 ] ],
  [ 2, CF(3), 12, [ 2, 7, 0 ] ] ]

# wedderga/doc/decomp.xml:317-330

gap> QG:=GroupRing(Rationals,SmallGroup(240,89));
<algebra-with-one over Rationals, with 2 generators>
gap> WedderburnDecompositionInfo(QG);
Wedderga: Warning!!! 
Some of the Wedderburn components displayed are FRACTIONAL MATRIX ALGEBRAS!!!

[ [ 1, Rationals ], [ 1, Rationals ], [ 1, Rationals, 10, [ 4, 3, 5 ] ],
  [ 4, Rationals ], [ 4, Rationals ], [ 5, Rationals ], [ 5, Rationals ],
  [ 6, Rationals ], [ 1, NF(12,[ 1, 11 ]), 10, [ 4, 3, 5 ] ],
  [ 3/2, NF(8,[ 1, 7 ]), 10, [ 4, 3, 5 ] ] ]

# wedderga/doc/decomp.xml:387-408

gap> A5 := AlternatingGroup(5);
Alt( [ 1 .. 5 ] )
gap> SimpleAlgebraByCharacter( GroupRing( Rationals , A5 ) , Irr( A5 ) [3] );
( NF(5,[ 1, 4 ])^[ 3, 3 ] )
gap> SimpleAlgebraByCharacter( GroupRing( GF(7) , A5 ) , Irr( A5 ) [3] );
( GF(7^2)^[ 3, 3 ] )
gap> G:=SmallGroup(128,100);               
<pc group of size 128 with 7 generators>
gap> chi4:=Filtered(Irr(G),x->Degree(x)=4);;
gap> List(chi4,x->SimpleAlgebraByCharacter(GroupRing(Rationals,G),x));
[ ( <crossed product with center NF(8,[ 1, 3 ]) over AsField( NF(8,
    [ 1, 3 ]), CF(8) ) of a group of size 2>^[ 2, 2 ] ), 
  ( <crossed product with center NF(8,[ 1, 3 ]) over AsField( NF(8,
    [ 1, 3 ]), CF(8) ) of a group of size 2>^[ 2, 2 ] ), 
  ( <crossed product with center NF(8,[ 1, 3 ]) over AsField( NF(8,
    [ 1, 3 ]), CF(8) ) of a group of size 2>^[ 2, 2 ] ), 
  ( <crossed product with center NF(8,[ 1, 3 ]) over AsField( NF(8,
    [ 1, 3 ]), CF(8) ) of a group of size 2>^[ 2, 2 ] ) ]

# wedderga/doc/decomp.xml:438-450

gap> G:=SmallGroup(144,11);
<pc group of size 144 with 6 generators>
gap> QG:=GroupRing(Rationals,G);
<algebra-with-one over Rationals, with 6 generators>
gap> SimpleAlgebraByCharacter( QG , Irr(G)[40] );
<crossed product with center NF(36,[ 1, 17 ]) over AsField( NF(36,
[ 1, 17 ]), CF(36) ) of a group of size 2>
gap> SimpleAlgebraByCharacterInfo( QG , Irr(G)[48] );
[ 1, NF(9,[ 1, 8 ]), 18, [ 2, 17, 9 ] ]

# wedderga/doc/decomp.xml:519-533

gap> F:=FreeGroup("a","b");; a:=F.1;; b:=F.2;;
gap> G:=F/[ a^16, b^2*a^8, b^-1*a*b*a^9 ];; a:=G.1;; b:=G.2;;
gap> K:=Subgroup(G,[a]);; H:=Subgroup(G,[]);;
gap> QG:=GroupRing( Rationals, G );;
gap> FG:=GroupRing( GF(7), G );;
gap> SimpleAlgebraByStrongSP( QG, K, H );
<crossed product over CF(16) of a group of size 2>
gap> SimpleAlgebraByStrongSP( FG, K, H, [1,7] );
( GF(7)^[ 2, 2 ] )
gap> SimpleAlgebraByStrongSP( FG, K, H, 1 );
( GF(7)^[ 2, 2 ] )

# wedderga/doc/decomp.xml:593-609

gap> F:=FreeGroup("a","b");; a:=F.1;; b:=F.2;;
gap> G:=F/[ a^16, b^2*a^8, b^-1*a*b*a^9 ];; a:=G.1;; b:=G.2;;
gap> K:=Subgroup(G,[a]);; H:=Subgroup(G,[]);; 
gap> QG:=GroupRing( Rationals, G );;
gap> FG:=GroupRing( GF(7), G );;
gap> SimpleAlgebraByStrongSP( QG, K, H );
<crossed product over CF(16) of a group of size 2>
gap> SimpleAlgebraByStrongSPInfo( QG, K, H );
[ 1, NF(16,[ 1, 7 ]), 16, [ [ 2, 7, 8 ] ], [  ] ]
gap> SimpleAlgebraByStrongSPInfo( FG, K, H, [1,7] );
[ 2, 7 ]
gap> SimpleAlgebraByStrongSPInfo( FG, K, H, 1 );
[ 2, 7 ]

gap> STOP_TEST("wedderga02.tst", 1 );
