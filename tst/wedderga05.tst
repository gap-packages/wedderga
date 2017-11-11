# wedderga, chapter 5
gap> START_TEST( "wedderga05.tst");

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 91, 122 ]

gap> R := GaussianRationals;
GaussianRationals
gap> G := Units( ZmodnZ(4) );
<group of size 2 with 1 generators>
gap> act := function(RG,g)
> return ANFAutomorphism( LeftActingDomain(RG), Int(g) );
> end;
function( RG, g ) ... end
gap> twist1 := function( RG, g, h )
> if IsOne(g) or IsOne(h) then
>    return One(LeftActingDomain(RG));
> else
>    return -One(LeftActingDomain(RG));
> fi;
> end;
function( RG, g, h ) ... end
gap> RG := CrossedProduct( R, G, act, twist1 );
<crossed product over GaussianRationals of a group of size 2>
gap> i := E(4) * One(G)^Embedding(G,RG); 
(ZmodnZObj( 1, 4 ))*(E(4))
gap> j := ZmodnZObj(3,4)^Embedding(G,RG); 
(ZmodnZObj( 3, 4 ))*(1)
gap> i^2;
(ZmodnZObj( 1, 4 ))*(-1)
gap> j^2;
(ZmodnZObj( 1, 4 ))*(-1)
gap> i*j+j*i;  
<zero> of ...

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 131, 154 ]

gap> twist2:=function(RG,g,h)
> if IsOne(g) or IsOne(h) then
>     return One(LeftActingDomain( RG ));
> else
>     return -3*One(LeftActingDomain( RG ));
> fi;
> end;
function( RG, g, h ) ... end
gap> RG := CrossedProduct( R, G, act, twist2 );  
<crossed product over GaussianRationals of a group of size 2>
gap> i := E(4) * One(G)^Embedding(G,RG); 
(ZmodnZObj( 1, 4 ))*(E(4))
gap> j := ZmodnZObj(3,4)^Embedding(G,RG);  
(ZmodnZObj( 3, 4 ))*(1)
gap> i^2;                           
(ZmodnZObj( 1, 4 ))*(-1)
gap> j^2;                                
(ZmodnZObj( 1, 4 ))*(-3)
gap> i*j+j*i;                       
<zero> of ...

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 160, 192 ]

gap> C2 := CyclicGroup(2);
<pc group of size 2 with 1 generators>
gap> G := DirectProduct(C2,C2);
<pc group of size 4 with 2 generators>
gap> act := function(RG,a)
>     return IdentityMapping( LeftActingDomain(RG));
> end;
function( RG, a ) ... end
gap> twist := function( RG, g , h )
> local one,g1,g2,h1,h2,G;
> G := UnderlyingMagma( RG );
> one := One( C2 );
> g1 := Image( Projection(G,1), g );
> g2 := Image( Projection(G,2), g );
> h1 := Image( Projection(G,1), h );
> h2 := Image( Projection(G,2), h );
> if g = One( G ) or h = One( G ) then return 1;
>   elif IsOne(g1) and not IsOne(g2) and not IsOne(h1) and not IsOne(h2)
>     then return 1;
>   elif not IsOne(g1) and IsOne(g2) and IsOne(h1) and not IsOne(h2)
>     then return 1;
>   elif not IsOne(g1) and not IsOne(g2) and not IsOne(h1) and IsOne(h2)
>     then return 1;
>   else return -1;
> fi;
> end;
function( RG, g, h ) ... end
gap> HQ := CrossedProduct( Rationals, G, act, twist );
<crossed product over Rationals of a group of size 4>

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 196, 211 ]

gap> HZ := CrossedProduct( Integers, G, act, twist );
<crossed product over Integers of a group of size 4>
gap> i := GeneratorsOfGroup(G)[1]^Embedding(G,HZ); 
(f1)*(1)
gap> j := GeneratorsOfGroup(G)[2]^Embedding(G,HZ);
(f2)*(1)
gap> i^2;
(<identity> of ...)*(-1)
gap> j^2; 
(<identity> of ...)*(-1)
gap> i*j+j*i;                                      
<zero> of ...

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 224, 240 ]

gap> LeftActingDomain(HZ);
Integers
gap> G:=UnderlyingMagma(HZ);
<pc group of size 4 with 2 generators>
gap> ac := ActionForCrossedProduct(HZ);
function( RG, a ) ... end
gap> List( G , x -> ac( HZ, x ) );
[ IdentityMapping( Integers ), IdentityMapping( Integers ),
  IdentityMapping( Integers ), IdentityMapping( Integers ) ]
gap> tw := TwistingForCrossedProduct( HZ );
function( RG, g, h ) ... end
gap> List( G, x -> List( G , y -> tw( HZ, x, y ) ) );
[ [ 1, 1, 1, 1 ], [ 1, -1, -1, 1 ], [ 1, 1, -1, -1 ], [ 1, -1, 1, -1 ] ]  

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 246, 284 ]

gap> G := SmallGroup(32,50);
<pc group of size 32 with 5 generators>
gap> A := SimpleAlgebraByCharacter( GroupRing(Rationals,G), Irr(G)[17]) ;
( <crossed product with center Rationals over GaussianRationals of a group of \
size 2>^[ 2, 2 ] )
gap> SimpleAlgebraByCharacterInfo( GroupRing(Rationals,G), Irr(G)[17]) ;
[ 2, Rationals, 4, [ 2, 3, 2 ] ]
gap> B := LeftActingDomain(A);
<crossed product with center Rationals over GaussianRationals of a group of si\
ze 2>
gap> L := LeftActingDomain(B);
GaussianRationals
gap> H := UnderlyingMagma( B );
<group of size 2 with 2 generators>
gap> Elements(H);
[ ZmodnZObj( 1, 4 ), ZmodnZObj( 3, 4 ) ]
gap> i := E(4) * One(H)^Embedding(H,B);
(ZmodnZObj( 1, 4 ))*(E(4))
gap> j := ZmodnZObj(3,4)^Embedding(H,B);
(ZmodnZObj( 3, 4 ))*(1)
gap> i^2;
(ZmodnZObj( 1, 4 ))*(-1)
gap> j^2;
(ZmodnZObj( 1, 4 ))*(-1)
gap> i*j+j*i;
<zero> of ...
gap> ac := ActionForCrossedProduct( B );
function( RG, a ) ... end
gap> tw := TwistingForCrossedProduct( B );
function( RG, a, b ) ... end
gap> List( H , x -> ac( B, x ) );
[ IdentityMapping( GaussianRationals ), ANFAutomorphism( GaussianRationals,
    3 ) ]
gap> List( H , x -> List( H , y -> tw( B, x, y ) ) );
[ [ 1, 1 ], [ 1, -1 ] ]

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 287, 342 ]

gap> QG:=GroupRing( Rationals, SmallGroup(24,3) );;
gap> WedderburnDecomposition(QG);
[ Rationals, CF(3), ( Rationals^[ 3, 3 ] ),
  <crossed product with center Rationals over GaussianRationals of a group of \
size 2>, <crossed product with center CF(3) over AsField( CF(3), CF(
    12) ) of a group of size 2> ]
gap> R:=WedderburnDecomposition( QG )[4];
<crossed product with center Rationals over GaussianRationals of a group of si\
ze 2>
gap> IsCrossedProduct(R);
true
gap> IsAlgebra(R);
true
gap> IsRing(R);       
true
gap> LeftActingDomain( R );
GaussianRationals
gap> AsList( UnderlyingMagma( R ) );
[ ZmodnZObj( 1, 4 ), ZmodnZObj( 3, 4 ) ]
gap> Print( ActionForCrossedProduct( R ) ); Print("\n");
function ( RG, a )
    local  cond, redu;
    cond := OperationRecord( RG ).cond;
    redu := OperationRecord( RG ).redu;
    return
     ANFAutomorphism( CF( cond ), Int( PreImagesRepresentative( redu, a ) ) );
end
gap> Print( TwistingForCrossedProduct( R ) ); Print("\n");                     
function ( RG, a, b )
    local  orderroot, cocycle;
    orderroot := OperationRecord( RG ).orderroot;
    cocycle := OperationRecord( RG ).cocycle;
    return E( orderroot ) ^ Int( cocycle( a, b ) );
end
gap> IsAssociative(R);
true
gap> IsFinite(R);           
false
gap> IsFiniteDimensional(R);
true
gap> AsList(Basis(R));
[ (ZmodnZObj( 1, 4 ))*(1), (ZmodnZObj( 3, 4 ))*(1) ] 
gap> GeneratorsOfLeftOperatorRingWithOne(R);
[ (ZmodnZObj( 1, 4 ))*(1), (ZmodnZObj( 3, 4 ))*(1) ]
gap> One(R);
(ZmodnZObj( 1, 4 ))*(1)
gap> Zero(R);
<zero> of ...
gap> Characteristic(R);
0
gap> CenterOfCrossedProduct(R);
Rationals

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 348, 402 ]

gap> Quat := function(R,a,b)
> local G,act,twist;
> if not(a in R and b in R and a <> Zero(R) and b <> Zero(R) ) then
> Error("<a>  and <b> must be non zero elements of <R>!!!");
> fi;
> G := SmallGroup(4,2);
> act := function(RG,a)
>     return IdentityMapping( LeftActingDomain(RG));
> end;
> twist := function( RG, g , h )
> local one,g1,g2;
> one := One(G);
> g1 := G.1;
> g2 := G.2;
> if   g = one or h = one then
>   return One(R);
> elif g = g1 then
>   if h = g2 then
>     return One(R);
>   else
>     return a;
>   fi;
> elif g = g2 then
>   if h = g1 then
>     return -One(R);
>   elif h=g2 then
>     return b;
>   else
>     return -b;
>   fi;
> else
>   if h = g1 then
>     return -b;
>   elif h=g2 then
>     return b;
>   else
>     return -a*b;
>   fi;
> fi;
> end;
> return CrossedProduct(R,G,act,twist);
> end;
function( R, a, b ) ... end
gap> HQ := Quat(Rationals,2,3);
<crossed product over Rationals of a group of size 4>
gap> G := UnderlyingMagma(HQ);
<pc group of size 4 with 2 generators>
gap> tw := TwistingForCrossedProduct( HQ );
function( RG, g, h ) ... end
gap> List( G, x -> List( G, y -> tw( HQ, x, y ) ) );
[ [ 1, 1, 1, 1 ], [ 1, 3, -1, -3 ], [ 1, 1, 2, 2 ], [ 1, 3, -3, -6 ] ]

# [ "/Users/alexk/GITREPS/gap-stable/pkg/wedderga/doc/crossed.xml", 448, 480 ]

gap> QG := GroupRing( Rationals, SmallGroup(24,3) );
<algebra-with-one over Rationals, with 4 generators>
gap> R := WedderburnDecomposition( QG )[4];
<crossed product with center Rationals over GaussianRationals of a group of si\
ze 2>
gap> H := UnderlyingMagma( R );;
gap> fam := ElementsFamily( FamilyObj( R ) );;
gap> g := ElementOfCrossedProduct( fam, 0, [ 1, E(4) ], AsList(H) );
(ZmodnZObj( 1, 4 ))*(1)+(ZmodnZObj( 3, 4 ))*(E(4))
gap> CoefficientsAndMagmaElements( g );    
[ ZmodnZObj( 1, 4 ), 1, ZmodnZObj( 3, 4 ), E(4) ]
gap> t := List( H, x -> x^Embedding( H, R ) );
[ (ZmodnZObj( 1, 4 ))*(1), (ZmodnZObj( 3, 4 ))*(1) ]
gap> t[1] + t[2]*E(4);  
(ZmodnZObj( 1, 4 ))*(1)+(ZmodnZObj( 3, 4 ))*(E(4))
gap> g = t[1] + E(4)*t[2];
false
gap> g = t[1] + t[2]*E(4);
true
gap> h := ElementOfCrossedProduct( fam, 0, [ E(4), 1 ], AsList(H) );     
(ZmodnZObj( 1, 4 ))*(E(4))+(ZmodnZObj( 3, 4 ))*(1)
gap> g+h;
(ZmodnZObj( 1, 4 ))*(1+E(4))+(ZmodnZObj( 3, 4 ))*(1+E(4))
gap> g*E(4);
(ZmodnZObj( 1, 4 ))*(E(4))+(ZmodnZObj( 3, 4 ))*(-1)
gap> E(4)*g;     
(ZmodnZObj( 1, 4 ))*(E(4))+(ZmodnZObj( 3, 4 ))*(1)
gap> g*h;
(ZmodnZObj( 1, 4 ))*(2*E(4))

gap> STOP_TEST("wedderga05.tst", 1 );
