#############################################################################
##
#W  auxiliar.gi           The Wedderga package            Osnel Broche Cristo
#W                                                         Olexandr Konovalov
#W                                                            Aurora Olivieri
#W                                                           Gabriela Olteanu
#W                                                              Ángel del Río
#W                                                          Inneke Van Gelder
##
#############################################################################


#############################################################################
##
#P IsSemisimpleRationalGroupAlgebra( FG )
##  
## The function checks whether a group ring is a rational group algebra 
## of a finite group
##
InstallImmediateMethod( IsSemisimpleRationalGroupAlgebra,
    IsGroupRing,
    0,
R -> IsRationals(LeftActingDomain(R)) and IsFinite(UnderlyingMagma(R))); 


#############################################################################
##
#P IsSemisimpleZeroCharacteristicGroupAlgebra( FG )
##  
## The function checks whether a group ring is a group algebra 
## of a finite group over the field of characteristic zero
##
InstallImmediateMethod( IsSemisimpleZeroCharacteristicGroupAlgebra,
    IsGroupRing, 
    0,
R -> Characteristic(LeftActingDomain(R))=0 and IsFinite(UnderlyingMagma(R))); 


#############################################################################
##
#P IsCFGroupAlgebra( FG )
##  
## The function checks whether a group ring is a group algebra of a finite
## group over a cyclotomic field
##
InstallImmediateMethod( IsCFGroupAlgebra,
    IsGroupRing, 
    0,
R ->  IsField(LeftActingDomain(R)) and IsCyclotomicField(LeftActingDomain(R)) and IsFinite(UnderlyingMagma(R))); 

#############################################################################
##
#P IsSemisimpleANFGroupAlgebra( FG )
##  
## The function checks whether a group ring is a group algebra of a finite
## group over an abelian number field
##
InstallImmediateMethod( IsSemisimpleANFGroupAlgebra,
    IsGroupRing, 
    0,
R -> IsField(LeftActingDomain(R)) and IsAbelianNumberField(LeftActingDomain(R)) and IsFinite(UnderlyingMagma(R))); 


#############################################################################
##
#P IsSemisimpleFiniteGroupAlgebra( FG )
##  
## The function checks whether a group ring is a semisimple group algebra
## of a finite group over a finite field
##
InstallImmediateMethod( IsSemisimpleFiniteGroupAlgebra,
    IsGroupRing, 
    0,
function( FG )
	local   F,      # Field
    	    G;      # Group

	F := LeftActingDomain( FG );
	G := UnderlyingMagma( FG );
       
	return IsField( F ) and IsFinite( F ) and IsFinite( G ) and 
    	   Gcd( Size( F ), Size( G ) ) = 1;
end); 


#############################################################################
##
#M IsCompleteSetOfPCIs( R, ListPCIs )
##
## The function checks if the sum of given idempotents of a ring R is the
## identity element of R. It is supposed to be used to check if a given 
## list of PCIs of R is complete.
##
InstallMethod( IsCompleteSetOfPCIs,
    "for list of idempotents", 
    true, 
    [ IsRing, IsList ], 
    0,
function( R, ListPCIs )
    local x;
    if not ForAll( ListPCIs, x -> x in R ) then
        Error("Wedderga: An element of <ListPCIs> does not belong to <R>!!!\n");
    else
        return Sum( ListPCIs ) = One( R ) and ForAll( ListPCIs, x -> x=x^2 );
    fi;
end);

#############################################################################
##
#M IsCompleteSetOfOrthogonalIdempotents( R, List )
##
## The function checks if List is a complete set of orthogonal central 
## idempotents of a ring R.
##
InstallMethod( IsCompleteSetOfOrthogonalIdempotents,
    "for list of idempotents", 
    true, 
    [ IsRing, IsList ], 
    0,
function( R, ListPCIs )
    if not ForAll( ListPCIs, x -> x in R ) then
        Error("Wedderga: An element of <ListPCIs> does not belong to <R> !!!\n");
    elif ForAny( ListPCIs, IsZero ) then
        Error("Wedderga: Zero element in  <ListPCIs> !!!\n");
    else
        return Sum( ListPCIs ) = One( R ) and 
                    ForAll( [1..Length(ListPCIs)], i -> 
                        ListPCIs[i]= ListPCIs[i]^2 and 
                        ForAll ( [(i+1)..Length(ListPCIs)], j -> 
                            ListPCIs[i]* ListPCIs[j] = Zero(R) ) ) ;
    fi;
end);




#############################################################################
##
#F IsStrongShodaPair( G, K, H )
##
## The function IsStrongShodaPair checks if (H,K) is a strong Shoda pair of G
##
InstallMethod( IsStrongShodaPair,
    "for a group and two subgroups", 
    true,
    [ IsGroup, IsGroup, IsGroup ], 
    0,
function( G, K, H )
local   QG,
        NH,
        Eps,
        NdK,
        eGKH1,
        RTNH,
        i,
        g,
        RTNdK,
        nRTNdK,
        Epi,
        NHH,
        KH,
		zero;

# First verifies if H, K are subgroups of G and K is a normal subgroup of K
if not ( IsSubgroup( G, K ) and IsSubgroup( K, H ) ) then
    Error("Wedderga: <G> must contain <K> and <K> must contain <H> !!!\n");
fi;

if not IsNormal( K, H ) then
    Info(InfoWedderga, 2, "Wedderga: IsSSP: <H> is not normal in <K>");
    return false;
fi;

# Now, if K/H is the maximal abelian subgroup in N/H,
# where N is the normalizer of H in G

NH:=Normalizer(G,H);

if not(IsNormal( NH, K ) ) then
    Info(InfoWedderga, 2, "Wedderga: IsSSP: <K> is not normal in N_<G>(<H>)");
    return false;
fi;

Epi:=NaturalHomomorphismByNormalSubgroup( NH, H ) ;
NHH:=Image( Epi, NH ); #It is isomorphic to the factor group NH/H.
KH:=Image( Epi, K ); #It is isomorphic to the factor group K/H.

if not(IsCyclic(KH)) then
    Info(InfoWedderga, 2, "Wedderga: IsSSP: <K>/<H> is not cyclic");
    return false;
fi;

if Centralizer( NHH, KH ) <> KH then
    Info(InfoWedderga, 2, "Wedderga: IsSSP: <K>/<H> is not maximal ",
                     "abelian in N_<G>(<H>)");
    return false;
fi;

#Now (SSS3)
QG := GroupRing( Rationals, G );
zero := Zero( QG );
Eps := IdempotentBySubgroups( QG, K, H );
NdK := Normalizer( G, K );

if NdK<>G then
    RTNH := RightTransversal( NdK, NH );
    eGKH1 := Sum( List( RTNH, g -> Eps^g ) ); 
    RTNdK := RightTransversal( G, NdK );
    nRTNdK:=Length(RTNdK);
    for i in [ 2 .. nRTNdK ] do
        g:=RTNdK[i];
        if eGKH1*eGKH1^g <> zero then
            Info(InfoWedderga, 2, 
                 "Wedderga: IsSSP: The conjugates of epsilon are not orthogonal");
            return  false;
        fi;
    od;
fi;

return true;

end);


#############################################################################
##
#M Centralizer( G, a )
##
## The function Centralizer computes the centralizer of an
## element of a group ring in a subgroup of the underlying group
##
InstallMethod( Centralizer,
    "Wedderga: for a subgroup of an underlying group and a group ring element",
    function( F1, F2 )    
      return IsBound( F2!.familyMagma ) and
             IsIdenticalObj( F1, F2!.familyMagma);
    end,  
    [ IsGroup, IsElementOfFreeMagmaRing ], 
    0,
function( G, a ) 
local   C,
        leftover,
        cosrep,
        g, 
        h;

C := TrivialSubgroup( G );
cosrep := [ One(G) ];
leftover := Difference( G, C );
while leftover <> [] do
    g := leftover[1];
    if a^g = a then 
        C := Subgroup( G, Union( GeneratorsOfGroup( C ), [ g ] ) );
    else 
        Add( cosrep, g );
    fi;
    leftover := Difference( leftover, 
    Flat( List ( Set( List( cosrep, h -> RightCoset( C, h ) ) ), AsList ) ) );
od;
return C;
end);


#############################################################################
##
#M OnPoints( a, g )
##
## The function OnPoints(a,g) computes the conjugate a^g where a is an 
## element of the group ring FG and g an element of G. You can use the
## notation a^g to compute it as well.
##
InstallMethod( \^,
    "Wedderga: for a group ring element and a group element",
    function( F1, F2 )
      return IsBound( F1!.familyMagma ) and 
             IsIdenticalObj( ElementsFamily( F1!.familyMagma ), F2 );
    end,
    [ IsElementOfFreeMagmaRing, IsMultiplicativeElementWithInverse ], 
    0,
function( a, g )
local   coeffsupp,
        coeff,
        supp,
        lsupp;

coeffsupp := CoefficientsAndMagmaElements(a);
lsupp := Size(coeffsupp)/2;
supp := List( [ 1 .. lsupp ], i -> coeffsupp [ 2*i-1 ]^g );
coeff := List([ 1 .. lsupp ], i -> coeffsupp[2*i] );

return ElementOfMagmaRing( FamilyObj( a ) ,
                                Zero( a ),
                                coeff,
                                supp);
end);


#############################################################################
##
## CyclotomicClasses( q, n )
##
## The function CyclotomicClasses computes the set of the Cyclotomic Classes
## of q module n 
##
InstallMethod( CyclotomicClasses,
    "for pairs of positive integers", 
    true, 
    [ IsPosInt, IsPosInt ], 
    0,
function(q, n)
local   cc,     # List of cyclotomic classes
        ccc,    # Cyclotomic Class
        i,      # Representative of cyclotomic class
        leftover,  # List of integers
        j;      # Integer
        
# Initialization     
if Gcd( q, n ) <> 1 then
    Error("Wedderga: <q> and <n> should be relatively prime"); 
fi;

#Program
cc := [ [0] ];
leftover:=[ 1 .. n-1 ];
while leftover <> [] do
    i := leftover[ 1 ];
    ccc := [ i ];
    j:=q*i mod n;
    while j <> i do
        Add( ccc, j );
        j:=j*q mod n;
    od;
    Add( cc, ccc );
    leftover := Difference( leftover, ccc );
od;
return cc;
end);


#############################################################################
##
## BigPrimitiveRoot( q )
##
## The function BigPrimitiveRoot computes a primitive root of the finite 
## field of order q.
##
InstallMethod( BigPrimitiveRoot,
    "for a prime power", 
    true, 
    [ IsPosInt ], 
    0,
function(q)
local   Fq,      # The finite field of order q
        factors, # prime factors of q
        p,       # The only prime divisor of q
        o,       # q = p^o
        cp,      # The conway polynomial 
        pr;      # A primitive root of Fq 

# Initialization     
if not IsPrimePowerInt( q ) then
    Error("Wedderga: The input must be a prime power!!!"); 
fi;

#Program
if q<=2^16 then
    Fq := GF(q);
    pr := PrimitiveRoot(Fq);
else
    factors := FactorsInt(q);
    p:=factors[1];
    o:=Size(factors);
    # If q^o is too big then gap never finish to compute the ConwayPolynomial.
    # If there is not cheap ConwayPolynomial, random primitive will be enough
    if IsCheapConwayPolynomial(p,o) then
      cp := ConwayPolynomial(p,o);
    else
      cp := RandomPrimitivePolynomial(p,o);  
    fi;  
    Fq := GF(p, cp);
    pr := RootOfDefiningPolynomial(Fq);
fi;
return pr;
end);


#############################################################################
##
## BigTrace( o, Fq, a )
##
## The function BigTrace returns the trace of the element a in the field 
## extension Fq^o/Fq.
##
InstallMethod( BigTrace,
    "for elements of finite fields", 
    true, 
    [ IsPosInt, IsField, IsObject ], 
    0,
function( o, Fq, a )
local   q,      # The order of the field Fq
        t, y,   # Elements of finite field
        i;      # Integer
        
# Program               
q := Size(Fq);
t := a;
y := a;
for i in [ 1 .. o-1 ] do
    y := y^q;
    t := t + y;
od;
if not(IsFFE(t)) then
    t := ExtRepOfObj(t)[1];
fi;
  return t;
end);


#############################################################################
##
#P IsStronglyMonomial( G )
## ## The property checks whether a group is strongly monomial
##
InstallMethod( IsStronglyMonomial,
    "for finite groups",
    true,
  [ IsGroup ],
  0,
function( G )

local QG ;
 if IsFinite(G) then
#  if IsSupersolvable(G) or IsAbelian(DerivedSubgroup(G)) then
   if IsAbelian(SupersolvableResiduum(G)) then
      return true;
  elif IsMonomial(G) then
      QG := GroupRing( Rationals, G );
      # enforce setting of IsStronglyMonomial
      SSPNonESSPAndTheirIdempotents( QG );
      return IsStronglyMonomial(G);
  else
      return false;
  fi;
else
  Error("Wedderga: The input should be a finite group\n");
fi;
end);


#############################################################################
##
#M IsCyclotomicClass( q, n, C )
##
## The function IsCyclotomicClass checks if C is a q-cyclotomic class module n
##
InstallMethod( IsCyclotomicClass,
	"for two coprime positive integers and a list of integers", 
	true, 
	[ IsPosInt, IsPosInt, IsList ], 
	0,
function( q, n, C )
    local c,C1,j;

c:=C[1];

if n=1 then 
    return C=[0];
elif Gcd(q,n)<> 1 then
    Error("Wedderga: <q> and <n> should be coprime");
elif c >= n then
    return false;
else 
    C1:=[c];
    j:=q*c mod n;
    while j <> c do
      Add( C1, j );
      j:=j*q mod n;
    od;  
    return Set(C)=Set(C1);
fi;

end);


#############################################################################
##
#P IsCyclGroupAlgebra( FG )
##  
## The function checks whether a group is strongly monomial
##
InstallMethod( IsCyclGroupAlgebra, 
    "for semisimple group algebras", 
    true, 
    [ IsGroupRing ], 
    0,
function( FG )

local G ;

G := UnderlyingMagma(FG);

if not(IsSemisimpleFiniteGroupAlgebra(FG) or IsSemisimpleZeroCharacteristicGroupAlgebra(FG)) then 
  Error("Wedderga: The input should be a semisimple group algebra \n",
        "over a finite or a zero characteristic field \n");
fi;

if IsSemisimpleFiniteGroupAlgebra(FG) or 
   IsStronglyMonomial(G) or 
   ForAll( WeddDecomp(FG), x-> not IsList(x) ) then 
    return true;
else
    return false;
fi;
end);


#############################################################################
##
## SizeOfSplittingField( char, p )  
##
## The function SizeOfSplittingField returns the SizeOfFieldOfDefinition 
## for the character char and the prime p
##
InstallMethod( SizeOfSplittingField,
    "for character of finite group and prime number", 
    true, 
    [ IsCharacter, IsPosInt ], 
    0,
function( char, p )

local size,    #  The power of p
      cc,      #  Conjugacy Classes
      value,   #  Element in ValuesOfClassFunction(char)
      m,       #  counter
      image;   #  Galois Group

size := 1;
cc := ConjugacyClasses( UnderlyingCharacterTable( char ) );
for value in ValuesOfClassFunction(char) do
   m:= 1;
   image:= GaloisCyc( value, p );
   while image <> value do
       m:= m+1;
       image:= GaloisCyc( image, p );
   od;      
   size := Lcm( size , m );
od;

return p^size;

end);


#############################################################################
##
## SquareRootMod( a, q )  
##
## The function SquareRootMod solves the equation x^2 = a (mod q) 
## using the algorithm of Tonelli-Shanks
##
InstallMethod( SquareRootMod,
    "for positive integers", 
    true, 
    [ IsPosInt, IsPosInt ], 
    0,
function( a, q )

local x,Q,z,bool,c,i,t,M,b;

# Initialization  
if not IsPrimePowerInt(q)
	then Error("Wedderga: The input needs to be a prime power integer","\n");
fi;   
if IsEvenInt(q) 
	then Error("Wedderga: The input needs to be an odd integer","\n"); 
fi;
if Legendre(a,q)<1
	then Error("Wedderga: The input needs to be a square mod q","\n"); 
fi;

# Program
if (q mod 4) = 3
	then x := a^((q+1)/4) mod q;
else
	Q:=Product(Filtered(Factors(q-1),IsOddInt)); #take the odd part of q-1
	
	# take integer z such that Jacobi(z,q)=-1
	bool := false;
	while bool = false do
		z:=Random([1..q-1]);
		if Jacobi(z,q)=-1 then bool:=true; fi;
	od;
	c := z^Q mod q;

	x := a^((Q+1)/2) mod q;
	t := a^Q mod q;
	M := LogInt((q-1)/Q,2);

	while (t mod q) <> 1 do
		# find lowest o<i<M such that t^(2^i) moq q = 1
		i := 1;
		while (t^(2^i) mod q) <> 1 do
			i := i+1;			
		od;

		b := c^(2^(M-i-1)) mod q;
		x := x*b mod q;
		t := t*b^2 mod q;
		c := b^2 mod q;
		M := i;
	od;	
fi;

return x;
end);


#############################################################################
##
## SquaresMod( q )  
##
## The function SquaresMod returns an integer a such that both a and -1-a 
## are squares modulo q (odd prime power)
##
InstallMethod( SquaresMod,
    "for a positive integer", 
    true, 
    [ IsPosInt], 
    0,
function( q )

local a,bool;

# Initialization
if not IsPrimePowerInt(q)
	then Error("Wedderga: input needs to be a prime power integer","\n");
fi;
if IsEvenInt(q) 
	then Error("Wedderga: input needs to be odd","\n");
fi;

# Program 
bool := false;

while bool = false do
	a:=Random([1..q-1]);
	if Jacobi(a,q)=1 and Jacobi(-1-a,q)=1 then bool:=true; fi;
od;

return a;
end);

 
#############################################################################
##
## SolveEquation2@wedderga ( q )  
##
## The function SolveEquation2@wedderga returns a in GF(q) satisfying 
## a^((q-1)/2)=-1=-Z(q)^0 when q is odd
##
InstallMethod( SolveEquation2@,
    "for a positive integer", 
    true, 
    [ IsPosInt], 
    0,
function( q )

local a,F,bool;

# Initialization

if not IsPrimePowerInt(q) 
	then Error("Wedderga: input needs to be a prime power integer","\n");
fi;
if IsEvenInt(q) 
	then Error("Wedderga: input needs to be odd","\n"); 
fi;

# Program

bool := false;
F:=GF(q);

while bool = false do
	a:=Random(F);
	if a^((q-1)/2) = -Z(q)^0 then bool:=true; fi;
od;

return a;
end);


#############################################################################
##
## SolveEquation3@wedderga( q )  
##
## The function SolveEquation3@wedderga returns b in GF(q^2) satisfying 
## b^((q^2-1)/2)=-1=-Z(q^2)^0 when q is odd
##
InstallMethod( SolveEquation3@,
    "for a positive integer", 
    true, 
    [ IsPosInt], 
    0,
function( q )

local b,F,bool;

# Initialization
if not IsPrimePowerInt(q)
	then Error("Wedderga: input needs to be a prime power integer","\n");
fi;
if IsEvenInt(q) 
	then Error("Wedderga: input needs to be odd","\n"); 
fi;

# Program
bool := false;
F:=GF(q^2);

while bool = false do
	b:=Random(F);
	if b^((q^2-1)/2) = -Z(q^2)^0 then bool:=true; fi;
od;

return b;
end);


#############################################################################
##
## SolveEquation@wedderga( F )  
##
## The function SolveEquation@wedderga returns x and y<>0 in a finite field 
## (of odd characteristic q) F satisfying x^2+y^2=-1=-Z(q^m)^0 
##
InstallMethod( SolveEquation@,
    "for a field", 
    true, 
    [ IsField], 
    0,
function( F )

local q,m,x,y,a,b;

# Initialization
if not IsFinite(F) 
	then Error("Wedderga: input needs to be finite","\n"); 
fi;

q := Characteristic(F);
m := LogInt(Size(F),q);
# F is finite field GF(q^m)

if IsEvenInt(q) 
	then Error("Wedderga: input needs to be of odd characteristic","\n"); 
fi;

# Program
if (q mod 4) = 1
	then 
			x := 0*Z(q^m); 
			a := SolveEquation2@(q);
			y := a^((q-1)/4);
elif (m mod 2) = 0
	then 
			x := 0*Z(q^m); 
			b := SolveEquation3@(q);			
			y := b^((q^2-1)/4);
else
			a := SquaresMod(q);
			x := (Z(q^m)^0)*SquareRootMod(a,q);
			y := (Z(q^m)^0)*SquareRootMod((-1-a) mod q,q);
fi;

return [x,y];

end);



#############################################################################
##
## PrimRootOfUnity( F,n )  
##
## The function PrimRootOfUnity returns returns a n-th primitive root of unity in F
InstallMethod( PrimRootOfUnity,
    "for a field and an positive integer", 
    true, 
    [ IsField, IsPosInt], 
    0,
function( F,n)
  local m,q;
	q := Size(F);
	m := First([1..1000],m->RemInt(q^m,n)=1);
  return BigPrimitiveRoot(q^m)^((q^m-1)/n);
end);


#############################################################################
##
## MakeMatrixByBasis( f,B )  
##
## The function MakeMatrixByBasis represents a linear mapping f with respect 
## to a basis B
InstallMethod( MakeMatrixByBasis,
    "for a mapping and a basis", 
    true, 
    [ IsMapping, IsBasis], 
    0,
function( f,B)
	local M,b;
	M:=[];
	for b in B do
		Add(M,Coefficients(B,Image(f,b)));
	od;	
	return TransposedMat(M);
end);



#############################################################################
##
## ReturnGalElement( e,E,H,K,F1,xi )  
##
## We know that E/H = Gal(F1/F2). And that F1=F(xi).
## The action is defined by the induced conjugation action of E on H/K.
## This function returns an element e in E as element in the Galois group.
InstallMethod( ReturnGalElement,
    "for subgroups and a field", 
    true, 
    [ IsObject , IsGroup, IsGroup, IsGroup, IsField, IsObject], 
    0,
function( e,E,H,K,F1,xi)
	local f,epi,h,hK,xK,i,image;

	epi := NaturalHomomorphismByNormalSubgroup(H,K);
	hK:=MinimalGeneratingSet(Image(epi,H))[1];
	for h in H do
		if Image(epi,h) = hK then break; fi;
	od;
	xK := Image(epi,h^e);
	i:=1;
	while xK <> hK^i do
		i:=i+1;			
	od;	
	
	image := [xi^i];
	f := AlgebraHomomorphismByImages(F1,F1,[xi],image);

	return f;
end);


#############################################################################
##
## LeftMultiplicationBy(x,F1)
##
## Returns the mapping F1->F1 associated with left multiplication by x
InstallMethod( LeftMultiplicationBy,
    "for fields", 
    true, 
    [IsObject, IsField], 
    0,
function(x,F1)
	return MappingByFunction(F1,F1,y->x*y);
end);


#############################################################################
##
## MakeLinearCombination( FE, coef, supp)
##
## Makes the linear combination of coef and supp in FE
InstallMethod( MakeLinearCombination,
    "for algebra, a list and a list", 
    true, 
    [IsAlgebra, IsList, IsList], 
    0,
function( FE, coef, supp)
	local x,i;

	x := Zero(FE);

	for i in [1..Size(coef)] do
		x := x+coef[i]*supp[i];
	od;

	return x;
end);

#############################################################################
##
## Product3Lists( L )  
##
## The function Product3Lists returns the product [a_i*b_j*c_k] of 3 lists 
## [a_1,a_2,...], [b_1,b_2,...] and [c_1,c_2,...]
##
InstallMethod( Product3Lists,
    "for a list", 
    true, 
    [ IsList], 
    0,
function( L )

local I,i,j,k;

I := [];
for i in L[1] do
	for j in L[2] do
		for k in L[3] do
			Add(I,i*j*k);
		od;
	od;
od;

return I;
end);


#############################################################################
##
## IsTwistingTrivial(G,H,K) 
##
## The function IsTwistingTrivial checks if twisting of the simple algebra 
## associated with the strong Shoda Pair (H,K) is trivial
## This can be done in the rational group ring QG
InstallMethod( IsTwistingTrivial,
    "for subgroups", 
    true, 
    [IsGroup, IsGroup, IsGroup], 
    0,
function(G,H,K)
	local A,QG;

	QG:=GroupRing(Rationals,G);
	A:=SimpleAlgebraByStrongSPInfo(QG,H,K);
	if Size(A)=4 and A[4][3]<>0 then return false;
	elif Size(A)=5 and ( not(IsZero(A[5])) or not(IsZero(List(A[4],x->x[3]))) ) then return false;
	fi;
	return true;
end);


#############################################################################
##
#E
##
