#############################################################################
##
#W  auxiliar.gi           The Wedderga package            Osnel Broche Cristo
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                              Ángel del Río
##
#H  $Id$
##
#############################################################################


#############################################################################
##
#P IsSemisimpleRationalGroupAlgebra( FG )
##  
## The function checks whether a group ring is a rational group algebra
##
InstallImmediateMethod( IsSemisimpleRationalGroupAlgebra,
                        IsGroupRing, 
                        0,
    R -> IsRationals(LeftActingDomain(R)) and IsFinite(UnderlyingMagma(R))
); 


#############################################################################
##
#P IsSemisimpleFiniteGroupAlgebra( FG )
##  
## The function checks whether a group ring is a semisimple finite group algebra
##
InstallImmediateMethod( IsSemisimpleFiniteGroupAlgebra,
                        IsGroupRing, 
                        0,
                
function( FG )
local   F,      # Field
        G;      # Group

F := LeftActingDomain( FG );
G := UnderlyingMagma( FG );
       
return IsField( F ) and IsFinite( F ) and IsFinite( G ) and Gcd( Size( F ), Size( G ) )=1;
end); 


#############################################################################
##
#M IsCompleteSetOfPCIs( QG, ListPCIs )
##
## The function IsCompleteSetOfPCIs checks if the sum of the elements of QG is 1
## It is supposed to be used to check if a list of PCIs of QG is complete.
##
InstallMethod( IsCompleteSetOfPCIs,"for list of primitive central idempotents", true, 
[IsFreeMagmaRing,IsList ], 0,
function( QG, ListPCIs )
    local x;
    if not IsSemisimpleRationalGroupAlgebra( QG ) then
        Error("Wedderga: The first argument must be a rational group algebra!!!\n");
    elif not ForAll( ListPCIs, x -> x in QG ) then
        Error("Wedderga: An element of the list of PCIs do not belong to the 1st argument!!!\n");
    else
        return Sum(ListPCIs)=One(QG);
    fi;
end);


#############################################################################
##
#F IsStronglyShodaPair( G, K, H )
##
## The function IsStronglyShodaPair verifies if (H,K) is a SSP 
##
InstallMethod( IsStronglyShodaPair,
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
        KH;

# First verifies if H, K are subgroups of G and K is a normal subgroup of K
if not ( IsSubgroup( G, K ) and IsSubgroup( K, H ) ) then
    Error("Wedderga: Each argument should contain the next one!!!\n");
fi;

if not IsNormal( K, H ) then
    Info(InfoPCI, 2, "Wedderga: The 3rd subgroup is not normal in the 2nd");
    return false;
fi;

# Now, if K/H is the maximal abelian subgroup in N/H,
# where N is the normalizer of H in G

NH:=Normalizer(G,H);

if not(IsNormal( NH, K ) ) then
    Info(InfoPCI, 2, "Wedderga: The 2nd is not normal in the normalizer of 3rd one in the 1st");
    return false;
fi;

Epi:=NaturalHomomorphismByNormalSubgroup( NH, H ) ;
NHH:=Image( Epi, NH ); #It is isomorphic to the factor group NH/H.
KH:=Image( Epi, K ); #It is isomorphic to the factor group K/H.

if not(IsCyclic(KH)) then
    Info(InfoPCI, 2, "Wedderga: The 2nd subgroup over the 3rd one is not cyclic");
    return false;
fi;

if Centralizer( NHH, KH ) <> KH then
    Info(InfoPCI, 2, "Wedderga: The 2nd subgroup over the 3rd one is not cyclic");
    Info(InfoPCI, 2, "Wedderga: The factor group (2nd over 3rd) is not maximal ",
                     "abelian in the normalizer of the 3rd in the 1st");
    return false;
fi;

#Now (SSS3)
QG := GroupRing( Rationals, G );
Eps := IdempotentBySubgroups( QG, K, H );
NdK := Normalizer( G, K );

if NdK<>G then
    RTNH := RightTransversal( NdK, NH );
    eGKH1 := Sum( List( RTNH, g -> Conjugate( QG, Eps, g ) ) );
    RTNdK := RightTransversal( G, NdK );
    nRTNdK:=Length(RTNdK);
    for i in [ 2 .. nRTNdK ] do
        g:=RTNdK[i];
        if not IsZero( eGKH1*Conjugate(QG,eGKH1,g) ) then
            Info(InfoPCI, 2, "Wedderga: The conjugates of epsilon are not orthogonal");
            return  false;
        fi;
    od;
fi;

return true;

end);


#############################################################################
##
#M CentralizerInUnderlyingGroup( FG, a )
##
##
## The function CentralizerInUnderlyingGroup computes the centralizer of an
## element of a group ring in the underlying group
##
InstallMethod(  CentralizerInUnderlyingGroup,
    "for a group ring and a group ring element",
    IsCollsElms, 
    [ IsGroupRing, IsElementOfFreeMagmaRing ], 
    0,
function( FG, a ) 
local   G,
        C,
        Leftover,
        Classes,
        one,
        ElemG,
        g;

G := UnderlyingMagma( FG );
one := Identity( G );
C := Subgroup( G, [one] );
Classes := [ one ];
ElemG := Elements(G);
Leftover := Difference(G,C);
while Leftover <> [] do
    g := Leftover[1];
    if Conjugate( FG, a, g ) = a then 
        C := Subgroup( G, Union( C, [g] ) );
    else 
        Add( Classes, g );
    fi;
    Leftover := Difference( G, Union( List( Classes, i -> RightCoset( C, i ) ) ) );
od;
return C;
end);


#############################################################################
##
#M Conjugate( FG, a, g )
##
## The function Conjugate computes the conjugate a^g where a is an element 
## of the group ring FG and g an element of G.
##
InstallMethod(  Conjugate,
    "for a group ring, group ring element and a group element",
    true, 
    [ IsGroupRing, IsElementOfFreeMagmaRing, IsObject ], 
    0,
function( FG, a, g )
local   coeffsupp,
        coeff,
        supp,
        lsupp;

# may be these conditions can be checked in the filter, but there is another
# option to install two-argument version for a^g via OnPoints later
    if not a in FG then
        Error("The second argument must be an element of the first one !!!\n");
    fi;

    if not g in UnderlyingMagma( FG ) then
        Error("The last argument should belong to the group basis !!!\n");
    fi;

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
if Gcd( q, n ) <>1 then
    Error("The inputs should be coprime!!!"); 
fi;

#Program
cc := [[0]];
leftover:=[1..n-1];
while leftover <> [] do
    i := leftover[1];
    ccc := [i];
    j:=q*i mod n;
    while j <> i do
        Add(ccc,j);
        j:=j*q mod n;
    od;
    Add(cc,ccc);
    leftover:=Difference(leftover,ccc);
od;
return cc;
end);


#############################################################################
##
## BigPrimitiveRoot( q )
##
## The function BigPrimitiveRoot computes a primitive root of the finite field
## of order q.
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
    if IsCheapConwayPolynomial(p,o) then
      cp := ConwayPolynomial(p,o);
    # If q^o is too big then gap never finish to compute the ConwayPolynomial
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
function(o, Fq, a)
local   q,      # The order of the field Fq
        t, y,   # Elements of finite field
        i;      # Integer
        
# Program               
q := Size(Fq);
t := a;
y := a;
for i in [1..o-1] do
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
#E
##
