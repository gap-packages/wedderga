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
#P  IsSemisimpleRationalGroupAlgebra( FG )
##  
##  The function checks whether a group ring is a rational group algebra
InstallImmediateMethod( IsSemisimpleRationalGroupAlgebra,
                        IsGroupRing, 
                        0,
    R -> IsRationals(LeftActingDomain(R)) and IsFinite(UnderlyingMagma(R))
); 

#############################################################################
##
#P  IsSemisimpleFiniteGroupAlgebra( FG )
##  
##  The function checks whether a group ring is a semisimple finite group algebra
InstallImmediateMethod( IsSemisimpleFiniteGroupAlgebra,
                        IsGroupRing, 
                        0,
                
function( FG )
local   F,      # Field
        G;      # Group

F := LeftActingDomain( FG );
G := UnderlyingMagma( FG );
       
return IsField( F ) and IsFinite( F ) and IsFinite( G ) and Gcd( Size( F ), Size( G ) ) =1;
end); 

#############################################################################
##
#M  IsCompleteSetOfPCIs( QG, ListPCIs )
##
##
##  The function IsCompleteSetOfPCIs checks if the sum of the elements of QG is 1
##  It is suppose to be used to check if a list of PCIs of QG is complete.
##
InstallMethod( IsCompleteSetOfPCIs,"for list of primitive central idempotents", true, 
[IsFreeMagmaRing,IsList ], 0,
function( QG, ListPCIs )
    local x;
    if not IsSemisimpleRationalGroupAlgebra( QG ) then
        Error("The first argument must be a rational group algebra!!!");
    elif not ForAll(ListPCIs, x -> x in QG) then
        Error("An element of the list of PCIs do not belong to the 1st argument!!!");
    else
        return Sum(ListPCIs)=One(QG);
    fi;
end);

#############################################################################
##
##  The function IsStronglyShodaPair verifies if (H,K) is a SSP 
##
#F IsStronglyShodaPair( G, K, H )
##
InstallMethod( IsStronglyShodaPair,
                "for a group and two subgroups", 
                true,
                [ IsGroup, IsGroup, IsGroup ], 
                0,
function( G, K, H )
local   QG,
        zero,
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
if not ( IsSubgroup( G, K ) and IsSubgroup( K, H ) and IsNormal( K, H ) ) then
    Error("Each input should contain the next one and the last one should be normal ",
            "in the second!!!\n");
fi;

# Now, if K/H is the maximal abelian subgroup in N/H,
# where N is the normalizer of H in G

NH:=Normalizer(G,H);

if not(IsNormal( NH, K ) ) then
    Print("The second input must be a normal subgroup in the normalizer of third one ",
            "in the first\n");
    return false;
fi;

Epi:=NaturalHomomorphismByNormalSubgroup( NH, H ) ;
NHH:=Image( Epi, NH ); #It is isomorphic to the factor group NH/H.
KH:=Image( Epi, K ); #It is isomorphic to the factor group K/H.

if not(IsCyclic(KH)) then
    Print("The second input over the third one should be cyclic\n");
    return false;
fi;

if Centralizer( NHH, KH ) <> KH then
    Print("The factor group (second input over third one) is not maximal abelian ",
            "on the normalizer of the third one in the first\n");
    return false;
fi;

#Now (SSS3)
QG := GroupRing( Rationals, G );
zero := Zero( QG );
Eps := Epsilon( QG, K, H );
NdK := Normalizer( G, K );

if NdK<>G then
    RTNH := RightTransversal( NdK, NH );
    eGKH1 := Sum(List(RTNH,g->Conjugate(QG,Eps,g)));
#    eGKH := eGKH1;
    RTNdK := RightTransversal( G, NdK );
    nRTNdK:=Length(RTNdK);
    for i in [ 2 .. nRTNdK ] do
        g:=RTNdK[i];
        if eGKH1*Conjugate(QG,eGKH1,g) <> zero then
            Print("The conjugates of epsilon are not orthogonal \n");
            return  false;
        fi;
    od;
fi;

return true;

end);

#############################################################################
##
#M  CentralizerG( FG, a )
##
##
##  The function CentralizerG computes the centralizer of an element of a group ring in the 
##  defining underlying group
##
InstallMethod(  CentralizerG,
                "for an group ring element",
                true, 
                [ IsFreeMagmaRing, IsElementOfFreeMagmaRing ], 
                0,
function( FG, a ) 
local   G,
        C,
        Leftover,
        Classes,
        one,
        ElemG,
        g;

    if not IsGroupRing( FG ) then
        Error("The first argument must be a group ring !!!\n");
    fi;
    
    if not a in FG then
        Error("The second argument must be an element of the first one !!!\n");
    fi;

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
#M  Conjugate( FG, a, g )
##
##
##  The function Conjugate computes the conjugate a^g where a is an element of the group ring
##  FG and g an element of G.
##
InstallMethod(  Conjugate,
                "for a group ring element and a group element",
                true, 
                [ IsFreeMagmaRing, IsElementOfFreeMagmaRing, IsObject ], 
                0,
function( FG, a, g )
local   coeffsupp,
        coeff,
        supp,
        lsupp;

    if not IsGroupRing(FG) then
        Error("The first argument must be a group ring !!!\n");
    fi;
    
    if not a in FG then
        Error("The second argument must be an element of the first one !!!\n");
    fi;

    if not g in UnderlyingMagma( FG ) then
        Error("The last argument should belong to the group basis !!!\n");
    fi;

coeffsupp := CoefficientsAndMagmaElements(a);
lsupp := Size(coeffsupp)/2;
supp := List([1..lsupp],i->coeffsupp[2*i-1]^g);
coeff := List([1..lsupp],i->coeffsupp[2*i]);

return ElementOfMagmaRing( FamilyObj(Zero(FG)),
                               Zero(FG),
                               coeff,
                               supp);
end);

#############################################################################
##
## The function CyclotomicClasses computes the set of the Cyclotomic Classes
## of q module n 
##
## CyclotomicClasses( q, n )
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
## The function BigPrimitiveRoot computes a primitive root of the finite field
## of order q.
##
## BigPrimitiveRoot( q )
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
    Error("The input must be a prime power!!!"); 
fi;

#Program
if q<=2^16 then
    Fq := GF(q);
    pr := PrimitiveRoot(Fq);
else
    factors := FactorsInt(q);
    p:=factors[1];
    o:=Size(factors);
    cp := ConwayPolynomial(p,o);
# If q^o is too big then gap never finish to compute the ConwayPolynomial
    Fq := GF(p, cp);
    pr := RootOfDefiningPolynomial(Fq);
fi;
return pr;
end);

#############################################################################
##
## The function BigTrace returns the trace of the element a in the field extension Fq^o/Fq.
##
## BigTrace( o, Fq, a)
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
