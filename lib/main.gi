#############################################################################
##
#W  main.gi               The Wedderga package            Osnel Broche Cristo
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                              Ángel del Río
##
#H  $Id$
##
#############################################################################


#############################################################################
##
#A WedderburnDecomposition( FqG )
##
## The function WedderburnDecomposition computes the Wedderburn 
## decomposition of the group algebra FG over the finite field 
## or the field of rationals
##
InstallMethod( WedderburnDecomposition, 
    "for semisimple rational or finite group algebra", 
    true, 
    [ IsGroupRing ], 
    0,
function( FG )
local   A,      # Simple algebra
        i,      # Counter
        output;
output := [];

if IsSemisimpleFiniteGroupAlgebra( FG ) then
  if not(StronglyShodaPairsAndIdempotents(FG).StronglyMonomial) then 
    Print("Warning!!\nThe direct product of the output is a PROPER direct factor of the input! \n");
  fi;
  for i in StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs do
    A := SimpleAlgebraByStronglySPNC( FG, i[ 1 ], i[ 2 ], i[ 3 ][ 1 ]);
    Append(output, List(i[3], j -> A ) );
  od;
  return output;
  
elif IsSemisimpleRationalGroupAlgebra( FG ) then
  if not(StronglyShodaPairsAndIdempotents(FG).StronglyMonomial) then 
    Print("Warning!!\nThe direct product of the output is a PROPER direct factor of the input! \n");
  fi;
  for i in StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs do
    A := CrossedProductBySSP( UnderlyingGroup( FG ), i[ 1 ], i[ 2 ] );
    Add(output, A );
  od;
  return output;
  
else
  Error("Wedderga: <FG> must be a semisimple group algebra over rationals or over finite field!!!");
fi;  
end);


#############################################################################
##
#A WedderburnDecompositionInfo( FG ) 
##
## The function WedderburnDecompositionInfo compute the data describing 
## Wedderburn Decomposition of the group algebra FG
##
InstallMethod( WedderburnDecompositionInfo , 
    "for semisimple rational or finite group algebra", 
    true, 
    [ IsGroupRing ], 
    0,
function( FG )
local   G,      # Group
        pairs,  # Strongly Shoda pairs of G
        A,      # Simple algebra
        i,      # Counter
        output;
        
G := UnderlyingGroup(FG);
output := [];

if IsSemisimpleRationalGroupAlgebra( FG ) then

  if not(StronglyShodaPairsAndIdempotents(FG).StronglyMonomial) then 
    Print("Warning!!\nThe direct product of the output is a PROPER direct factor of the input! \n");
  fi;
    
#    if HasStronglyShodaPairs( G ) then
#        pairs := StronglyShodaPairs( G );
#    else
#        pairs := StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs;
#    fi;
      
    for i in StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs do
        Add(output, SimpleAlgebraByStronglySPInfoNC( FG, i[ 1 ], i[ 2 ] ) );
    od;
    return output;
    
elif IsSemisimpleFiniteGroupAlgebra( FG ) then

  if not(StronglyShodaPairsAndIdempotents(FG).StronglyMonomial) then 
    Print("Warning!!\nThe direct product of the output is a PROPER direct factor of the input! \n");
  fi;
    for i in StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs do
        A := SimpleAlgebraByStronglySPInfoNC( FG, i[ 1 ], i[ 2 ], i[ 3 ][ 1 ]);
        Append(output, List(i[3], j -> A ) );
    od;
    return output;
    
else

    Error("Wedderga: <FG> must be a semisimple algebra over rationals or over finite field!!!");

fi;

end); 


#############################################################################
##
#O SimpleAlgebraByStronglySP( QG, K, H )
##
## The function SimpleAlgebraByStronglySP computes the simple algebras 
## QG*e( G, K, H) if ( K, H ) is a SSP of G 
## This version does not check the input
##
InstallMethod( SimpleAlgebraByStronglySP, 
    "for semisimple rational group algebras", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function( QG, K, H)
if  IsStronglyShodaPair( UnderlyingGroup( QG ), K, H ) then
    return SimpleAlgebraByStronglySPNC( QG, K, H );
else
    return fail;
fi;
end);




#############################################################################
##
#O SimpleAlgebraByStronglySPNC( QG, K, H )
##
## The function SimpleAlgebraByStronglySPNC computes simple algebras 
## QG*e( G, K, H), for ( K, H ) a SSP of G 
## This version does not check the input
##
InstallMethod( SimpleAlgebraByStronglySPNC, 
    "for semisimple rational group algebras", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function( QG, K, H)
local   G,          # Underlying group
        N,          # Normalizer of H in G
        ind,        # Index of N in G
        NH,         # NH/H
        KH,         # K/H
        NdK,        # N/K
        k,          # Generator of K/H
        ok,         # Order of k
        Potk,       # List of powers of k
        Epi,        # N --> N/H
        Epi2,       # NH --> NH/KH
        i,          # Loop controller
        act,        # Action for the crossed product
        coc;        # Twisting for the crossed product
        
G := UnderlyingGroup( QG );
N   := Normalizer(G,H);
ind := Index(G,N);
if N=K then
    ok := Index( K, H );
    if ind=1 then # G=N
        Info( InfoPCI, 2, "N_G(H) = K = G, returning CF(", ok, ")");
        return CF(ok);
    else
        Info( InfoPCI, 2, "N_G(H) = K <> G, returning M_", 
              ind, "( CF(", ok, ") )");
        return FullMatrixAlgebra( CF(ok), ind );
    fi;                          
else # if N_G(H) <> K
    Epi := NaturalHomomorphismByNormalSubgroup( N, H ) ;
    NH  := Image(Epi,N);
    KH  := Image(Epi,K);
    repeat
        k  := Random(KH);
        ok := Order(k);
    until ok = Size(KH);
    Potk:= [ k ];
    for i in [ 2 .. ok ] do
        Potk[i] := Potk[i-1]*k; 
    od;
    Epi2:=NaturalHomomorphismByNormalSubgroup( NH, KH ) ;
    NdK:=Image(Epi2,NH);
        
      act := function(a) 
             return MappingByFunction( CF(ok), CF(ok), x -> 
               GaloisCyc(x, Position(Potk,k^PreImagesRepresentative(Epi2,a))));
             end;
               
      coc := function(a,b)
             return E(ok)^Position( Potk,
                                    PreImagesRepresentative(Epi2,a*b)^-1 *
                                    PreImagesRepresentative(Epi2,a) *
                                    PreImagesRepresentative(Epi2,b) );
             end;        
    if ind=1 then
      Info( InfoPCI, 2, "N_G(H) <> K, returning crossed product");
      return CrossedProduct(CF(ok), NdK, act, coc);
    else
      Info( InfoPCI, 2, 
        "N_G(H) <> K, returning matrix algebra over crossed product");
      return FullMatrixAlgebra( CrossedProduct(CF(ok), NdK, act, coc), ind );
    fi;  
fi;      
end);



#############################################################################
##
#O SimpleAlgebraByStronglySP( FqG, K, H, c ) 
##
## The function SimpleAlgebraByStronglySP verifies if ( H, K ) is a SSP of G and
## c is a cyclotomic class of q=|Fq| module n=[K:H] containing generators
## of K/H, and in that case computes the simple algebra  FqG*e( G, K, H, c)
##
InstallMethod( SimpleAlgebraByStronglySP, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
    0,
function( FqG, K, H, c )
local   G,      # Group
        Fq,     # Field
        n;      # Index of H in K

G := UnderlyingGroup( FqG );
Fq := LeftActingDomain( FqG );
n := Index( K, H );

if Gcd( c[ 1 ], n ) = 1 and c in CyclotomicClasses( Size( Fq ), n ) and 
                            IsStronglyShodaPair(G, K, H ) then
    return SimpleAlgebraByStronglySPNC( FqG, K, H, c );
    
else

   Error("Wedderga: The input is not appropriate!!!\n");

fi;

end);


#############################################################################
##
#O SimpleAlgebraByStronglySPNC( FqG, K, H, c )
##
## The function SimpleAlgebraByStronglySPNC computes simple algebras 
## FqG*e( G, K, H, c), for ( H, K ) a SSP of G and c a cyclotomic class 
## of q=|Fq| module n=[K:H] containing generators of K/H.
## This version does not check the input
##
InstallMethod( SimpleAlgebraByStronglySPNC, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
    0,
function( FqG, K, H, c )
local   G,          # Group
        N,          # Normalizer of H in G
        ind,        # index of K in G
        epi,        # N -->N/H
        QNH,        # N/H
        QKH,        # K/H
        gq,         # Generator of K/H
        C1,         # Cyclotomic class of q module n in N/H
        St,         # Stabilizer of C1 in N/H
        Fq,F,       # Fields
        q,          # Order of Fq
        ord,        # Integer
        factors,    # prime factors of q
        p,          # The only prime divisor of q
        o,          # q = p^o
        E;          # Stabilizer of C1 in G

G := UnderlyingGroup( FqG );
Fq := LeftActingDomain( FqG );
q := Size( Fq );

if G = H then
    return Fq;
fi;

N := Normalizer( G, H );
epi := NaturalHomomorphismByNormalSubgroup( N, H );
QNH := Image( epi, N );
QKH := Image( epi, K );
gq := MinimalGeneratingSet( QKH )[ 1 ];
C1 := Set( List( c, i -> gq^i ) );
St := Stabilizer( QNH, C1, OnSets );
E := PreImage( epi, St );
ord := Size( c )/Index( E, K ) ;
if q^ord <= 2^16 then
    F := GF(q^ord);
else
    factors := FactorsInt(q);
    p:=factors[1];
    o:=Size(factors);
    if IsCheapConwayPolynomial(p,o*ord) then
      F := GF( p, ConwayPolynomial(p,o*ord) );
    else
      F := GF( p, RandomPrimitivePolynomial(p,o*ord) );  
    fi;  
fi;

ind := Index( G, K );
if ind=1 then
    return F;
else
    return FullMatrixAlgebra( F, ind );
fi;  

end);



#############################################################################
##
#O SimpleAlgebraByStronglySPInfo( QG, K, H ) 
##
## The function SimpleAlgebraByStronglySPInfo compute the data describing simple algebras
## QG*e( G, K, H ), for ( H, K ) a SSP of G, but first verify the inputs 
##
InstallMethod( SimpleAlgebraByStronglySPInfo, 
    "for semisimple rational group algebras", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function( QG, K, H )
if  IsStronglyShodaPair( UnderlyingGroup( QG ), K, H ) then
    return SimpleAlgebraByStronglySPInfoNC( QG, K, H );
else
    return fail;
fi;
end);


#############################################################################
##
#O SimpleAlgebraByStronglySPInfoNC( QG, K, H ) 
##
## The function SimpleAlgebraByStronglySPInfoNC compute the data describing simple 
## algebras QG*e( G, K, H ), for ( H, K ) a SSP of G 
##
InstallMethod( SimpleAlgebraByStronglySPInfoNC, 
    "for semisimple rational group algebras", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function( QG, K, H )
local   G,          # Group
        N,          # Normalizer of H in G
        NH,         # NH/H
        KH,         # K/H
        NdK,        # N/K
        k,          # Generator of K/H
        ok,         # Order of k
        Potk,       # List of powers of k
        Epi,        # N --> N/H
        Epi2,       # NH --> NH/KH
        PrimGen,    # Primary set of independent of generators of N/K
        l,          # Length of PrimGen
        Gen,        # Elementary set of independent of generators of
        o,          # Orders of the elemnets of PrimGen
        p,          # Prime divisors of the elements of o
        primes,     # The different elements of p
        lp,         # Length of primes,
        first,      # First Positions of the elements of primes in p,
        next,       # Next Positions of the elements of primes in p,
        g,          # An element of PrimGen
        plus,       # Counter
        newpos,     # A component of next
        gen,        # Preimage of Gen in N/H
        i,ll;       # Controlers
        
    G := UnderlyingGroup( QG );
    if G = H then
        return [ 1, 1, [], [] ];
    fi;
    
    # First one computes an idependent set PrimGen of generators 
    # of a Primary decomposition of N/K
    N   := Normalizer(G,H);
    if N=K then
        ok := Index( K, H );
        return [ Index(G,N), ok, [ ], [ ] ];
    else
        Epi := NaturalHomomorphismByNormalSubgroup( N, H ) ;
        NH  := Image(Epi,N);
        KH  := Image(Epi,K);
        k   := Product(IndependentGeneratorsOfAbelianGroup(KH));
        ok  := Order(k);
        Potk:= [ k ];
        for i in [ 2 .. ok ] do
            Potk[i] := Potk[i-1]*k; 
        od;
        Epi2:=NaturalHomomorphismByNormalSubgroup( NH, KH ) ;
        NdK:=Image(Epi2,NH);
        PrimGen:=IndependentGeneratorsOfAbelianGroup(NdK);
        # Using PrimGen one computes an independent set Gen of
        # generators of an invariant decomposition of N/K
        l := Length( PrimGen );
        o := List( [ 1 .. l ], i -> Order( PrimGen[i] ) );
        p := List( [ 1 .. l ], i -> FactorsInt( o[i] )[1] );
        primes := DuplicateFreeList( p );
        lp:= Length( primes );
        first:=List( [ 1 .. lp ], i -> Position( p, primes[i] ) );
        g := Product( List( first, i -> PrimGen[i] ) );
        Gen:=[ g ];
        ll:=lp;
        plus:=0;
        while ll<l do
            next:=[];
            for i in [ 1 .. lp ] do
                newpos := Position( p, primes[i], first[i]+plus );
                if newpos <> fail then
                    Add( next, newpos );
                fi;
            od;
            g:=Product( List( next, i -> PrimGen[i] ) );
            Add( Gen, g );
            ll:=ll+Length(next);
            plus:=plus+1;
        od;
        gen:=List( [ 1 .. Length(Gen) ], i -> PreImagesRepresentative(Epi2,Gen[i]) );
        return [ Index(G,N), 
                 ok, 
                 List( [1..Length(Gen)],
                   i->[ Order(Gen[i]),
                        # we have a list Potk of powers of k and find the 
                        # position of k^gen[i] in it. Is there better way
                        # to determine j such that k^gen[i] = k^j ?
                        RemInt(Position(Potk,k^gen[i]),ok),
                        RemInt(Position(Potk,gen[i]^Order(Gen[i])),ok) ]),
                 List( [1..Length(Gen)-1], i -> 
                   List( [i+1..Length(Gen)], j -> 
                     RemInt(Position(Potk,Comm(gen[j],gen[i])),ok))) ];
    fi;
end);


#############################################################################
##
#O SimpleAlgebraByStronglySPInfo( FqG, K, H, c )
##
## The function SimpleAlgebraByStronglySPInfo compute the data describing simple algebra 
## FqG*e( G, K, H, c) for ( H, K ) a SSP of G and c a cyclotomic class 
## of q=|Fq| module n=[K:H], containing generators of K/H, 
## but first verify the inputs 
##
InstallMethod( SimpleAlgebraByStronglySPInfo, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
    0,
function( FqG, K, H, c )  
local   G,      # Group
        Fq,     # Field
        n;      # Index of H in K

G := UnderlyingGroup( FqG );
Fq := LeftActingDomain( FqG );
n := Index( K, H );

if  Gcd( c[ 1 ], n ) = 1 and c in CyclotomicClasses( Size( Fq ), n ) and 
                             IsStronglyShodaPair(G, K, H ) then
    return SimpleAlgebraByStronglySPInfoNC( FqG, K, H, c );
else
    Error("Wedderga: The input is not appropriate!!!\n");
fi;
end);


#############################################################################
##
#O SimpleAlgebraByStronglySPInfoNC( FqG, K, H, c )
##
## The function SimpleAlgebraByStronglySPInfoNC compute the data describing simple 
## algebra FqG*e( G, K, H, c) for ( H, K ) a SSP of G and c a cyclotomic 
## class of q=|Fq| module n=[K:H], containing generators of K/H.
##
InstallMethod( SimpleAlgebraByStronglySPInfoNC, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
    0,
function( FqG, K, H, c )  
local   G,          # Group
        Fq,         # Finite field
        q,          # Order of Fq
        N,          # Normalizer of H in G
        epi,        # N -->N/H
        QNH,        # N/H
        QKH,        # K/H
        gq,         # Generator of K/H
        C1,         # Cyclotomic class of q module n in N/H
        St,         # Stabilizer of C1 in N/H
        ord,        # Integer
        E;          # Stabilizer of C1 in G

G := UnderlyingGroup( FqG );
Fq := LeftActingDomain( FqG );
q := Size( Fq );

if G = H then
  return [ 1, q ];
fi;

N := Normalizer( G, H );
epi := NaturalHomomorphismByNormalSubgroup( N, H );
QNH := Image( epi, N );
QKH := Image( epi, K );
# We guarantee that QKH is cyclic so we can randomly obtain its generator
repeat
  gq := Random(QKH);
until Order(gq) = Size(QKH);
C1 := Set( List( c, ii -> gq^ii ) );
St := Stabilizer( QNH, C1, OnSets );
E := PreImage( epi, St );
ord := q^( Size( c )/Index( E, K ) );

return [ Index( G, K ), ord ];

end);


#############################################################################
##
#A StronglyShodaPairs( G )
##
## The function StronglyShodaPairs computes a list of strongly Shoda pairs 
## of the group G that covers the complete set of primitive central 
## idempotents of the rational group algebra QG realizable by strongly 
## Shoda pairs
##
InstallMethod( StronglyShodaPairs, 
    "for finite group ", 
    true, 
    [IsGroup and IsFinite], 
    0,
function( G )
local   QG;     # Rational Group Algebra
       
QG:=GroupRing( Rationals, G ); 
        
return StronglyShodaPairsAndIdempotents(QG).StronglyShodaPairs;

end);


#############################################################################
##
#A StronglyShodaPairsAndIdempotents( QG )
##
## The attribute StronglyShodaPairsAndIdempotents of the rational group algebra QG 
## returns a record with components StronglyShodaPairs, PrimitiveCentralIdempotents and 
## StronglyMonomial where 
## StronglyShodaPairs = list of SSP that covers the complete set of primitive 
##       central idempotents of QG realizable by SSPs, 
## PrimitiveCentralIdempotents = list of PCIs of QG realizable by SSPs,
## StronglyMonomial := Yes if PrimitiveCentralIdempotents is a complete set of PCIs of QG, No otherwise 
##
InstallMethod( StronglyShodaPairsAndIdempotents, 
    "for rational group algebra", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra ], 
    0,
function( QG )
local   G,          # Group
        DG,         # Derived subgroup of G
        CCS,        # The conjugacy classes of subgroups
        LCCS,       # The length of CCS 
        eGKHs,      # The list of primitive central idempotents 
        SeGKHs,     # The sum of the elements of eGKHs
        H,          # Subgroup of G
        NH,         # Normalizer of H in G
        NHH,        # NH/H
        K,          # Subgroup of G 
        i, j,       # Counters
        idempeGKH,  # A primitive central idempotent of the form e(G,K,H) for (H,K) a SSP
        KH,         # K/H
        KHs;        # The list of SSP

G := UnderlyingGroup(QG);

if HasStronglyShodaPairs( G ) then
    eGKHs := List( StronglyShodaPairs( G ), i -> 
                                     CentralElementBySubgroups( QG, i[1], i[2] ) );
    return rec( 
    StronglyShodaPairs := StronglyShodaPairs( G ), 
    PrimitiveCentralIdempotents := eGKHs,
    StronglyMonomial := IsCompleteSetOfPCIs( QG , eGKHs ) ); 

else

  CCS:=ConjugacyClassesSubgroups(G);
  LCCS:=Length(CCS);
  DG:=DerivedSubgroup(G); 
  KHs:=[];
  eGKHs:=[];
  SeGKHs:=Zero(QG);
  if Size(G)=1 then
    return rec( 
      StronglyShodaPairs := [ [ G, G ] ], 
      PrimitiveCentralIdempotents := [ One(QG)] );
  fi;
  for j in [ LCCS, LCCS-1 .. 1 ] do
    H:=Representative(CCS[j]);        
    if IsSubset(H,DG) then
      if IsCyclic(FactorGroup(G,H)) then 
        idempeGKH:=eG(QG,G,H)[2]; 
        SeGKHs:= SeGKHs + idempeGKH;
        Add(KHs,[G,H]);
        Add(eGKHs,idempeGKH);
      fi;
    else 
      idempeGKH:=SearchingKForSSP(QG,H);
      if idempeGKH<>fail and idempeGKH[2]*SeGKHs=Zero(QG) then 
        SeGKHs:= SeGKHs + idempeGKH[2];
        Add(KHs,idempeGKH[1]);
        Add(eGKHs,idempeGKH[2]);
      fi;
    fi;
    if SeGKHs=One(QG) then
      break;
    fi;
  od;

  #Here finish the main loop

#  if SeGKHs<>One(QG) then 
#    Print(  "Warning!!! Some primitive central idempotents are not realizable ", 
#            "by strongly Shoda pairs!!!\n");
#  fi;

  SetStronglyShodaPairs( G , KHs ); 

  return rec( 
    StronglyShodaPairs := KHs, 
    PrimitiveCentralIdempotents := eGKHs,
    StronglyMonomial := SeGKHs=One(QG) );

fi;

end);

#############################################################################
##
#A StronglyShodaPairsAndIdempotents( FqG )
##
## The attribute StronglyShodaPairsAndIdempotents of the semisimple finite group algebra FqG 
## returns a record with components StronglyShodaPairs, PrimitiveCentralIdempotents and 
## StronglyMonomial where 
## StronglyShodaPairs = list of SSP and cyclotomic classes that covers the set of PCIs of FqG 
##        realizable by SSPs, 
## PrimitiveCentralIdempotents = list of PCIs of FqG realizable by SSPs and cyclotomic classes,
## StronglyMonomial := Yes if PrimitiveCentralIdempotents is a complete set of PCIs of FqG, No otherwise 
##

## The function StronglyShodaPairsAndIdempotents computes the record [SSPs, PCIs], 
## where SSPs is a list of the SSP and cyclotomic classes that covers 
## the complete set of primitive central idempotents,PCIs, 
## of the finite group algebra FqG
##
InstallMethod( StronglyShodaPairsAndIdempotents, 
    "for semisimple finite group algebra", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra ], 
    0,
function( FqG )                
local   G,          # Group
        SSPsG,      # List of strongly Shoda pairs of G
        H,K,        # Subgroups of G
        n,          # Index of H in K
        Fq,         # Field (finite)
        q,          # Order of Fq
        idemp,      # Idempotent eGKHc 
        cc,         # Set of cyclotomic classes of q module n
        lcc,        # Set of cc's 
        i,          # Cyclotomic class of q module
        e,          # The list of primitive central idempotents
        list,       # List SSP and cyclotomic classes
        p,          # Integer
        a,          # Primitive n-th root of 1 in an extension of Fq
        ltrace,     # List of traces of a^c over Fq for c in representatives of cc
        lltrace,    # List of ltrace's for n in setind
        setind,     # Set of n's 
        pos,        # Positions
        j,          # Counter
        o,          # The  multiplicative order of q module n
        lorders,    # Set of o's for various n's
        pr,         # Primitive root of the field of order q^o
        lprimitives,# Set of pr's for o in lorders
        etemp,      # List of idempotents eGKHc for different classes c and fixed K and H
        templist,   # List of some cyclotomic classes
        elmsG,      # Elements of G
        F,          # Family of elements of FqG 
        zero;       # Zero of Fq

# Program
G := UnderlyingGroup( FqG  );
Fq := LeftActingDomain( FqG );
F := FamilyObj(Zero(FqG));
elmsG := Elements(G);
q := Size( Fq );
zero := Zero(Fq);
e := [AverageSum(FqG,G)];
SSPsG := StronglyShodaPairs(G);
list := [ [ SSPsG[ 1 ][1], SSPsG[ 1 ][2], [ [ 0 ] ] ] ];
setind := [];
lltrace := [];
lcc := [];
lorders := [];
lprimitives := [];
for p in [ 2 .. Size(SSPsG) ] do
    H :=SSPsG[p][2];
    K := SSPsG[p][1];
    n := Index(K,H);
    if n in setind then
    # If n is in setind then we just take Cyclotomic Classes and traces 
    # from lcc and lltrace
        pos := Position(setind,n);
        cc := lcc[pos];
        ltrace  := lltrace[pos];
    else
    # Otherwise we compute traces and cyclotomic classes and store them
    # in lltrace and lcc
        cc := CyclotomicClasses(q,n);
        o:=Size(cc[2]);
        if o in lorders then
        # If o is in lorders then a primitive root of 1 is stored in lprimitives
            pr := lprimitives[Position(lorders,o)];
        else
        # Otherwise we compute the primitive root and store it in lprimitives
            pr := BigPrimitiveRoot(q^o);
            Add(lorders,o);
            Add(lprimitives,pr);
        fi;
        a:=pr^((q^o-1)/n);
        ltrace := [];
        for i in cc do
            ltrace[i[1]+1] := BigTrace(o,Fq, a^i[1]);
            for j in i do
                ltrace[j+1] := ltrace[i[1]+1];
            od;
        od;
        Add(lltrace,ltrace);
        Add(lcc,cc);
        Add(setind,n);
    fi;
    etemp := [];
    templist := [];
    for i in cc do
        if Gcd(i[1],n)=1 then
            idemp := CentralElementBySubgroups(FqG, K, H, i, ltrace);
            if not(idemp in etemp) then
                Add(etemp, idemp);
                Add( templist, i );
            fi;
        fi;
    od;
    Append( e, etemp );
    Add( list, [ K, H, templist ] );
od;
return rec( StronglyShodaPairs := list, 
            PrimitiveCentralIdempotents := e, 
            StronglyMonomial := IsCompleteSetOfPCIs ( FqG, e ));
end);



#############################################################################
## 
#F SearchingKForSSP(QG,H)
##
## The following function search an element K such that (K,H) is a SSP
## and returns [ [ K, H ], e( G, K, H ) ] or returns fail, if
## such K doesn't exist
##
InstallGlobalFunction( SearchingKForSSP, function(QG,H)
    local   
        G,          # underlying group of QG
        NH,         # Normalizer of H in G
        Epi,        # NH --> NH/H        
        NHH,        # NH/H
        L,          # <NHH',Z(NHH)>
        Cen,        # Centralizer of L in NHH
        K,          # The subgroup searched
        e,          # e(G,K,H) for some of the searched K
        KH,         # K/H
        X;          # a subset of Cen

        G:=UnderlyingGroup(QG);
        NH:=Normalizer(G,H);
        Epi:=NaturalHomomorphismByNormalSubgroup( NH, H ) ;
        NHH:=Image(Epi,NH);
        L:=ClosureSubgroup( DerivedSubgroup(NHH), Centre(NHH) );
        if IsCyclic(L) then 
            Cen:=Centralizer(NHH,L);
            if IsAbelian(Cen) then
                if IsCyclic(Cen) and Centralizer(NHH,Cen)=Cen then
                    K:=PreImages(Epi,Cen);
                    return eG(QG,K,H);
                else 
                    return fail;
                fi;
            else 
                X:=Difference(Elements(Cen),Elements(L));
                while X<>[] do
                    KH:=ClosureSubgroup( L, [X[1]] );
                    if IsCyclic(KH) and Centralizer(NHH,KH)=KH then
                        K:=PreImages(Epi,KH);
                        return eG(QG,K,H);
                    fi;
                    X:=Difference(X,KH);
                od;
            fi;
        fi;
    return fail;                  
    end);     


#############################################################################
##
#M eG( QG, K, H )
##
## The following function computes e(G,K,H)    
## Note that actually it returns a list of the form [ [K,H], eGKH ]
##
InstallMethod( eG,
    "for pairs of subgroups", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function(QG,K,H)
    local   
        G,      # underlying group of QG
        Eps,    # \varepsilon(K,H), 
        eGKH,   # is the final return that takes partial values 
        NH,     # Normalizer of H in G
        NdK,    # Normalizer of K in G
        RTNH,   # Right transveral of NH in NdK
        nRTNH,  # Cardinal de RTNH
        eGKH1,  # e(NdK,K,H)
        eGKH1g, # eGKH1^g
        i,      # counter 
        g,      # element of G
        RTNdK,  # Right transversal of G/NdK
        nRTNdK, # Cardinal of RTNdK
        zero;   # zero of QG

        Eps:=IdempotentBySubgroups(QG,K,H);
        G:=UnderlyingGroup(QG);
	zero := Zero( QG );
        NH:=Normalizer(G,H);
        if NH=G then
            return [ [ K, H ], Eps ];
        else
            NdK:=Normalizer(G,K);
            RTNH:=RightTransversal(NdK,NH);
            eGKH1:=Sum( List( RTNH,g->Eps^g ) );
            eGKH:=eGKH1;
            if NdK<>G then
                RTNdK:=RightTransversal(G,NdK); 
                nRTNdK:=Length(RTNdK);  
                for i in [ 2 .. nRTNdK ] do
                    g:=RTNdK[i];
                    eGKH1g:=eGKH1^g;
                    if eGKH1*eGKH1g <> zero then 
                        return  fail;
                    else
                        eGKH:= eGKH + eGKH1g;
                    fi;
                od;                    
            fi;
            return [ [ K, H ], eGKH ];
        fi;       
    end);




#############################################################################
##
#O PrimitiveCentralIdempotentsByStronglySP( FG )
##
## The function PrimitiveCentralIdempotentsByStronglySP computes the set of 
## primitive central idempotents of the group algebra FG, realizable by 
## strongly Shoda pairs, where FG is either a rational or finite group algebra
##
InstallGlobalFunction( PrimitiveCentralIdempotentsByStronglySP, 
function( FG )

if not(StronglyShodaPairsAndIdempotents(FG).StronglyMonomial) then 
   Print("Warning!!\nThe output is a NON-COMPLETE list of prim. central idemp.s of the input! \n");
fi;

return StronglyShodaPairsAndIdempotents( FG ).PrimitiveCentralIdempotents; 
end);


#############################################################################
##
#E
##
