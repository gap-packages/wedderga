##############################################################################
##
##  The function WedderburnDecomposition computes the Wedderburn 
##  Decomposition of the finite group algebra FqG
##
#A WedderburnDecomposition( FqG )
##
InstallMethod( WedderburnDecomposition, 
                "for semisimple finite group algebras", 
                true, 
                [ IsSemisimpleFiniteGroupAlgebra ], 
                0,
function( FqG )
local   A,      # Simple algebra
        i,      # Counter
        output;
        
output := [];
for i in StronglyShodaPairsAndIdempotents( FqG ).StronglyShodaPairs do
    A := SimpleAlgebraNC( FqG, i[ 1 ], i[ 2 ], i[ 3 ][ 1 ]);
    Append(output, List(i[3], j -> A ) );
od;

return output;

end);

#############################################################################
##
##  The function WedderburnDecompositionInfo compute the data describing 
##  Wedderburn Decomposition of the group algebra FG
##
#A WedderburnDecompositionInfo( FG ) 
##
InstallMethod( WedderburnDecompositionInfo , 
                "for semisimple rational or finite group algebra ", 
                true, 
                [ IsGroupRing ], 
                0,
function( FG )
local   G,      #Group
        pairs,  # Strongly Shoda pairs of G
        A,      # Simple algebra
        i,      # Counter
        output;
        
G:=UnderlyingMagma(FG);
output := [];

if IsSemisimpleRationalGroupAlgebra( FG ) then
    
    if Tester( StronglyShodaPairs )( G ) then
        pairs := StronglyShodaPairs( G );
    else
        pairs := StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs;
    fi;
      
    for i in pairs do
        Add(output, SimpleAlgebraInfoNC( FG, i[ 1 ], i[ 2 ] ) );
    od;
    
    return output;
    
elif IsSemisimpleFiniteGroupAlgebra( FG ) then

    for i in StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs do
        A := SimpleAlgebraInfoNC( FG, i[ 1 ], i[ 2 ], i[ 3 ][ 1 ]);
        Append(output, List(i[3], j -> A ) );
    od;

    return output;
    
fi;

Error("The input must be a semisimple algebra over rational or finite field!!!\n");

end); 

#############################################################################
##
##  The function SimpleAlgebra verifies if ( H, K ) is a SSP of G and
##  c is a cyclotomic class of q=|Fq| module n=[K:H] containing generators of K/H,
##  and in that case computes the simple algebra  FqG*e( G, K, H, c)
##
#O SimpleAlgebra( FqG, K, H, c ) 
##
InstallMethod( SimpleAlgebra, 
                "for semisimple finite group algebras", 
                true, 
                [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
                0,
function( FqG, K, H, c )
local   G,      # Group
        Fq,     # Field
        n;      # Index of H in K
        

G := UnderlyingMagma( FqG );
Fq := LeftActingDomain( FqG );
n := Index( K, H );

if  Gcd( c[ 1 ], n ) = 1 and  c in CyclotomicClasses( Size( Fq ), n ) and 
                                    IsStronglyShodaPair(G, K, H ) then
    return SimpleAlgebraNC( FqG, K, H, c );
fi;

Error("The input is not appropriate!!!\n");

end);

#############################################################################
##
#  The function SimpleAlgebraNC computes the simple algebras 
##  FqG*e( G, K, H, c), for ( H, K ) a SSP of G and c a cyclotomic class 
##  of q=|Fq| module n=[K:H] containing generators of K/H.
## This version  does not check the inputs
##  
#O SimpleAlgebraNC( FqG, K, H, c )
##
InstallMethod( SimpleAlgebraNC, 
                "for semisimple finite group algebras", 
                true, 
                [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
                0,
function( FqG, K, H, c )
local   G,          # Group
        N,          # Normalizer of H in G
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



G := UnderlyingMagma( FqG );
Fq := LeftActingDomain( FqG );
q := Size( Fq );

if G = H then
    return FullMatrixAlgebra( Fq, 1 );
fi;

N := Normalizer( G, H );
epi := NaturalHomomorphismByNormalSubgroup( N, H );
QNH := Image( epi, N );
QKH := Image( epi, K );
gq := MinimalGeneratingSet( QKH )[ 1 ];
C1 := Set( List( c, i -> gq^i ) );
St := Stabilizer( QNH, C1, OnSets );
E := PreImage( epi, St );
ord := q^( Size( c )/Index( E, K ) );
if ord <= 2^16 then
    F := GF(ord);
else
    factors := FactorsInt(q);
    p:=factors[1];
    o:=Size(factors);
    F := GF( p, ConwayPolynomial( p, o*ord ) );
    # If q^o is too big then gap never finish to compute the ConwayPolynomial
fi;

return FullMatrixAlgebra( F, Index( G, K ) );

end);

#############################################################################
##
## The function SimpleAlgebraInfo compute the data describing simple algebras
##  QG*e( G, K, H ), for ( H, K ) a SSP of G, but first verify the inputs 
##
#O SimpleAlgebraInfo( QG, K, H ) 
##
InstallMethod( SimpleAlgebraInfo, 
                "for semisimple rational group algebras", 
                true, 
                [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
                0,
function( QG, K, H )
local   G;      # Group

G := UnderlyingMagma( QG );

if  IsStronglyShodaPair(G, K, H ) then
    return SimpleAlgebraInfoNC( QG, K, H );
fi;

Error("The input is not appropriate!!!\n");

end);

#############################################################################
##
## The function SimpleAlgebraInfoNC compute the data describing simple algebras
##  QG*e( G, K, H ), for ( H, K ) a SSP of G 
##
#O SimpleAlgebrasInfoNC( QG, K, H ) 
##
InstallMethod( SimpleAlgebraInfoNC, 
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
        Epi,        # N--->N/H
        Epi2,       # NH--->NH/KH
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
        i,ll;     # Controlers
        
    G := UnderlyingMagma( QG );
    if G = H then
        return [1,1,[],[]];
    fi;
    
    #First one computes an idependent set PrimGen of generators of a Primary
    #decomposition of N/K
    N   := Normalizer(G,H);
    Epi := NaturalHomomorphismByNormalSubgroup( N, H ) ;
    NH  := Image(Epi,N);
    KH  := Image(Epi,K);
    k   := Product(IndependentGeneratorsOfAbelianGroup(KH));
    ok  := Order(k);
    Potk:= List([1..ok],x->k^x);
    if N=K then
        return [Size(G)/Size(N),ok,[ ],[ ]];
    else
        Epi2:=NaturalHomomorphismByNormalSubgroup( NH, KH ) ;
        NdK:=Image(Epi2,NH);
        PrimGen:=IndependentGeneratorsOfAbelianGroup(NdK);
        # Using PrimGen one computes an independent set Gen of
        # generators of an invariant decomposition of N/K
        l := Length(PrimGen);
        o := List( [ 1 .. l ], i -> Order(PrimGen[i]) );
        p := List( [ 1 .. l ], i -> FactorsInt(o[i])[1] );
        primes:=Unique(p);
        lp:= Length(primes);
        first:=List( [ 1 .. lp ], i -> Position(p,primes[i]) );
        g := Product(List(first,i->PrimGen[i]));
        Gen:=[g];
        ll:=lp;
        plus:=0;
        while ll<l do
            next:=[];
            for i in [1..lp] do
                newpos:=Position(p,primes[i],first[i]+plus);
                if newpos <> fail then
                    Add(next,newpos);
                fi;
            od;
            g:=Product(List(next,i->PrimGen[i]));
            Add(Gen,g);
            ll:=ll+Length(next);
            plus:=plus+1;
        od;
        gen:=List([1..Length(Gen)], i->PreImagesRepresentative(Epi2,Gen[i]));
        return [Size(G)/Size(N),ok, List([1..Length(Gen)],
                i->[Order(Gen[i]),
                RemInt(Position(Potk,k^gen[i]),ok),
                RemInt(Position(Potk,gen[i]^Order(Gen[i])),ok)]),
                List([1..Length(Gen)-1],i->List([i+1..Length(Gen)],
                j->RemInt(Position(Potk,Comm(gen[j],gen[i])),ok)))];
    fi;

end);

#############################################################################
##
##  The function SimpleAlgebraInfo compute the data describing simple algebra 
##  FqG*e( G, K, H, c) for ( H, K ) a SSP of G and c a cyclotomic class 
##  of q=|Fq| module n=[K:H], containing generators of K/H, 
##  but first verify the inputs 
##
#O SimpleAlgebraInfo( FqG, K, H, c )  
##
InstallMethod( SimpleAlgebraInfo, 
                "for semisimple finite group algebras", 
                true, 
                [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
                0,
function( FqG, K, H, c )  
local   G,      # Group
        Fq,     # Field
        n;      # Index of H in K
        

G := UnderlyingMagma( FqG );
Fq := LeftActingDomain( FqG );
n := Index( K, H );

if  Gcd( c[ 1 ], n ) = 1 and  c in CyclotomicClasses( Size( Fq ), n ) and 
                                    IsStronglyShodaPair(G, K, H ) then
    return SimpleAlgebraInfoNC( FqG, K, H, c );
fi;

Error("The input is not appropriate!!!\n");


end);

#############################################################################
##
##  The function SimpleAlgebraInfoNC compute the data describing simple algebra 
##  FqG*e( G, K, H, c) for ( H, K ) a SSP of G and c a cyclotomic class 
##  of q=|Fq| module n=[K:H], containing generators of K/H.
##
#O SimpleAlgebraInfoNC( FqG, K, H, c )  
##
InstallMethod( SimpleAlgebraInfoNC, 
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

G := UnderlyingMagma( FqG );
Fq := LeftActingDomain( FqG );
q := Size( Fq );

if G = H then
return [ 1, q ];
fi;

N := Normalizer( G, H );
epi := NaturalHomomorphismByNormalSubgroup( N, H );
QNH := Image( epi, N );
QKH := Image( epi, K );
gq := MinimalGeneratingSet( QKH )[ 1 ];
C1 := Set( List( c, ii -> gq^ii ) );
St := Stabilizer( QNH, C1, OnSets );
E := PreImage( epi, St );
ord := q^( Size( c )/Index( E, K ) );

return [ Index( G, K ), ord ];

end);

#############################################################################
##
##  The function StronglyShodaPairs computes a list of strongly Shoda pairs 
##  of the group G that covers the complete set of primitive central idempotents
##  of the rational group algebra QG realizable by strongly Shoda pairs
##
#A  StronglyShodaPairs( G )
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
## The function StronglyShodaPairsAndIdempotents computes a record (SSPs, PCIs),
## where SSPs is a list of SSP that covers the complete set of primitive 
## central idempotents, PCIs, of the rational group algebra QG 
##
#A StronglyShodaPairsAndIdempotents( QG )
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

G := UnderlyingMagma(QG);

if Tester( StronglyShodaPairs )( G ) then


        return rec( StronglyShodaPairs := StronglyShodaPairs( G ), 
        PrimitiveCentralIdempotents := List( StronglyShodaPairs( G ), i -> eGKH( QG, i[1], i[2] ) )); 
else

CCS:=ConjugacyClassesSubgroups(G);
LCCS:=Length(CCS);
DG:=DerivedSubgroup(G); 
KHs:=[];
eGKHs:=[];
SeGKHs:=Zero(QG);
if Size(G)=1 then
    return [[G,G], [One(QG)]];
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

if SeGKHs<>One(QG) then 
    Print(  "Caution! Some primitive central idempotents are not realizable ", 
            "by strongly Shoda pairs!!!\n");
fi;

Setter( StronglyShodaPairs )( G , KHs ); 

return rec( StronglyShodaPairs := KHs, PrimitiveCentralIdempotents := eGKHs );

fi;

end);

#############################################################################
## 
## The following function search an element K such that (K,H) is a SSP
## and returns [ [ K, H ], e( G, K, H ) ] or returns fail, if
## such K doesn't exist
##
#F SearchingKForSSP(QG,H)
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

        G:=UnderlyingMagma(QG);
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
## The following function computes e(G,K,H)    
## Note that actually it returns a list of the form [ [K,H], eGKH ]
##
#F eG( QG, K, H )
##
InstallGlobalFunction(eG, function(QG,K,H)
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
        nRTNdK; # Cardinal of RTNdK

        Eps:=Epsilon(QG,K,H);
        G:=UnderlyingMagma(QG);
        NH:=Normalizer(G,H);
        if NH=G then
            return [ [ K, H ], Eps ];
        else
            NdK:=Normalizer(G,K);
            RTNH:=RightTransversal(NdK,NH);
            eGKH1:=Sum(List(RTNH,g->Conjugate(QG,Eps,g)));
            eGKH:=eGKH1;
            if NdK<>G then
                RTNdK:=RightTransversal(G,NdK); 
                nRTNdK:=Length(RTNdK);  
                for i in [ 2 .. nRTNdK ] do
                    g:=RTNdK[i];
                    eGKH1g:=Conjugate(QG,eGKH1,g);
                    if eGKH1*eGKH1g <> Zero(QG) then    
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
## The function StronglyShodaPairsAndIdempotents computes the record [SSPs, PCIs], 
## where SSPs is a list of the SSP and cyclotomic classes that covers 
## the complete set of primitive central idempotents,PCIs, 
## of the finite group algebra FqG
##
#A StronglyShodaPairsAndIdempotents( FqG )
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
G := UnderlyingMagma( FqG  );
Fq := LeftActingDomain( FqG );
F := FamilyObj(Zero(FqG));
elmsG := Elements(G);
q := Size( Fq );
zero := Zero(Fq);
e := [Hat(FqG,G)];
SSPsG := StronglyShodaPairs(G);
list := [ [ SSPsG[ 1 ][1], SSPsG[ 1 ][2], [ [ 0 ] ] ] ];
setind := [];
lltrace := [];
lcc := [];
lorders := [];
lprimitives := [];
for p in [2..Size(SSPsG)] do
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
        # If q^o is too big then gap never finish to compute the ConwayPolynomial
        # Is there a method to avoid this calculation if this is the case?
        #if o > 67 then ###### !!! for q=2 ( for q=3 is 37 and for q=5 is 31 ...)
        #    Print("fallo","\n");
        #    return fail;
        #else
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
        #fi;
    fi;
    etemp := [];
    templist := [];
    for i in cc do
        if Gcd(i[1],n)=1 then
            idemp := eGKHc(FqG, K, H, i, ltrace);
            if not(idemp in etemp) then
                Add(etemp, idemp);
                Add( templist, i );
            fi;
        fi;
    od;
    Append( e, etemp );
    Add( list, [ K, H, templist ] );
od;
return rec( StronglyShodaPairs := list, PrimitiveCentralIdempotents := e );
end);



#############################################################################
##
##  The function PrimitiveCentralIdempotentsFromSSP computes the set of primitive central 
##  idempotents of the group algebra FG, realizable by strongly Shoda pairs, where
##  FG is either a ratioanl group algebra or finite group algebra
##
#O  PrimitiveCentralIdempotentsFromSSP( FG )
##
InstallGlobalFunction( PrimitiveCentralIdempotentsFromSSP, function( FG )

return StronglyShodaPairsAndIdempotents( FG ).PrimitiveCentralIdempotents; 

end);