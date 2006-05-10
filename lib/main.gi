F#############################################################################
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
##                                                                         ##
##                   WEDDERBURN DECOMPOSITION                              ##
##                                                                         ##
#############################################################################


#############################################################################
##
#O WedderburnDecomposition( FG )
##
## The function WeddDecomp computes the Wedderburn components realizable by
## strongly Shoda pairs of the underlying group, of the semisimple group algebra 
## FG over the finite field or the field of rationals as matrix algebras over 
## cyclotomic algebras and stores the result as an attribute of FG. 
## It uses the attribute WeddDecomp and IsStronglyMonomial to display a warning.
## The reason for such combination of operation 'WedderburnDecomposition' and
## attribute 'WeddDecomp' was in the necessity of displaying the warning each
## time when we refer to this information
##
InstallMethod( WedderburnDecomposition, 
    "for semisimple zero characteristic or finite group algebra", 
    true, 
    [ IsGroupRing ], 
    0,
function( FG )
local G;   # Underlying group
G := UnderlyingMagma( FG );
if not IsStronglyMonomial( G ) and IsSemisimpleFiniteGroupAlgebra( FG ) then 
    Print("Wedderga: Warning!!!\nThe direct product of the output is a PROPER direct factor of the input! \n");
fi;
return WeddDecomp( FG );
end);


#############################################################################
##
#A WeddDecomp( FG )
##
## The function WeddDecomp computes the Wedderburn components realizable by
## strongly Shoda pairs of the underlying group, of the semisimple group algebra 
## FG over the finite field or the field of rationals as matrix algebras over 
## cyclotomic algebras and stores the result as an attribute of FG. 
## This is an auxiliar function not to be documented.
##
InstallMethod( WeddDecomp, 
    "for semisimple zero characteristic or finite group algebra", 
    true, 
    [ IsGroupRing ], 
    0,
function( FG )
local   A,      # Simple algebra
        x,      # description of current component
        i,      # Counter
        output;
output := [];

if IsSemisimpleFiniteGroupAlgebra( FG ) then
  for x in StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs do
    A := SimpleAlgebraByStronglySPNC( FG, x[ 1 ], x[ 2 ], x[ 3 ][ 1 ]);
    Append(output, List(x[3], i -> A ) );
  od;
  return output;
  
elif IsZeroCharacteristicGroupAlgebra( FG ) then
  for x in GenWeddDecomp( FG ) do
    A := SimpleAlgebraByData(x);
    Add(output, A );
  od;
  return output;
  
else
  Error("Wedderga: <FG> must be a semisimple group algebra over rationals or over finite field!!!");
fi;  

end);


#############################################################################
##
#O WedderburnDecompositionInfo( FG ) 
##
## The function WeddDecompInfo compute a list of numerical data describing 
## the Wedderburn components, realizable by strongly Shoda pairs of the 
## underlying group, of the semisimple group algebra FG over a finite field or 
## the field of rationals and stores the result as an attribute of FG. 
## It uses the attribute WeddDecomp and IsStronglyMonomial to display a warning
##
InstallMethod( WedderburnDecompositionInfo , 
    "for semisimple zero-characteristic or finite group algebra", 
    true, 
    [ IsGroupRing ], 
    0,
function( FG )
local   G;      # Underlying group

G := UnderlyingMagma( FG );

if IsSemisimpleFiniteGroupAlgebra( FG ) and not IsStronglyMonomial(G) then 
    Print("Wedderga: Warning!!!\nThe direct product of the output is a PROPER direct factor of the input! \n");
fi;

return WeddDecompInfo( FG );

end);


#############################################################################
##
#A WeddDecompInfo( FG ) 
##
## The function WeddDecompInfo compute a list of numerical data describing 
## the Wedderburn components, realizable by strongly Shoda pairs of the 
## underlying group, of the semisimple group algebra FG over a finite field or 
## the field of rationals and stores the result as an attribute of FG. 
## This is an auxiliar function not to be documented.
##
InstallMethod( WeddDecompInfo , 
    "for semisimple zero-characteristic or finite group algebra",
    true, 
    [ IsGroupRing ], 
    0,
function( FG )
local   G,      # Group
        pairs,  # Strongly Shoda pairs of G
        A,      # Simple algebra
        i,      # Counter
        exp,    # Exponent of G
        br,     # List of lists of strongly Shoda triples
        sst,    # an element of sst
        chi,    # an irreducible character
        cf,     # character field of chi
        output;
        
G := UnderlyingMagma(FG);
output := [];

if IsSemisimpleRationalGroupAlgebra( FG ) then

    if IsAbelian(G) then 
      return List( RationalClasses(G), x -> [ 1, Order(Representative(x)) ] );
    else
      if HasStronglyShodaPairs( G ) then
          pairs := StronglyShodaPairs( G );
      else
          pairs := StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs;
      fi;
      
      for i in pairs do
        Add(output, SimpleAlgebraByStronglySPInfoNC( FG, i[ 1 ], i[ 2 ] ) );
      od;
  
      if not IsStronglyMonomial(G) then 
        exp := Exponent(G);
        br:=BWNoStMon(G);
      
        for sst in br do
          chi:=sst[1];
          cf:=sst[2];
          if Length(sst)=2 then 
            Add(output,[chi[1],Conductor(cf)]);
          else
            Add(output,SimpleAlgebraByStronglySTInfo(exp,chi[1],cf,sst[4],sst[3]));
          fi;
        od;
      fi;
    fi;
  
  return output;

elif IsZeroCharacteristicGroupAlgebra(FG) then
    
    for i in GenWeddDecomp(FG) do
      A := SimpleAlgebraInfoByData(i);
      Append( output, [ A ] );
    od;  
    return output;
    
elif IsSemisimpleFiniteGroupAlgebra( FG ) then

    for i in StronglyShodaPairsAndIdempotents( FG ).StronglyShodaPairs do
        A := SimpleAlgebraByStronglySPInfoNC( FG, i[ 1 ], i[ 2 ], i[ 3 ][ 1 ]);
        Append(output, List(i[3], j -> A ) );
    od;
    return output;
    
else

    Error("Wedderga: <FG> must be a zero-characteristic or finite semisimple group algebra!!!");

fi;

end); 


#############################################################################
##
#O GenWeddDecomp( KG )
## 
InstallMethod( GenWeddDecomp,
"for semisimple infinite group algebras",
true,
[ IsGroupAlgebra ],
0,
function(KG)
local 
K,          # Coefficient Field
G,          # Underlying Group
wdd,        # Wedderburn Decomposition data for QG
output,     # the output
x,          # an element of wdd
z,          # Centre of a Wedderburn component of QG
F,          # Centre of a Wedderburn component of KG
a,          # The number of Wedderburn components of KG associated to a simple 
            # component of QG
i,          # counter
n,          # Matrix size of a Wederburn component of QG
ok,         # Order of root of unity
Gal,        # The group of a crossed product of a simple component of QG
coc,        # The cocycle of a crossed product of a simple component of QG and KG
Fxi,        # F(ok)
d,          # Factor of increase of matrix size
condK,      # Conductor of K
m,          # Lcm(condK,ok)
redmok,     # Reduction Z_m --> Z_ok
redmcondK,  # Reduction Z_m --> Z_ok
gal;        # The group of a crossed product of a simple component of KG

K := LeftActingDomain(KG);
G := UnderlyingMagma(KG);
wdd := WeddDecompData(G);
if K=Rationals then
    return wdd;
else

    output := [];
    
    for x in wdd do
      n := x[1];
      z := x[2];
      if Length(x) = 2 then 
        z := x[2];
        F :=  Field(Union(GeneratorsOfField(z),GeneratorsOfField(K)));
        a := Dimension(z)*Dimension(K)/Dimension(F);
        for i in [1..a] do
          Add(output,[n,F]);
        od;
      else 
        ok := x[3];
        Gal := x[4];
        coc := x[5];
        F :=  Field(Union(GeneratorsOfField(z),GeneratorsOfField(K)));
        a := Dimension(z)*Dimension(K)/Dimension(F);
        Fxi := Field(F,[E(ok)]);
        d := Dimension(Field(z,[E(ok)]))/Dimension(Fxi);
        condK := Conductor(K);
        m := Lcm(condK,ok);
        redmok := ReductionModnZ(m,ok);
        redmcondK := ReductionModnZ(m,condK);
        gal := Subgroup(Units(ZmodnZ(ok)),
                    Filtered(Gal,y->
                            Size(
                                Intersection(
                                    GaloisStabilizer(K),
                                    List(PreImages(redmok,y),w->Int(w^redmcondK))
                                            )
                                )<>0
                            )
                        );
        for i in [1..a] do
            Add(output,[n*d,F,ok,gal,coc]);
        od;
      fi;
    od;
fi;

return output;
end);


#############################################################################
##
#A WeddDecompData( G )
##
InstallMethod( WeddDecompData,
"for numerical data for decomposition of semisimple infinite group algebras",
true,
[ IsGroup ],
0,
function(G)

local output,exp,br,sst,chi,cf;

output :=  List(StronglyShodaPairs(G),x->AddCrossedProductBySSP(G,x[1],x[2]));

if not IsStronglyMonomial(G) then
  exp := Exponent(G);
  br:=BWNoStMon(G);

  for sst in br do
    chi:=sst[1];
    cf:=sst[2];
    if Length(sst)=2 then 
      Add(output,[chi[1],cf]);
    else
      Add(output,AddCrossedProductBySST(exp,chi[1],cf,sst[4],sst[3]));
    fi;
  od;  
fi;

return output;

end);


#############################################################################
##
#O AddCrossedProductBySST( exp, n, cf , Gal , LSST )
##
InstallMethod(AddCrossedProductBySST,
"for semisimple infinite group algebras",
true,
[ IsInt, IsInt, IsCyclotomicField, IsGroup, IsList ],
0,
function( exp, n, cf , Gal , LSST )
local
  Galnum,       # Numeric version of Gal(Q(exp)/cf)
  pp,           # Maximum prime power divisors of n = Degree of the character
  LC,           # List of cocycles
  x,            # An element of LSST or of LC
  primes,       # List of primes covered by x, an SST
  a,            # The products of the elements of pp corresponding to primes
  Cond,         # The conductor of the coefficient field of the output
  GalCond,      # The reduction of Galnum module Cond (the grading group)
  coc,          # The cocycle of the algebra
  out,          # The output of definition of coc
  redu;         # Reduction Cond to the conductor corresponding to a partial cocycle
  

  Galnum := Image(GalToInt(Gal));
  
  if Gcd(exp,4)=2 then
    Galnum := Subgroup(Units(ZmodnZ(exp)),PreImage(ReductionModnZ(exp,exp/2),Galnum));
  fi;
  
  pp := PrimePowersInt(n);
  
  LC:=[];
  for x in LSST do
    primes := x[4];
    a:= Product(primes,p->p^pp[Position(pp,p)+1]);
    Add(LC,CocycleByData(exp,Galnum,cf,x[1],x[2],x[3],a));
  od;
  
  Cond := Lcm(List(LC,x->x[1]));
  GalCond := Subgroup(Units(ZmodnZ(Cond)),Image(ReductionModnZ(exp,Cond),Galnum));
  coc := function(a,b)
    local out,x,redu;  
      out := Zero(ZmodnZ(Cond));
      for x in LC do
        redu := ReductionModnZ(Cond,x[1]);
        out := out + (Cond/x[1])*ZmodnZObj(Int(x[2](a^redu,b^redu)),Cond);
      od;
      return out;
  end;

if Size(GalCond)=1 then 
    return [n,cf];
else
    return [ n/Size(GalCond), cf, Cond, GalCond, coc ];
fi;
   
end);



#############################################################################
##                                                                         ##
##                       SIMPLE ALGEBRA                                    ##
##                                                                         ##
#############################################################################


#############################################################################
##
#O SimpleAlgebraByStronglySP( QG, K, H )
##
## The function SimpleAlgebraByStronglySP computes the simple algebras 
## QG*e( G, K, H) if ( K, H ) is a SSP of G 
## This version does not check the input
##
InstallOtherMethod( SimpleAlgebraByStronglySP, 
    "for semisimple rational group algebras", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function( QG, K, H)
if IsStronglyShodaPair( UnderlyingMagma( QG ), K, H ) then
    return SimpleAlgebraByStronglySPNC( QG, K, H );
else
    Error("Wedderga: <(K,H)> should be a strongly Shoda pair of the underlying group of <QG>\n");
fi;
end);


#############################################################################
##
#O SimpleAlgebraByStronglySP( FqG, K, H, C ) 
##
## The function SimpleAlgebraByStronglySP verifies if ( H, K ) is a SSP of G and
## C is a cyclotomic class of q=|Fq| module n=[K:H] containing generators
## of K/H, and in that case computes the simple algebra  FqG*e( G, K, H, C)
##
InstallMethod( SimpleAlgebraByStronglySP, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
    0,
function( FqG, K, H, C )
local   G,      # Group
        n,      # Index of H in K
        j,      # Integer
        q,      # Size of Fq
        C1;     # Cyclotomic Class

G := UnderlyingMagma( FqG );
q := Size( LeftActingDomain( FqG ) );
n := Index( K, H );

if not(IsStronglyShodaPair(G, K, H )) then
    Error("Wedderga: (<K>,<H>) should be a strongly Shoda pair of the underlying group of <FqG>\n");
elif IsCyclotomicClass( q, n, C) and Gcd(n,C[1]) =1 then
    return SimpleAlgebraByStronglySPNC( FqG, K, H, C );
else Error("Wedderga: <C> should be a generating cyclotomic class module the index of <H> in <K>\n");
fi;

end);

#############################################################################
##
#O SimpleAlgebraByStronglySP( FqG, K, H, c ) 
##
## The function SimpleAlgebraByStronglySP verifies if ( H, K ) is a SSP of G and
## c is an integer coprime with n=[K:H]. 
## If the answer is positive then returns SimpleAlgebraByStronglySP(FqG, K, H, C) where
## C is the cyclotomic class of q=|Fq| module n=[K:H] containing c.
##
InstallOtherMethod( SimpleAlgebraByStronglySP, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsPosInt ], 
    0,
function( FqG, K, H, c )
local   G,      # Group
        n;      # Index of H in K        
        
G := UnderlyingMagma( FqG );
n := Index( K, H );

if  IsStronglyShodaPair(G, K, H ) then
  if Gcd( c, n ) = 1 then
    return SimpleAlgebraByStronglySPNC( FqG, K, H, c  mod n);
  else
    Error("Wedderga: <c> should be coprime with the index of <H> in <K>");   
  fi;
else
   Error("Wedderga: (<K>,<H>) should be a strongly Shoda pair of the underlying group of <FqG>\n");
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
InstallOtherMethod( SimpleAlgebraByStronglySPNC, 
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
        
G := UnderlyingMagma( QG );
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
             local x;
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
#O AddCrossedProductBySSP( G, K, H )
##
InstallMethod( AddCrossedProductBySSP,
"for semisimple infinite group algebras",
true,
[ IsGroup, IsGroup, IsGroup ],
0,
function( G, K, H )
local   N,          # Normalizer of H in G
        ind,        # Index of N in G
        ok,         # Order of k
        Epi,        # N --> N/H
        NH,         # NH/H
        KH,         # K/H
        k,          # Generator of K/H
        Epi2,       # NH --> NH/KH        
        NdK,        # N/K
        bij,bijunit,
        coc,        # Twisting for the crossed product over NdK
        Uok,        # Units(ZmodnZ(ok))
        funNdK,     # Embedding of NdK in Uok,
        GalSSP,     # Subgroup(Uok,Image(funNdK))
        cocSSP,     # cocycle in Z^2(GalSSP,<E(ok)>)
        chi,        # Monomial character of G induced the SSP (K,H)
        cf;         # Fields of character values of chi = Centre
  
N   := Normalizer(G,H);
ind := Index(G,N);
ok := Index( K, H );
if N=K then

    return [ ind, CF(ok) ];

else # if N_G(H) <> K

    Epi := NaturalHomomorphismByNormalSubgroup( N, H ) ;
    NH  := Image(Epi,N);
    KH  := Image(Epi,K);
    repeat
        k  := Random(KH);
    until Order(k) = ok;
    Epi2:=NaturalHomomorphismByNormalSubgroup( NH, KH ) ;
    NdK:=Image(Epi2,NH);
    bij := MappingByFunction(ZmodnZ(ok),KH,i->k^Int(i));
    
    # The cocycle in Z^2(NdK,<E(ok)>)
    coc := function(a,b)
       return PreImagesRepresentative(bij,
                              PreImagesRepresentative(Epi2,a*b)^-1 *
                              PreImagesRepresentative(Epi2,a) *
                              PreImagesRepresentative(Epi2,b) );
       end;   

    # The cocycle in Z^2(GalSSP,<E(ok)>)
    Uok:=Units(ZmodnZ(ok));
    bijunit := MappingByFunction(Uok,KH,i->k^Int(i));

    funNdK := MappingByFunction(NdK,Uok,
        function(n) 
            return PreImagesRepresentative(bijunit,
                                 k^PreImagesRepresentative( Epi2 , n ) );
              end
              );
    GalSSP := Subgroup(Uok,Image(funNdK));
    cocSSP := function(a,b)
                return 
        coc(PreImagesRepresentative(funNdK,a),PreImagesRepresentative(funNdK,b));
                end;
    
    chi := LinCharByStronglySP(K,H)^G;
    cf := Field( chi );
                    
    return [ ind, cf, ok , GalSSP , cocSSP ];
fi;      
end);



#############################################################################
##
#O SimpleAlgebraByData( x )
##
## x is a list of length 5.
##
InstallMethod( SimpleAlgebraByData,
"for semisimple infinite group algebras",
true,
[ IsList ],
0,
function(x)

local 
L,     #Field 
cond,  #
redu,  #The reduction 
act,   #The action
coc;   #The cocycle

if Length(x) = 2 or Size(x[4])=1 then
    if x[1] = 1 then 
        return x[2];
    else
        return FullMatrixAlgebra( x[2], x[1] );
    fi;
else
    L := Field(x[2],[E(x[3])]);
    cond := Lcm(Conductor(L),x[3]);
    redu := ReductionModnZ(cond,x[3]);
    
    act := function(a) 
             return ANFAutomorphism(CF(cond),Int(PreImagesRepresentative(redu,a)));
             end;
             
    coc := function(a,b)
            return E(x[3])^x[5](Int(a),Int(b));
            end;
    
    if x[1] = 1 then 
        return CrossedProduct(L,x[4],act,coc);
    else
        return FullMatrixAlgebra( CrossedProduct(L,x[4],act,coc) , x[1]/Size(x[4]) );
    fi;
fi;

end);


#############################################################################
##
#O SimpleAlgebraByCharacter( FG, chi )
##
InstallMethod( SimpleAlgebraByCharacter,
"for semisimple infinite group algebras",
true,
[ IsGroupRing, IsCharacter ],
0,
function( FG, chi )

end);


#############################################################################
##
#O SimpleAlgebraByStronglySPNC( FqG, K, H, C )
##
## The function SimpleAlgebraByStronglySPNC computes simple algebras 
## FqG*e( G, K, H, C), for ( H, K ) a SSP of G and C a cyclotomic class 
## of q=|Fq| module n=[K:H] containing generators of K/H.
## This version does not check the input
##
InstallMethod( SimpleAlgebraByStronglySPNC, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
    0,
function( FqG, K, H, C )
local   G,          # Group
        Fq,F,       # Fields
        q,          # Order of Fq
        N,          # Normalizer of H in G
        epi,        # N -->N/H
        QNH,        # N/H
        QKH,        # K/H
        gq,         # Generator of K/H
        C1,         # Cyclotomic class of q module [K:H] in N/H
        St,         # Stabilizer of C1 in N/H
        E,          # Stabilizer of C1 in G
        ord,        # Integer
        factors,    # prime factors of q
        p,          # The only prime divisor of q
        o,          # q = p^o
        ind;        # index of K in G        
        
G := UnderlyingMagma( FqG );
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
C1 := Set( List( C, i -> gq^i ) );
St := Stabilizer( QNH, C1, OnSets );
E := PreImage( epi, St );
ord := Size( C )/Index( E, K ) ;

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
#O SimpleAlgebraByStronglySPNC( FqG, K, H, c ) 
##
## The function SimpleAlgebraByStronglySP verifies if ( H, K ) is a SSP of G and
## c is an integer coprime with n=[K:H]. 
## In the answer is positive then return SimpleAlgebraByStronglySP(FqG, K, H, C) where
## C is the cyclotomic class of q=|Fq| module n=[K:H] containing c.
##
InstallOtherMethod( SimpleAlgebraByStronglySPNC, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsPosInt ], 
    0,
function( FqG, K, H, c )
local   G,      # Group
        n,      # Index of H in K
        q,      # Size of Fq
        j,      # integer module n
        C;      # q-cyclotomic class module [K,H] containing c

G := UnderlyingMagma( FqG );
n := Index( K, H );
q:=Size( LeftActingDomain( FqG ) );
C := [ c mod n];
j:=q*c mod n;
while j <> C[1] do
  Add( C, j );
  j:=j*q mod n;
od;  
    return SimpleAlgebraByStronglySPNC( FqG, K, H, C );

end);


#############################################################################
##
#O SimpleAlgebraByStronglySPInfo( QG, K, H ) 
##
## The function SimpleAlgebraByStronglySPInfo compute the data describing simple algebras
## QG*e( G, K, H ), for ( H, K ) a SSP of G, but first verify the inputs 
##
InstallOtherMethod( SimpleAlgebraByStronglySPInfo, 
    "for semisimple rational group algebras", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function( QG, K, H )
if  IsStronglyShodaPair( UnderlyingMagma( QG ), K, H ) then
    return SimpleAlgebraByStronglySPInfoNC( QG, K, H );
else
    Error("Wedderga: <(K,H)> should be a strongly Shoda pair of the underlying group of <QG>\n");
fi;
end);


#############################################################################
##
#O SimpleAlgebraInfoByData( x )
##
## x is a list of length 5.
##
InstallMethod( SimpleAlgebraInfoByData,
"for semisimple infinite group algebras",
true,
[ IsList ],
0,
function(x)
local 
Cond,           # Positive integer 
coc,            # cocycle
PrimGen,        # Indenpendent Generators of GalCond  
l, o, p,        # Positive integer and lists of integers
primes2,        # Duplicate of p
lp,             # Length of primes2
first,          # Positions,
g,              # One generator
Gen,            # Generators of GalCond
ll, plus, next, 
i, j, newpos,   # Counters
genF,           # E(Cond)
powgenF,        # Powers of genF
beta,           # numerical value of cyclic cocycle 
h,              # Group element
c;              # Value of cocycle 

if Length(x) = 2 then 
    return x;
elif Size(x[4])=1 then
    return [x[1],x[2],x[3]];
else
    Cond := x[3];
    coc := x[5];
    
# Computing a set Gen of generators of the canonical decomposition of x[4]
# an abelian group

    PrimGen:=IndependentGeneratorsOfAbelianGroup(x[4]);
    l := Length( PrimGen );
    o := List( [ 1 .. l ], i -> Order( PrimGen[i] ) );
    p := List( [ 1 .. l ], i -> FactorsInt( o[i] )[1] );
    primes2 := DuplicateFreeList( p );
    lp:= Length( primes2 );
    first:=List( [ 1 .. lp ], i -> Position( p, primes2[i] ) );
    g := Product( List( first, i -> PrimGen[i] ) );
    Gen:=[ g ];
    ll:=lp;
    plus:=0;
    while ll<l do
        next:=[];
        for i in [ 1 .. lp ] do
            newpos := Position( p, primes2[i], first[i]+plus );
            if newpos <> fail then
                Add( next, newpos );
            fi;
        od;
        g:=Product( List( next, i -> PrimGen[i] ) );
        Add( Gen, g );
        ll:=ll+Length(next);
        plus:=plus+1;
    od;
     
    o:=List(Gen,x->Order(x));
    beta := [];
    for i in [1..Length(Gen)] do
        g:=Gen[i];
        h:=g;
        c:=Zero(ZmodnZ(Cond));
        for j in [1..o[i]-1] do
            c:=c+coc(g,h);
            h:=h*g;
        od;
        Add(beta, Int(c));
    od;
         
    if Size(Gen)=1 then
        return [x[1]/Size(x[4]),x[2],Cond, [o[1], Int(Gen[1]) , beta[1]]]; 
                   
    else
        return [x[1]/Size(x[4]),x[2],Cond,
            List([1..Length(Gen)],i->
                [o[i], Int(Gen[i]) , beta[i]]),
            List( [1..Length(Gen)-1], i -> 
                       List( [i+1..Length(Gen)], 
                       j -> 
                            Int(coc(Gen[j],Gen[i])-coc(Gen[i],Gen[j]))
                            )
                )
                ];        
    fi;       

fi;

end);


#############################################################################
##
#O SimpleAlgebraInfoByCharacter( FG, chi )
##
InstallMethod( SimpleAlgebraInfoByCharacter,
"for semisimple infinite group algebras",
true,
[ IsGroupRing, IsCharacter ],
0,
function( FG, chi )

end);


#############################################################################
##
#O SimpleAlgebraByStronglySPInfo( FqG, K, H, C )
##
## The function SimpleAlgebraByStronglySPInfo cheks that (K,H) is a strongly 
## Shoda pair of G, the underlying group of the semisimple finite group algebra
## FqG with coefficients in the field of order q and if C is a generating 
## q-cyclotomic class module n=[K:H]. In that case computes the data describing 
## the simple algebra FqG*e( G, K, H, C)
##

InstallMethod( SimpleAlgebraByStronglySPInfo, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
    0,
function( FqG, K, H, C )  
local   G,      # Group
        C1,     # Cyclotomic class,
        j,      # integer
        q,      # Size of Fq
        n;      # Index of H in K

G := UnderlyingMagma( FqG );
q := Size( LeftActingDomain( FqG ) );
n := Index( K, H );

if not(IsStronglyShodaPair(G, K, H )) then
    Error("Wedderga: (<K>,<H>) should be a strongly Shoda pair of the underlying group of <FqG>\n");
elif IsCyclotomicClass( q, n, C) and Gcd(n,C[1]) =1 then
    return SimpleAlgebraByStronglySPInfoNC( FqG, K, H, C );
else Error("Wedderga: <C> should be a generating cyclotomic class module the index of <H> in <K>\n");
fi;

end);


#############################################################################
##
#O SimpleAlgebraByStronglySPInfo( FqG, K, H, c )
##
## The function SimpleAlgebraByStronglySPInfo cheks that (K,H) is a strongly 
## Shoda pair of G, the underlying group of the semisimple finite group algebra
## FqG with coefficients in the field of order q and in that c is a positive
## integer coprime with n=[K:H]. In that case computes the data describing the 
## simple algebra FqG*e( G, K, H, C) for C the q-cyclotomic class module n
## containing c
##
InstallOtherMethod( SimpleAlgebraByStronglySPInfo, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsPosInt ], 
    0,
function( FqG, K, H, c )  
local   G,      # Group
        n;      # Index of H in K

G := UnderlyingMagma( FqG );

if IsStronglyShodaPair(G, K, H ) then
  n := Index( K, H );
  if c<n and Gcd( c, n ) = 1 then
    Print("pasé");
    return SimpleAlgebraByStronglySPInfoNC( FqG, K, H, c );
  else 
    Error("Wedderga: <c> should be coprime with the index of <H> in <K>\n");
  fi;  
else
   Error("Wedderga: (<K>,<H>) should be a strongly Shoda pair of the underlying group of <FqG>\n");
fi;

end);


#############################################################################
##
#O SimpleAlgebraByStronglySPInfoNC( QG, K, H ) 
##
## The function SimpleAlgebraByStronglySPInfoNC compute the data describing simple 
## algebras QG*e( G, K, H ), for ( H, K ) a SSP of G 
##
InstallOtherMethod( SimpleAlgebraByStronglySPInfoNC, 
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
        
    G := UnderlyingMagma( QG );
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
#O SimpleAlgebraByStronglySPInfoNC( FqG, K, H, C )
##
## The function SimpleAlgebraByStronglySPInfo computes the data describing 
## the algebra FqG*e( G, K, H, C) without checking conditions on the input
##
InstallMethod( SimpleAlgebraByStronglySPInfoNC, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ], 
    0,
function( FqG, K, H, C )  
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
# We guarantee that QKH is cyclic so we can randomly obtain its generator
repeat
  gq := Random(QKH);
until Order(gq) = Size(QKH);
C1 := Set( List( C, ii -> gq^ii ) );
St := Stabilizer( QNH, C1, OnSets );
E := PreImage( epi, St );
ord := q^( Size( C )/Index( E, K ) );

return [ Index( G, K ), ord ];

end);


#############################################################################
##
#O SimpleAlgebraByStronglySPInfoNC( FqG, K, H, c )
##
## The function SimpleAlgebraByStronglySPInfo computes the data describing 
## the algebra FqG*e( G, K, H, C), where C is the q=|Fq|-cyclotomic class module
## [K:H] containing c, without checking conditions on the input
##
InstallOtherMethod( SimpleAlgebraByStronglySPInfoNC, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsPosInt ], 
    0,
function( FqG, K, H, c )  

local   G,      # Group
        n,      # Index of H in K
        q,      # Size of Fq
        j,      # integer module n
        C;      # q-cyclotomic class module [K,H] containing c

q := Size( LeftActingDomain( FqG ) );

G := UnderlyingMagma( FqG );
n := Index( K, H );
q:=Size( LeftActingDomain( FqG ) );
C := [ c ];
j:=q*c mod n;
while j <> c do
  Add( C, j );
  j:=j*q mod n;
od;  
    return SimpleAlgebraByStronglySPInfoNC( FqG, K, H, C );

end);


#############################################################################
##                                                                         ##
##            STRONGLY SHODA PAIRS AND IDEMPOTENTS                         ##
##                                                                         ##
#############################################################################


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
    [ IsGroup and IsFinite ], 
    0,
function( G )
local   QG;     # Rational Group Algebra
       
QG := GroupRing( Rationals, G ); 
        
return StronglyShodaPairsAndIdempotents(QG).StronglyShodaPairs;

end);


#############################################################################
##
#A StronglyShodaPairsAndIdempotents( QG )
##
## The attribute StronglyShodaPairsAndIdempotents of the rational group algebra QG 
## returns a record with components StronglyShodaPairs and PrimitiveCentralIdempotents where 
## StronglyShodaPairs = list of SSP that covers the complete set of primitive 
##       central idempotents of QG realizable by SSPs, 
## PrimitiveCentralIdempotents = list of PCIs of QG realizable by SSPs.
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

if HasStronglyShodaPairs( G ) then
    eGKHs := List( StronglyShodaPairs( G ), i -> 
                                     CentralElementBySubgroups( QG, i[1], i[2] ) );
    return rec( 
    StronglyShodaPairs := StronglyShodaPairs( G ), 
    PrimitiveCentralIdempotents := eGKHs); 

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
    if IsSubset( H, DG ) then
      if IsCyclic( FactorGroup( G, H ) ) then 
        idempeGKH:=eGsum(QG,G,H)[2]; 
        SeGKHs:= SeGKHs + idempeGKH;
        Add( KHs, [ G, H ] );
        Add( eGKHs, idempeGKH );
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

  SetStronglyShodaPairs( G , KHs ); 
  SetIsStronglyMonomial( G , SeGKHs=One(QG) );

  return rec( 
    StronglyShodaPairs := KHs, 
    PrimitiveCentralIdempotents := eGKHs);

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
G := UnderlyingMagma( FqG  );
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
            PrimitiveCentralIdempotents := e);
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
                    return eGsum(QG,K,H);
                else 
                    return fail;
                fi;
            else 
                X:=Difference(Elements(Cen),Elements(L));
                while X<>[] do
                    KH:=ClosureSubgroup( L, [X[1]] );
                    if IsCyclic(KH) and Centralizer(NHH,KH)=KH then
                        K:=PreImages(Epi,KH);
                        return eGsum(QG,K,H);
                    fi;
                    X:=Difference(X,KH);
                od;
            fi;
        fi;
    return fail;                  
    end);     


#############################################################################
##
#M eGsum( QG, K, H )
##
## The following function computes e(G,K,H)    
## Note that actually it returns a list of the form [ [K,H], eGKH ]
##
InstallMethod( eGsum,
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
G:=UnderlyingMagma(QG);
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

local G;

G := UnderlyingMagma( FG );
if not IsStronglyMonomial(G)  then 
   Print("Wedderga: Warning!!!\nThe output is a NON-COMPLETE list of prim. central idemp.s of the input! \n");
fi;

return StronglyShodaPairsAndIdempotents( FG ).PrimitiveCentralIdempotents; 
end);


#############################################################################
##
#E
##
