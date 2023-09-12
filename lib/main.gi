#############################################################################
##
#W  main.gi               The Wedderga package            Osnel Broche Cristo
#W                                                         Olexandr Konovalov
#W                                                            Aurora Olivieri
#W                                                           Gabriela Olteanu
#W                                                              Ángel del Río
#W                                                          Inneke Van Gelder
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
## The function WeddDecomp computes the Wedderburn components of the semisimple
## group algebra FG over a cyclotomic field F and for G an arbitrary 
## finite group, as matrix algebras over cyclotomic algebras and stores the 
## result as an attribute of FG. WedderburnDecomposition uses the attributes 
## WeddDecomp and IsCyclGroupAlgebra to display a warning.
## The reason for such combination of operation 'WedderburnDecomposition' and
## attribute 'WeddDecomp' was in the necessity of displaying the warning each
## time when we refer to this information.
##
InstallMethod( WedderburnDecomposition, 
    "for semisimple group algebra over cyclotomic fields", 
    true, 
    [ IsSemisimpleANFGroupAlgebra ], 
    0,
function( FG )

if not IsCyclGroupAlgebra( FG ) then  #IsCyclotomicAlgebra
    Print("Wedderga: Warning!!! \n", 
    "Some of the Wedderburn components displayed are FRACTIONAL MATRIX ALGEBRAS!!!\n\n");
fi;

Info( InfoWedderga, 2, "Info version : ", WedderburnDecompositionInfo( FG ) );

return WeddDecomp( FG );
end);


#############################################################################
##
#A WeddDecomp( FG )
##
## The function WeddDecomp computes the Wedderburn components of the semisimple 
## group algebra FG over a cyclotomic field F and for G an arbitrary 
## finite group, as matrix algebras over cyclotomic algebras and stores the 
## result as an attribute of FG. This is an auxiliar function not to be 
## documented.
InstallMethod( WeddDecomp, 
    "for semisimple group algebra over cyclotomic fields", 
    true, 
    [ IsSemisimpleANFGroupAlgebra ], 
    0,
function( FG )
local   A,      # Simple algebra
        x,      # description of current component
        output;
        
output := [];

if IsSemisimpleANFGroupAlgebra( FG ) then
 
  for x in GenWeddDecomp( FG ) do
    A := SimpleAlgebraByData(x);
    Add(output, A );
  od;
  
  return output;
  
else
  Error("Wedderga: <FG> must be a semisimple group algebra over a cyclotomic field!!!");
fi;  

end);

#############################################################################
##
#A WedderburnDecomposition( FG )
##
## The function WeddDecomp computes the Wedderburn components of the semisimple
## finite group algebra FG as matrix algebras over cyclotomic algebras and 
## stores the result as an attribute of FG. 
##
InstallMethod( WedderburnDecomposition, 
    "for semisimple finite group algebra", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra ], 
    0,
function(FG)

local G,      # Underlying group of FG
      F,      # Coefficient field of FG
      p,      # Characteristic of the field F
      m,      # Power of p in the size of the field F
      irr,    # Irreducible characters of G
      lex,    # lexicographical ordering function
      data,   #  list of 2-tuples, the first is the degree of the character x and 
              #  the second is the lcm(m, the power of p in the field where x can 
              #  be realized
      cdata,  #  list of the form [x,n], where x is an irredicible character and
              #  n is the number of times is appears in data
      x,n,i,  #  counters
      wd,     #
      L;      #

G := UnderlyingMagma(FG);
F := LeftActingDomain(FG);
p := Characteristic(F);
m := Log(Size(F),p);
irr := Irr(G);

lex:=function(x,y) 
 return x[1]<y[1] or (x[1]=y[1] and x[2]<y[2]); 
 end;
 
data := List(irr,x->[x[1],Lcm(m,Log(SizeOfSplittingField(x,p),p))]);
Sort(data,lex);
cdata := [];
while data <> [] do
    x:=data[1];
    n:=0;
    while data<>[] and x=data[1] do
        n:=n+1;
        Remove(data,1);
    od;
    Add(cdata,[x,n]);
od;

wd := [];
for x in cdata do;
    for i in [1..m*x[2]/x[1][2]] do
        if IsCheapConwayPolynomial(p,x[1][2]) then
            L := GF( p, ConwayPolynomial(p,x[1][2]) );
        else
            L := GF( p, RandomPrimitivePolynomial(p,x[1][2]) ); 
        fi;
        Add(wd, FullMatrixAlgebra(L, x[1][1]));
    od;
od;

Info( InfoWedderga, 2, "Info version : ", WedderburnDecompositionInfo( FG ) );

return wd;

end);    


#############################################################################
##
#A WedderburnDecompositionInfo( FG ) 
##
## The function WedderburnDecompositionInfo compute a list of numerical data 
## describing the Wedderburn components of the semisimple group algebra FG over 
## a cyclotomic field, and stores the result as an attribute of FG. 
##
InstallMethod( WedderburnDecompositionInfo , 
    "for semisimple group algebra over cyclotomic fields", 
    true, 
    [ IsSemisimpleANFGroupAlgebra ], 
    0,
function( FG )

local   G,      # Group
        F,      # Coefficient field
        pairs,  # Strong Shoda pairs of G
        A,      # Simple algebra
        i,      # Counter
        exp,    # Exponent of G
        br,     # List of lists of strongly Shoda triples
        sst,    # an element of sst
        chi,    # an irreducible character
        cf,     # character field of chi
        output,
        x;
        
G := UnderlyingMagma(FG);
F:=LeftActingDomain(FG);
output := [];

if IsSemisimpleANFGroupAlgebra(FG) then
    
    for i in GenWeddDecomp(FG) do
      A := SimpleAlgebraInfoByData(i);
      Append( output, [ A ] );
    od;  

    if ForAny( output, x -> not IsInt(x[1]) ) then
        Print("Wedderga: Warning!!! \n", 
        "Some of the Wedderburn components displayed are FRACTIONAL MATRIX ALGEBRAS!!!\n\n");
    fi;
    
    return output;
    
else

    Error("Wedderga: <FG> must be a group algebra over a cyclotomic field!!!");

fi;

end); 

#############################################################################
##
#A  WedderburnDecompositionInfo( FG ) 
##
InstallMethod( WedderburnDecompositionInfo , 
    "for semisimple finite group algebra",
    true, 
    [ IsSemisimpleFiniteGroupAlgebra ], 
    0,
function( FG )

local G,      # Underlying group of FG
      F,      # Coefficient field of FG
      p,      # Characteristic of the field F
      m,      # Power of p in the size of the field F
      irr,    # Irreducible characters of G
      lex,    # lexicographical ordering function
      data,   #  list of 2-tuples, the first is the degree of the character x and 
              #  the second is the lcm(m, the power of p in the field where x can 
              #  be realized
      cdata,  #  list of the form [x,n], where x is an irredicible character and
              #  n is the number of times is appears in data
      x,n,i,  #  counters
      output;     

G := UnderlyingMagma(FG);
F := LeftActingDomain(FG);
p := Characteristic(F);
m := Log(Size(F),p);
irr := Irr(G);

lex:=function(x,y) 
 return x[1]<y[1] or (x[1]=y[1] and x[2]<y[2]); 
 end;
 
data := List(irr,x->[x[1],Lcm(m,Log(SizeOfSplittingField(x,p),p))]);
Sort(data,lex);
cdata := [];
while data <> [] do
    x:=data[1];
    n:=0;
    while data<>[] and x=data[1] do
        n:=n+1;
        Remove(data,1);
    od;
    Add(cdata,[x,n]);
od;

output := [];
for x in cdata do;
    for i in [1..m*x[2]/x[1][2]] do
        Add(output,[x[1][1],p^x[1][2]]);
    od;
od;

if ForAny( output, x -> not IsInt(x[1]) ) then
    Print("Wedderga: Warning!!! \n", 
    "Some of the Wedderburn components displayed are FRACTIONAL MATRIX ALGEBRAS!!!\n\n");
fi;
    
return output;

end);


#############################################################################
##
#O GenWeddDecomp( KG )
## 
## The function returns information about the Wedderburn decomposition of
## zero-characteristic group algebra KG in the form of list of 2-tuples
## or 5-tuples, where each tuple contains the following information:
## (in the case of a 2-tuple we consider only the first two entries of 
## a 5-tuple):
## 1st position = the size of the matrices
## 2nd position = the centre of the simple component
## 3rd position = integer that is the order of the root of unity
## 4th position = Galois group of a crossed product
## 5th position = the cocycle
## The function uses WeddDecompData(G), and in the case of K=Rationals
## this is the output. 
##
InstallMethod( GenWeddDecomp,
"for semisimple infinite group algebras",
true,
[ IsSemisimpleANFGroupAlgebra ],
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
                                    List(PreImagesNC(redmok,y),w->Int(w^redmcondK))
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
## The attribute stores data for a group G using the function 
## AddCrossedProductBySSP and, in a non strongly monomial case, also using
## the function AddCrossedProductBySST. The input of AddCrossedProductBySST
## uses the function BWNoStMon. The output is a list, each entry of which 
## is either 2-tuples or 5-tuple.
## The 2-tuple contains the following data:
## 1st position = the size of the matrices
## 2nd position = the cyclotomic field = the center of the simple component
## The 5 tuple contains the following data:
## 1st position = the size of the matrices
## 2nd position = the cyclotomic field = the center of the simple component
## 3rd position = an integer that is the index (K:H) in a strongly monomial case
##                or the conductor in the other case, and it gives us the order
##                of the root of unity to be used
## 4th position = the Galois group of the cyclotomic extension
## 5th position = the cocycle
##
InstallMethod( WeddDecompData,
"for numerical data for decomposition of semisimple infinite group algebras",
true,
[ IsGroup ],
0,
function(G)

local output, # the output
         exp, # the exponent of G
          br, # the information given by BWNoStMon(G)
         sst, # current element from br
         chi, # character that is the 1st entry of sst
          cf; # cyclotomic field that is the 2nd entry of sst

if IsAbelian(G) then 
      return List( RationalClasses(G), x -> [ 1, CF(Order(Representative(x))) ] );
fi;

output :=  List( StrongShodaPairs(G), x -> 
             AddCrossedProductBySSP(G,x[1],x[2]));

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
## The arguments are:
##  exp = the exponent of the group
##    n = an integer ...
##   cf = cyclotomic field
##  Gal = the Galois group of the cyclotomic extension
## LSST = a list of strongly Shoda triples needed to
##        describe the simple component
## Returns a list, each entry of which is either a 2-tuple or a 5-tuple
## containing the following data:
## for the 2-tuple:
## 1st position = the size of the matrices
## 2nd position = the cyclotomic field = the center of the simple component
## for the 5-tuple:
## 1st position = the size of the matrices
## 2nd position = the cyclotomic field = the center of the simple component
## 3rd position = an integer that is the index (K:H) in a strongly monomial case
##                or the conductor in the other case, and it gives us the order
##                of the root of unity to be used
## 4th position = the Galois group of the cyclotomic extension
## 5th position = the cocycle
##
InstallMethod(AddCrossedProductBySST,
"for semisimple infinite group algebras",
true,
[ IsInt, IsInt, IsField, IsGroup, IsList ],
0,
function( exp, n, cf , Gal , LSST )
local
	N,Epi,NH,KH,k,ok,
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
		ok := Index(x[2],x[3]);
		N := Normalizer(x[1],x[3]);
		Epi := NaturalHomomorphismByNormalSubgroup( N, x[3] ) ;
    NH  := Image(Epi,N);
    KH  := Image(Epi,x[2]);
	if Size(KH)=1 then
  		k:=One(KH);
	else
  		if IsCyclic(KH) then
    		k := MinimalGeneratingSet(KH)[1];
  		else
    		Error("One of the entries of the fifth input is not a strong Shoda triple!");
  		fi;
	fi;
	
    Add( LC, CocycleByData(exp,Galnum,cf,x[1],x[2],x[3],a,N,Epi,NH,KH,ok,k) );
  od;
  
  Cond := Lcm(List(LC,x->x[1]));
  GalCond := Subgroup(Units(ZmodnZ(Cond)),Image(ReductionModnZ(exp,Cond),Galnum));
  coc := function(a,b)
    local out, x, redu;  
      out := Zero( ZmodnZ( Cond ) );
      for x in LC do
        redu := ReductionModnZ( Cond, x[1] );
        out := out + (Cond/x[1]) * ZmodnZObj(Int(x[2](a^redu,b^redu)),Cond);
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
#O SimpleAlgebraByStrongSP( QG, K, H )
##
## The function SimpleAlgebraByStrongSP computes the simple algebras 
## QG*e( G, K, H) if ( K, H ) is a SSP of G 
## This version does not check the input
##
InstallOtherMethod( SimpleAlgebraByStrongSP, 
    "for semisimple rational group algebras", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function( QG, K, H)
if IsStrongShodaPair( UnderlyingMagma( QG ), K, H ) then
    return SimpleAlgebraByStrongSPNC( QG, K, H );
else
    Error("Wedderga: <(K,H)> should be a strongly Shoda pair of the underlying group of <QG>\n");
fi;
end);


#############################################################################
##
#O SimpleAlgebraByStrongSP( FqG, K, H, C ) 
##
## The function SimpleAlgebraByStrongSP verifies if ( H, K ) is a SSP of G and
## C is a cyclotomic class of q=|Fq| module n=[K:H] containing generators
## of K/H, and in that case computes the simple algebra  FqG*e( G, K, H, C)
##
InstallMethod( SimpleAlgebraByStrongSP, 
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

if not(IsStrongShodaPair(G, K, H )) then
    Error("Wedderga: (<K>,<H>) should be a strongly Shoda pair of the underlying group of <FqG>\n");
elif IsCyclotomicClass( q, n, C) and Gcd(n,C[1]) =1 then
    return SimpleAlgebraByStrongSPNC( FqG, K, H, C );
else Error("Wedderga: <C> should be a generating cyclotomic class module the index of <H> in <K>\n");
fi;

end);

#############################################################################
##
#O SimpleAlgebraByStrongSP( FqG, K, H, c ) 
##
## The function SimpleAlgebraByStrongSP verifies if ( H, K ) is a SSP of G and
## c is an integer coprime with n=[K:H]. 
## If the answer is positive then returns SimpleAlgebraByStrongSP(FqG, K, H, C)
## where C is the cyclotomic class of q=|Fq| module n=[K:H] containing c.
##
InstallOtherMethod( SimpleAlgebraByStrongSP, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsPosInt ], 
    0,
function( FqG, K, H, c )
local   G,      # Group
        n;      # Index of H in K        
        
G := UnderlyingMagma( FqG );
n := Index( K, H );

if  IsStrongShodaPair(G, K, H ) then
  if Gcd( c, n ) = 1 then
    return SimpleAlgebraByStrongSPNC( FqG, K, H, c  mod n);
  else
    Error("Wedderga: <c> should be coprime with the index of <H> in <K>");   
  fi;
else
   Error("Wedderga: (<K>,<H>) should be a strongly Shoda pair of the underlying group of <FqG>\n");
fi;
end);


#############################################################################
##
#O SimpleAlgebraByStrongSPNC( QG, K, H )
##
## The function SimpleAlgebraByStrongSPNC computes simple algebras 
## QG*e( G, K, H), for ( K, H ) a SSP of G 
## This version does not check the input
##
InstallOtherMethod( SimpleAlgebraByStrongSPNC, 
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
        R,          # Crossed product       
        act,        # Action for the crossed product
        coc;        # Twisting for the crossed product
        
G := UnderlyingMagma( QG );
N   := Normalizer(G,H);
ind := Index(G,N);
ok := Index( K, H );
if N=K then
    if ind=1 then # G=N
        Info( InfoWedderga, 2, "N_G(H) = K = G, returning CF(", ok, ")");
        return CF(ok);
    else
        Info( InfoWedderga, 2, "N_G(H) = K <> G, returning M_", 
              ind, "( CF(", ok, ") )");
        return FullMatrixAlgebra( CF(ok), ind );
    fi;                          
else # if N_G(H) <> K
    Epi := NaturalHomomorphismByNormalSubgroup( N, H ) ;
    NH  := Image(Epi,N);
    KH  := Image(Epi,K);
	if Size(KH)=1 then
  		k:=One(KH);
	else
  		if IsCyclic(KH) then
    		k := MinimalGeneratingSet(KH)[1];
  		else
    		Error("Second input modulo the third one must be a cyclic group!");
  		fi;
	fi;

	
    Potk:= [ k ];
    for i in [ 2 .. ok ] do
        Potk[i] := Potk[i-1]*k; 
    od;
    Epi2:=NaturalHomomorphismByNormalSubgroup( NH, KH ) ;
    NdK:=Image(Epi2,NH);
        
      act := function( RG, a ) 
             local x, ok, Potk, Epi2;
             ok   := OperationRecord(RG).ok;
             Potk := OperationRecord(RG).Potk;
             Epi2 := OperationRecord(RG).Epi2;
             return MappingByFunction( CF(ok), CF(ok), x -> 
               GaloisCyc(x, Position(Potk,k^PreImagesRepresentativeNC(Epi2,a))));
             end;
               
      coc := function( RG, a, b )
             local ok, Potk, Epi2;
             ok   := OperationRecord(RG).ok;
             Potk := OperationRecord(RG).Potk;
             Epi2 := OperationRecord(RG).Epi2;     
             return E(ok)^Position( Potk,
                                    PreImagesRepresentativeNC( Epi2, a*b )^-1 *
                                    PreImagesRepresentativeNC( Epi2, a ) *
                                    PreImagesRepresentativeNC( Epi2, b ) );
             end;   
      
    R := CrossedProduct(CF(ok), NdK, act, coc);       
                 
    SetOperationRecord( R, rec(ok:=ok, Potk:=Potk, Epi2:=Epi2) );                 
                  
    if ind=1 then
      Info( InfoWedderga, 2, "N_G(H) <> K, returning crossed product");
      return R;
    else
      Info( InfoWedderga, 2, 
        "N_G(H) <> K, returning matrix algebra over crossed product");
      return FullMatrixAlgebra( R, ind );
    fi;  
fi;      
end);


#############################################################################
##
#O AddCrossedProductBySSP( G, K, H )
## 
## Let G be a group and K,H be a strongly Shoda pair in G. The function
## returns the 2-tuple of the 5-tuple that will describe the structure 
## of the crossed product given by this SSP:
## for the 2-tuple ( if K=N, where N=N_G(H) ):
## 1st position = the size of the matrices = index (G:N)
## 2nd position = the cyclotomic field = the center of the simple component
## for the 5-tuple (if K<>N):
## 1st position = the size of the matrices
## 2nd position = the cyclotomic field = the center of the simple component
## 3rd position = an integer that is the index (K:H) in a strongly monomial case
##                or the conductor in the other case, and it gives us the order
##                of the root of unity to be used
## 4th position = the Galois group of the cyclotomic extension
## 5th position = the cocycle
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
	if Size(KH)=1 then
  		k:=One(KH);
	else
  		if IsCyclic(KH) then
    		k := MinimalGeneratingSet(KH)[1];
  		else
    		Error("Second input modulo the third one must be a cyclic group!");
  		fi;
	fi;
	
    Epi2:=NaturalHomomorphismByNormalSubgroup( NH, KH ) ;
    NdK:=Image(Epi2,NH);
    bij := MappingByFunction(ZmodnZ(ok),KH,i->k^Int(i));
    
    # The cocycle in Z^2(NdK,<E(ok)>)
    coc := function(a,b)
       return PreImagesRepresentativeNC(bij,
                              PreImagesRepresentativeNC(Epi2,a*b)^-1 *
                              PreImagesRepresentativeNC(Epi2,a) *
                              PreImagesRepresentativeNC(Epi2,b) );
       end;   

    # The cocycle in Z^2(GalSSP,<E(ok)>)
    Uok:=Units(ZmodnZ(ok));
    bijunit := MappingByFunction(Uok,KH,i->k^Int(i));

    funNdK := MappingByFunction(NdK,Uok,
        function(n) 
            return PreImagesRepresentativeNC(bijunit,
                                 k^PreImagesRepresentativeNC( Epi2 , n ) );
              end
              );
    GalSSP := Subgroup(Uok,Image(funNdK));
    cocSSP := function(a,b)
                return 
        coc(PreImagesRepresentativeNC(funNdK,a),PreImagesRepresentativeNC(funNdK,b));
                end;
    
    chi := LinCharByKernel(K,H)^G;
    cf := Field( chi );
                    
    return [ ind, cf, ok , GalSSP , cocSSP ];
fi;      
end);


#############################################################################
##
#O SimpleAlgebraByData( algdata )
##
## An argument is either a 2-tuple or a 5-tuple, with the following 
## components:
## 1st position = the size of the matrices
## 2nd position = the centre of the simple component
## 3rd position = integer that is the order of the root of unity
## 4th position = Galois group of a crossed product
## 5th position = the cocycle
##
## The output is a crossed product or the matrix algebra over the crossed 
## product, constructed using this input
##
InstallMethod( SimpleAlgebraByData,
"for semisimple infinite group algebras",
true,
[ IsList ],
0,
function( algdata )

local 
L,     # The field obtained by extension of the centre of the simple
       # component with the root of unity of degree algdata[3]
cond,  # Lcm( Conductor(L), algdata[3] );
redu,  # The reduction from cond to algdata[3]
act,   # The action
coc,   # The cocycle
R;     # The crossed product

if Length(algdata) = 2 or Size(algdata[4])=1 then
    if algdata[1] = 1 then 
        return algdata[2];
    else
        return FullMatrixAlgebra( algdata[2], algdata[1] );
    fi;
else
    L := Field(algdata[2],[E(algdata[3])]);
    cond := Lcm( Conductor(L),algdata[3] );
    redu := ReductionModnZ(cond,algdata[3]);
    
    act := function( RG, a ) 
             local cond, redu;
             cond := OperationRecord(RG).cond;
             redu := OperationRecord(RG).redu;
             return ANFAutomorphism(CF(cond),Int(PreImagesRepresentativeNC(redu,a)));
             end;
             
    coc := function( RG, a, b )
             local orderroot, cocycle;
             orderroot := OperationRecord(RG).orderroot;
             cocycle   := OperationRecord(RG).cocycle;             
              return E(orderroot)^Int(cocycle(a,b));
             end;

    R := CrossedProduct( L, algdata[4], act, coc );
    SetCenterOfCrossedProduct( R, algdata[2] ); 
    SetOperationRecord( R, rec( cond := cond, 
                                redu := redu, 
                           orderroot := algdata[3],
                             cocycle := algdata[5] ) );
            
    if algdata[1] = 1 then 
        return R;
    else
        if IsInt(algdata[1]) then
            return FullMatrixAlgebra( R, algdata[1] );
        else 
            # Print("wedderga: Warning!\nThe output is a FRACTIONAL MATRIX ALGEBRAS!!!!\n");
            return [ algdata[1], R ];
        fi;
    fi;
fi;

end);


#############################################################################
## 
#O SimpleAlgebraByCharacter( FG, chi ) 
#
# The input is a semisimple infinite group algebra and an irreducible character
# of the finite group G.
#
# The output is a crossed product or the matrix algebra over the crossed 
## product, the simple component of FG gven by the character chi.
##
InstallMethod( SimpleAlgebraByCharacter,
"for semisimple infinite group algebras",
true,
[ IsSemisimpleANFGroupAlgebra, IsCharacter ],
0,
function( FG, chi )
 local G,               # underlying group 
       ratchi,          # rationalized of chi
       L,          	    # Splitting Field of G
       sspsub,     	    # List of pairs [p,SST] where p is a set of primes and 
                        # SST is a strongly Shoda triple such that the simple 
                        # algebra associated to SST is the p-part of one 
                        # Wedderburn component of QG
        sylow,      	  # the list of Sylow subgroups of Gal
        i,          	  # counter
        cf,         	  # character field of chi, 
        Gal,        	  # Gal(L/cf)
        d,          	  # integer
        pr,         	  # prime divisors of d
        sub,        	  # Conjugacy Classes of subgroups of G
        nsub,       	  # Cardinality of sub
        subcounter, 	  # counter for sub
        M,          	  # subgroup of G
        ssp,        	  # strongly Shoda pairs of M
        m,          	  # Size of ssp
        sspcounter, 	  # counter for ssp
        K,H,        	  # strongly Shoda pair of M
        psi,        	  # the strongly monomial character of M given by M
        cfpsi,      	  # character field of psi
        gencfpsi,   	  # generators of character field of psi
        dropprimes,	    # list of primes to be drop from primes
        remainingprimes,# counter of remaining primes
        primecounter,   # primes counter
        p,          	  # element of primes[controlcounter]
        P,          	  # p-Sylow subgroup of GalList[controlcounter]
        genP,       	  # set of generators of P
        x,                # 5-tuples, output of AddCrossedProductBySST
        sprod,      	  # (chi_M,psi)
        F1,n,alg,ok,Gal1,coc,F2,F3,a1,b1,Fxi,d1,condK,redmok,redmcondK,gal,x1; # For bugfix


# if not IsSemisimpleZeroCharacteristicGroupAlgebra( FG ) then
#   Error("<FG> must be a zero-characteristic semisimple group algebra !!!");       
# fi;   
    
  G := UnderlyingMagma(FG);      
  cf := Field( chi );                          
  L := CF(Exponent(G));   
  sspsub:=[];
  Gal := GaloisGroup(AsField(cf,L));
  d:=Gcd(Size(Gal),chi[1]);
  
  if  d = 1 then 
      sspsub:=[chi,cf];
      pr:=[];
  else 
      pr := Set(FactorsInt(d));
      sspsub:=[chi,cf,[],Gal];
      sylow:=List(pr,p->SylowSubgroup(Gal,p));
  fi;    
  
  sub:=ConjugacyClassesSubgroups(G);
  if ForAny( [1 .. Length(sub)-1 ], i -> 
             Size(Representative(sub[i])) > Size(Representative(sub[i+1]))) 
             then
               sub:=ShallowCopy(ConjugacyClassesSubgroups(G));
               Sort(sub, function(v,w) return Size(Representative(v))<
               Size(Representative(w)); 
              end);
  fi;  
  
  nsub := Size(sub);
  subcounter := nsub;
  while Length(pr) > 0 do
        M:=Representative( sub[ subcounter ] );
        ssp := StrongShodaPairs(M);
        m := Length(ssp);
        sspcounter := 1;
        while sspcounter <= m and Length(pr) > 0 do
            K := ssp[sspcounter][1]; 
            H := ssp[sspcounter][2];
            psi := LinCharByKernel(K,H)^M;
            cfpsi := Field(psi);
            gencfpsi := GeneratorsOfField(cfpsi);
              dropprimes := [];
              remainingprimes := Length(pr);
            primecounter := 1;
               while primecounter <= remainingprimes do
                    p := pr[primecounter];
                    P := sylow[primecounter];
                    genP := GeneratorsOfGroup(P);
                    if ForAll(Cartesian(genP,gencfpsi), x -> x[2]^x[1]=x[2]) 
                        then 
                          sprod := ScalarProduct( Restricted(chi,M),
                                      ClassFunction(M,RationalizedMat([psi])[1]));
                          if sprod mod p <> 0 then
                             Add(dropprimes,p);
                          fi;
                    fi;
                    primecounter := primecounter+1;
                od;
                pr:= Difference(pr,dropprimes);
                if dropprimes <> [] then
                Add(sspsub[3],[M,K,H,dropprimes]);
                fi;

        sspcounter := sspcounter + 1;
        od;
        subcounter:=subcounter-1;
  od;
# bugfix for Field adjustment
    
  if Length(sspsub)=2 then
        F1:=LeftActingDomain(FG);
	a1:=PrimitiveElement(F1);
	b1:=PrimitiveElement(cf);
	F2:=Field([a1,b1]);
    return  SimpleAlgebraByData( [ sspsub[1][1], F2 ] );
    
  elif Size(sspsub[4])=1 then
    return [ sspsub[1][1], F2 ]; # sspsub[3] ];   
  else
    x:=AddCrossedProductBySST( Exponent(G), 
                              sspsub[1][1], 
                              sspsub[2], 
                              sspsub[4], 
                              sspsub[3]); 
    alg:=SimpleAlgebraByData(x);
# Field adjustment
    if IsField(alg) or ( IsMatrixFLMLOR(alg) and IsField(LeftActingDomain(alg)) ) then
        F1:=LeftActingDomain(FG);
	    a1:=PrimitiveElement(F1);
        b1:=PrimitiveElement(cf);
        F2 :=  Field([a1,b1]);
        
        return FullMatrixAlgebra(F2,sspsub[1][1]);
    else 
	F1:=LeftActingDomain(FG);
	a1:=PrimitiveElement(F1);
        b1:=PrimitiveElement(cf);
     if (a1 in cf) then 
        return SimpleAlgebraByData(x);
     else
	n := x[1];
      	#z := x[2]; z is cf
	ok := x[3];
        Gal1 := x[4];
        coc := x[5];
        b1:=PrimitiveElement(cf);
        F2 :=  Field([a1,b1]);
	F3 :=  Field([a1,b1,E(ok)]);  
        #a := Dimension(z)*Dimension(K)/Dimension(F); #a not needed, only one component
        Fxi := Field([b1,E(ok)]);
        d1 := (Dimension(Fxi)*Dimension(F2))/(Dimension(cf)*Dimension(F3));
	#  Dimension(Field([b1,E(ok)]))/Dimension(Fxi);
        condK := Conductor(F1);
        m := Lcm(condK,ok);
        redmok := ReductionModnZ(m,ok);
        redmcondK := ReductionModnZ(m,condK);
        gal := Subgroup(Units(ZmodnZ(ok)),
                    Filtered(Gal1,y->
                            Size(
                                Intersection(
                                    GaloisStabilizer(F2),
                                    List(PreImagesNC(redmok,y),w->Int(w^redmcondK))
                                            )
                                )<>0
                            )
                        );
        #for i in [1..a] do
            x1:=[n*d1,F2,ok,gal,coc];
        #od;
      #return SimpleAlgebraByData( x1 );  

##############################
  
# if Length(sspsub)=2 then
#    return SimpleAlgebraByData( [ sspsub[1][1], sspsub[2] ] );
# else
#    x:=AddCrossedProductBySST( Exponent(G), 
#                               sspsub[1][1], 
#                               sspsub[2], 
#                               sspsub[4], 
#                               sspsub[3]);
     if not IsInt(x1[1]) then 
      Print("Wedderga: Warning!\nThe output is a FRACTIONAL MATRIX ALGEBRA!!!\n\n");
     fi;                         
     return SimpleAlgebraByData(x1);  
 fi;
fi;
fi;
 
end);


#############################################################################
## 
#O SimpleAlgebraByCharacter( FG, chi ) 
##
## The input is a semisimple infinite group algebra and an irreducible character
## of the finite group G.
##
## The output is a crossed product or the matrix algebra over the crossed 
## product, the simple component of FG gven by the character chi.
##
InstallMethod( SimpleAlgebraByCharacter,
"for semisimple finite group algebras",
true,
[ IsSemisimpleFiniteGroupAlgebra, IsCharacter ],
0,
function( FG, chi )

local G,      # Underlying group of FG
      F,      # Coefficient field of FG
      p,      # Characteristic of the field F
      m,      # Power of p in the size of the field F
      power,  #  lcm(m, the power of p in the field where chi can be realized)
      alg;     #

G := UnderlyingMagma(FG);
F := LeftActingDomain(FG);
p := Characteristic(F);
m := Log(Size(F),p);

power := Lcm(m,Log(SizeOfSplittingField(chi,p),p));
alg := FullMatrixAlgebra(GF(p^power), chi[1]);

return alg;

end);


#############################################################################
##
#O SimpleAlgebraByStrongSPNC( FqG, K, H, C )
##
## The function SimpleAlgebraByStrongSPNC computes simple algebras 
## FqG*e( G, K, H, C), for ( H, K ) a SSP of G and C a cyclotomic class 
## of q=|Fq| module n=[K:H] containing generators of K/H.
## This version does not check the input
##
InstallMethod( SimpleAlgebraByStrongSPNC, 
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
#O SimpleAlgebraByStrongSPNC( FqG, K, H, c ) 
##
## The function SimpleAlgebraByStrongSP verifies if ( H, K ) is a SSP of G and
## c is an integer coprime with n=[K:H]. 
## In the answer is positive then return SimpleAlgebraByStrongSP(FqG, K, H, C) 
## where C is the cyclotomic class of q=|Fq| module n=[K:H] containing c.
##
InstallOtherMethod( SimpleAlgebraByStrongSPNC, 
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
    return SimpleAlgebraByStrongSPNC( FqG, K, H, C );

end);


#############################################################################
##
#O SimpleAlgebraByStrongSPInfo( QG, K, H ) 
##
## The function SimpleAlgebraByStrongSPInfo compute the data describing simple
## algebras QG*e( G, K, H ), for ( H, K ) a SSP of G, but first verify the input
##
InstallOtherMethod( SimpleAlgebraByStrongSPInfo, 
    "for semisimple rational group algebras", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ], 
    0,
function( QG, K, H )
if  IsStrongShodaPair( UnderlyingMagma( QG ), K, H ) then
    return SimpleAlgebraByStrongSPInfoNC( QG, K, H );
else
    Error("Wedderga: <(K,H)> should be a strongly Shoda pair of the underlying group of <QG>\n");
fi;
end);


#############################################################################
##
#O SimpleAlgebraInfoByData( x )
##
## An argument is either a 2-tuple or a 5-tuple, with the following 
## components:
## 1st position = the size of the matrices
## 2nd position = the centre of the simple component
## 3rd position = integer that is the order of the root of unity
## 4th position = Galois group of a crossed product
## 5th position = the cocycle
##
## The output is list of 2, 4 or 5 elements:
## 1st position = the size of the matrices
## 2nd position = the centre of the simple component
## 3rd position = integer that is the order of the root of unity
## 4th position = a list of 3 elements:
##                1st position 
##                2nd position
##                3rd position
##
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
c,              # Value of cocycle 
A;              # Info of the algebra before using GlobalSplittingOfCyclotomicAlgebra and SchurIndex

if Length(x) = 2 then 
    return x;
# elif Size(x[4])=1 then
#    return [ x[1], x[2] ];
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
        A:= [ x[1],                          # the size of matrices
                 x[2],                          # the centre of the simple component
                 Cond,                          # the order of the root of unity
                 [ o[1], Int(Gen[1]) , beta[1]] #
               ]; 
                   
    else
        A:= [ x[1],
                 x[2],
                 Cond,
                 List([1..Length(Gen)], i -> [ o[i], Int(Gen[i]) , beta[i] ] ),
                 List( [1..Length(Gen)-1], i -> 
                     List( [i+1..Length(Gen)], j -> 
                         Int(coc(Gen[j],Gen[i])-coc(Gen[i],Gen[j]))
                          )
                     )
                ];        
    fi;       
  return GlobalSplittingOfCyclotomicAlgebra(A);
  if Length(A) = 2 or SchurIndex(A)<>1 then 
    return A;
  else
    if Length(A)=4 then 
      A[1]:=A[1]*A[4][1];
    else 
      A[1]:=A[1]*Product(List(A[4],x->x[1]));
    fi;
    return [A[1],A[2]];
  fi;
fi;

end);


#############################################################################
##
#O SimpleAlgebraByCharacterInfo( FG, chi )
##
# The input is an infinite group algebra FG and chi an irreducible character of a 
# finite group G.
# 
# The output is a list of 2, 3, 4 or 5 elements that describe the simple 
# algebra given by the character chi, in the following form:
## 1st position = the size of the matrices
## 2nd position = the centre of the simple component
## 3rd position = integer that is the order of the root of unity
## 4th position = a list of 3 elements:
##                1st position 
##                2nd position
##                3rd position
#
InstallMethod( SimpleAlgebraByCharacterInfo,
"for semisimple infinite group algebras",
true,
[ IsSemisimpleANFGroupAlgebra, IsCharacter ],
0,
function( FG, chi )


local G,               # underlying group 
       ratchi,          # rationalized of chi
       L,          	    # Splitting Field of G
       sspsub,     	    # List of pairs [p,SST] where p is a set of primes and 
                        # SST is a strongly Shoda triple such that the simple 
                        # algebra associated to SST is the p-part of one 
                        # Wedderburn component of QG
        sylow,      	  # the list of Sylow subgroups of the elements in GalList
        i,          	  # counter
        cf,         	  # character field of chi
        Gal,        	  # Gal(L/cf)
        d,          	  # integer
        pr,         	  # prime divisors of d
        sub,        	  # Conjugacy Classes of subgroups of G
        nsub,       	  # Cardinality of sub
        subcounter, 	  # counter for sub
        M,          	  # subgroup of G
        ssp,        	  # strongly Shoda pairs of M
        m,          	  # Size of ssp
        sspcounter, 	  # counter for ssp
        K,H,        	  # strongly Shoda pair of M
        psi,        	  # the strongly monomial character of M given by M
        cfpsi,      	  # character field of psi
        gencfpsi,   	  # generators of character field of psi
        dropprimes,	    # list of primes to be drop from primes
        remainingprimes,# counter of remaining primes
        primecounter,   # primes counter
        p,          	  # element of primes[controlcounter]
        P,          	  # p-Sylow subgroup of GalList[controlcounter]
        genP,       	  # set of generators of P
        sprod,      	  # (chi_M,psi)
	x,F1,alg,n,ok,Gal1,coc,F2,F3,a1,b1,Fxi,d1,condK,redmok,redmcondK,gal,x1; # For bugfix

# if not IsSemisimpleZeroCharacteristicGroupAlgebra( FG ) then
#   Error("<FG> must be a zero-characteristic semisimple group algebra !!!");       
# fi;   
    
  G := UnderlyingMagma(FG);    
  ratchi:=RationalizedMat([chi])[1];
  cf := Field( chi );
  
  L := CF(Exponent(G));  
  Gal := GaloisGroup(AsField(cf,L));
  
  d:=Gcd(Size(Gal),chi[1]);
      if  d = 1 then 
          sspsub:=[chi,cf];
          pr:=[];
      else 
          pr := Set(FactorsInt(d));
          sspsub:=[chi,cf,[],Gal];
          sylow:=List(pr,p->SylowSubgroup(Gal,p));
       fi;    
  
  sub:=ConjugacyClassesSubgroups(G);
  if ForAny( [1 .. Length(sub)-1 ], i -> 
             Size(Representative(sub[i])) > Size(Representative(sub[i+1]))) then
    sub:=ShallowCopy(ConjugacyClassesSubgroups(G));
    Sort(sub, function(v,w) return Size(Representative(v))<
    Size(Representative(w)); 
    end);
  fi;  
  
  nsub := Size(sub);
  subcounter := nsub;
  while Length(pr) > 0 do
        M:=Representative( sub[ subcounter ] );
        ssp := StrongShodaPairs(M);
        m := Length(ssp);
        sspcounter := 1;
        while sspcounter <= m and Length(pr) > 0 do
            K := ssp[sspcounter][1]; 
            H := ssp[sspcounter][2];
            psi := LinCharByKernel(K,H)^M;
            cfpsi := Field(psi);
            gencfpsi := GeneratorsOfField(cfpsi);
              dropprimes := [];
              remainingprimes := Length(pr);
            primecounter := 1;
               while primecounter <= remainingprimes do
                    p := pr[primecounter];
                    P := sylow[primecounter];
                    genP := GeneratorsOfGroup(P);
                    if ForAll(Cartesian(genP,gencfpsi), x -> x[2]^x[1]=x[2]) 
                        then 
                          sprod := ScalarProduct( Restricted(chi,M),
                                      ClassFunction(M,RationalizedMat([psi])[1]));
                          if sprod mod p <> 0 then
                             Add(dropprimes,p);
                          fi;
                    fi;
                    primecounter := primecounter+1;
                od;
                pr:= Difference(pr,dropprimes);
                if dropprimes <> [] then
                Add(sspsub[3],[M,K,H,dropprimes]);
                fi;

        sspcounter := sspcounter + 1;
        od;
        subcounter:=subcounter-1;
    od;
    
  if Length(sspsub)=2 then
        F1:=LeftActingDomain(FG);
	a1:=PrimitiveElement(F1);
	b1:=PrimitiveElement(cf);
	F2:=Field([a1,b1]);
    return  [ sspsub[1][1], F2 ] ;
    
  elif Size(sspsub[4])=1 then
    return [ sspsub[1][1], F2 ]; # sspsub[3] ];   
  else
    x:=AddCrossedProductBySST( Exponent(G), 
                              sspsub[1][1], 
                              sspsub[2], 
                              sspsub[4], 
                              sspsub[3]); 
# Field adjustment 
    alg:=SimpleAlgebraByData(x);
    if IsField(alg) or ( IsMatrixFLMLOR(alg) and IsField(LeftActingDomain(alg)) ) then
	    F1:=LeftActingDomain(FG);
	    a1:=PrimitiveElement(F1);
        b1:=PrimitiveElement(cf);
        F2:=Field([a1,b1]);
        return [sspsub[1][1],F2]; 
    else
	F1:=LeftActingDomain(FG);
	a1:=PrimitiveElement(F1);
        b1:=PrimitiveElement(cf);
     if (a1 in cf) then 
        return SimpleAlgebraInfoByData(x);
     else
	n := x[1];
      	#z := x[2]; z is cf
	ok := x[3];
        Gal1 := x[4];
        coc := x[5];
        F2 :=  Field([a1,b1]);
	F3 :=  Field([a1,b1,E(ok)]);  
        #a := Dimension(z)*Dimension(K)/Dimension(F); #a not needed, only one component
        Fxi := Field([b1,E(ok)]);
        d1 := (Dimension(Fxi)*Dimension(F2))/(Dimension(cf)*Dimension(F3));
#  Dimension(Field([b1,E(ok)]))/Dimension(Fxi);
        condK := Conductor(F1);
        m := Lcm(condK,ok);
        redmok := ReductionModnZ(m,ok);
        redmcondK := ReductionModnZ(m,condK);
        gal := Subgroup(Units(ZmodnZ(ok)),
                    Filtered(Gal1,y->
                            Size(
                                Intersection(
                                    GaloisStabilizer(F2),
                                    List(PreImagesNC(redmok,y),w->Int(w^redmcondK))
                                            )
                                )<>0
                            )
                        );
        #for i in [1..a] do
            x1:=[n*d1,F2,ok,gal,coc];
        #od;
# End Field adjustment
      return SimpleAlgebraInfoByData( x1 );    
      fi; 
  fi;  
  fi;

end);

#############################################################################
##
#O SimpleAlgebraByCharacterInfo( FG, chi )
##
# The input is a finite group algebra FG and chi an irreducible character of a 
# finite group G.
# 
# The output is a 2-tuple with the first entry the degree of the character and
# the second entry the power of p

InstallMethod( SimpleAlgebraByCharacterInfo,
"for semisimple finite group algebras",
true,
[ IsSemisimpleFiniteGroupAlgebra, IsCharacter ],
0,
function( FG, chi )

local G,      # Underlying group of FG
      F,      # Coefficient field of FG
      p,      # Characteristic of the field F
      m,      # Power of p in the size of the field F
      power,  #  lcm(m, the power of p in the field where chi can be realized)
      alg;     #
      
     
G := UnderlyingMagma(FG);
F := LeftActingDomain(FG);
p := Characteristic(F);
m := Log(Size(F),p);

power := Lcm(m,Log(SizeOfSplittingField(chi,p),p));
alg := [chi[1], p^power];

return alg;

end);


#############################################################################
##
#O SimpleAlgebraByStrongSPInfo( FqG, K, H, C )
##
## The function SimpleAlgebraByStrongSPInfo cheks that (K,H) is a strongly 
## Shoda pair of G, the underlying group of the semisimple finite group algebra
## FqG with coefficients in the field of order q and if C is a generating 
## q-cyclotomic class module n=[K:H]. In that case computes the data describing 
## the simple algebra FqG*e( G, K, H, C)
##
InstallMethod( SimpleAlgebraByStrongSPInfo, 
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

if not(IsStrongShodaPair(G, K, H )) then
    Error("Wedderga: (<K>,<H>) should be a strongly Shoda pair of the underlying group of <FqG>\n");
elif IsCyclotomicClass( q, n, C) and Gcd(n,C[1]) =1 then
    return SimpleAlgebraByStrongSPInfoNC( FqG, K, H, C );
else Error("Wedderga: <C> should be a generating cyclotomic class module the index of <H> in <K>\n");
fi;

end);


#############################################################################
##
#O SimpleAlgebraByStrongSPInfo( FqG, K, H, c )
##
## The function SimpleAlgebraByStrongSPInfo cheks that (K,H) is a strongly 
## Shoda pair of G, the underlying group of the semisimple finite group algebra
## FqG with coefficients in the field of order q and in that c is a positive
## integer coprime with n=[K:H]. In that case computes the data describing the 
## simple algebra FqG*e( G, K, H, C) for C the q-cyclotomic class module n
## containing c
##
InstallOtherMethod( SimpleAlgebraByStrongSPInfo, 
    "for semisimple finite group algebras", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsPosInt ], 
    0,
function( FqG, K, H, c )  
local   G,      # Group
        n;      # Index of H in K

G := UnderlyingMagma( FqG );

if IsStrongShodaPair(G, K, H ) then
  n := Index( K, H );
  if c<n and Gcd( c, n ) = 1 then
    return SimpleAlgebraByStrongSPInfoNC( FqG, K, H, c );
  else 
    Error("Wedderga: <c> should be coprime with the index of <H> in <K>\n");
  fi;  
else
   Error("Wedderga: (<K>,<H>) should be a strongly Shoda pair of the underlying group of <FqG>\n");
fi;

end);


#############################################################################
##
#O SimpleAlgebraByStrongSPInfoNC( QG, K, H ) 
##
## The function SimpleAlgebraByStrongSPInfoNC compute the data describing simple 
## algebras QG*e( G, K, H ), for ( H, K ) a SSP of G 
##
InstallOtherMethod( SimpleAlgebraByStrongSPInfoNC, 
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
        return [ 1, Rationals ];
    fi;
    
    # First one computes an idependent set PrimGen of generators 
    # of a Primary decomposition of N/K
    N   := Normalizer(G,H);
    if N=K then
        ok := Index( K, H );
        return [ Index(G,N), CF(ok) ];
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
        gen:=List( [ 1 .. Length(Gen) ], i -> PreImagesRepresentativeNC(Epi2,Gen[i]) );
        return [ Index(G,N), 
                 NF(ok, List( [1..Length(Gen)],i->RemInt(Position(Potk,k^gen[i]),ok))),
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
                     RemInt(Position(Potk,Comm(gen[i],gen[j])),ok))) ];
    fi;
end);


#############################################################################
##
#O SimpleAlgebraByStrongSPInfoNC( FqG, K, H, C )
##
## The function SimpleAlgebraByStrongSPInfo computes the data describing 
## the algebra FqG*e( G, K, H, C) without checking conditions on the input
##
InstallMethod( SimpleAlgebraByStrongSPInfoNC, 
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
	if Size(QKH)=1 then
  		gq:=One(QKH);
	else
  		if IsCyclic(QKH) then
    		gq := MinimalGeneratingSet(QKH)[1];
  		else
    		Error("Second input modulo the third one must be a cyclic group!");
  		fi;
	fi;

C1 := Set( List( C, ii -> gq^ii ) );
St := Stabilizer( QNH, C1, OnSets );
E := PreImage( epi, St );
ord := q^( Size( C )/Index( E, K ) );

return [ Index( G, K ), ord ];

end);


#############################################################################
##
#O SimpleAlgebraByStrongSPInfoNC( FqG, K, H, c )
##
## The function SimpleAlgebraByStrongSPInfo computes the data describing 
## the algebra FqG*e( G, K, H, C), where C is the q=|Fq|-cyclotomic class module
## [K:H] containing c, without checking conditions on the input
##
InstallOtherMethod( SimpleAlgebraByStrongSPInfoNC, 
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
    return SimpleAlgebraByStrongSPInfoNC( FqG, K, H, C );

end);


#############################################################################
##                                                                         ##
##            STRONGLY SHODA PAIRS AND IDEMPOTENTS                         ##
##                                                                         ##
#############################################################################


#############################################################################
##
#A StrongShodaPairs( G )
##
## The function StrongShodaPairs computes a list of strongly Shoda pairs 
## of the group G that covers the complete set of primitive central 
## idempotents of the rational group algebra QG realizable by strongly 
## Shoda pairs
##
InstallMethod( StrongShodaPairs, 
    "for finite group ", 
    true, 
    [ IsGroup and IsFinite ], 
    0,
function(G)
local ESSPD,SSPD,QG;

QG:=GroupRing(Rationals,G);
ESSPD:=ExtSSPAndDim(G).ExtremelyStrongShodaPairs; 
SSPD:=SSPNonESSPAndTheirIdempotents(QG).NonExtremelyStrongShodaPairs; 
return Concatenation(ESSPD,SSPD);
end);

  
#############################################################################
##
#A StrongShodaPairsAndIdempotents( FqG )
##
## The attribute StrongShodaPairsAndIdempotents of the semisimple finite 
## group algebra FqG returns a record with components StrongShodaPairs
## and PrimitiveCentralIdempotents, where 
## StrongShodaPairs = list of SSP and cyclotomic classes that covers the 
##                    set of PCIs of FqG realizable by SSPs, 
## PrimitiveCentralIdempotents = list of PCIs of FqG realizable by SSPs 
##                      and cyclotomic classes
##
InstallMethod( StrongShodaPairsAndIdempotents, 
    "for semisimple finite group algebra", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra ], 
    0,
function( FqG )                
local   G,          # Group
        Fq,         # Field (finite)
        F,          # Family of elements of FqG 
        elmsG,      # Elements of G
        q,          # Order of Fq
        zero,       # Zero of Fq
        e,          # The list of primitive central idempotents
        SSPsG,      # List of strongly Shoda pairs of G
        list,       # List SSP and cyclotomic classes
        setind,     # Set of n's 
        lltrace,    # List of ltrace's for n in setind
        lcc,        # Set of cc's 
        lorders,    # Set of o's for various n's
        lprimitives,# Set of pr's for o in lorders
        p,          # Integer
        H,K,        # Subgroups of G
        n,          # Index of H in K
        N,          # Normalizer of H in G
        epi,        # N --> N/H
        QKH,        # K/H
        gq,         # Generator of K/H
        pos,        # Positions
        cc,         # Set of cyclotomic classes of q module n
        ltrace,     # List of traces of a^c over Fq for c in representatives of cc
        o,          # The  multiplicative order of q module n
        pr,         # Primitive root of the field of order q^o
        a,          # Primitive n-th root of 1 in an extension of Fq
        i,          # Cyclotomic class of q module
        j,          # Counter
        etemp,      # List of idempotents eGKHc for different classes c and fixed K and H
        templist,   # List of some cyclotomic classes
        idemp;      # Idempotent eGKHc 

G := UnderlyingMagma( FqG  );
Fq := LeftActingDomain( FqG );
F := FamilyObj(Zero(FqG));
elmsG := Elements(G);
q := Size( Fq );
zero := Zero(Fq);
e := [ AverageSum(FqG,G) ];
SSPsG := StrongShodaPairs(G);
list := [ [ SSPsG[1][1], SSPsG[1][2], [[0]] ] ];
setind := [];
lltrace := [];
lcc := [];
lorders := [];
lprimitives := [];
for p in [ 2 .. Size(SSPsG) ] do
    H := SSPsG[p][2];
    K := SSPsG[p][1];
    n := Index(K,H);
    N := Normalizer( G, H );
    epi := NaturalHomomorphismByNormalSubgroup( N, H );
    QKH := Image( epi, K );
	IsCyclic(QKH);
    gq := MinimalGeneratingSet(QKH)[1];

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
        Add( lltrace, ltrace );
        Add( lcc, cc );
        Add( setind, n );
    fi;
    etemp := [];
    templist := [];
    for i in cc do
        if Gcd(i[1],n)=1 then
            idemp := CentralElementBySubgroups( FqG, K, H, i, ltrace, [epi, gq] );
            if not idemp in etemp then
                Add( etemp, idemp );
                Add( templist, i );
            fi;
        fi;
    od;
    Append( e, etemp );
    Add( list, [ K, H, templist ] );
od;
return rec( StrongShodaPairs := list, 
            PrimitiveCentralIdempotents := e );
end);


RedispatchOnCondition( StrongShodaPairsAndIdempotents,
  true, [ IsGroupRing ], [ IsSemisimpleFiniteGroupAlgebra ], 0 );


#############################################################################
##
#M Idempotent_eGsum( QG, K, H )
##
## The following function computes e(G,K,H)    
## Note that actually it returns a list of the form [ [K,H], eGKH ]
##
InstallMethod( Idempotent_eGsum,
    "for group algebra and two subgroups of its underlying group", 
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
        g,      # element of G
        RTNdK,  # Right transversal of G/NdK
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
        for g in RTNdK do
            if not (g in NdK) then 
                eGKH1g:=eGKH1^g; 
                if eGKH1*eGKH1g <> zero then 
                    return fail;
                else
                    eGKH:= eGKH + eGKH1g;
                fi;
            fi;   
        od;                    
    fi;
    return [ [ K, H ], eGKH ];
fi;       
end);


#############################################################################
##
#O PrimitiveCentralIdempotentsByStrongSP( FG )
##
## The attribute PrimitiveCentralIdempotentsByStrongSP computes the set of 
## primitive central idempotents of the group algebra FG, realizable by 
## strongly Shoda pairs, where FG is either a rational or finite group algebra
##
InstallMethod( PrimitiveCentralIdempotentsByStrongSP,
    "for finite group algebra",
    true, 
    [ IsSemisimpleFiniteGroupAlgebra ], 
    0,
function( FG )

local G;

G := UnderlyingMagma( FG );
if not IsStronglyMonomial(G)  then 
   Print("Wedderga: Warning!!!\nThe output is a NON-COMPLETE list of prim. central idemp.s of the input! \n");
fi;

return StrongShodaPairsAndIdempotents( FG ).PrimitiveCentralIdempotents; 
end);



#############################################################################
##                                                                         ##
##            COMPLETE SET OF OTHOGONAL PRIMITIVE IDEMPOTENTS              ##
##                                                                         ##
#############################################################################



#############################################################################
##
#M PrimitiveIdempotentsNilpotent( FG,H,K,C,args )
##
## 
## The function PrimitiveIdempotentsNilpotent computes a complete set of 
## orthogonal primitive idempotents of FGe where F is a finite field, G is 
## a finite nilpotent group such that FG is semisimple and e=e_C(G,H,K) 
## is a primitive central idempotent of FG (i.e. (H,K) is a strong 
## Shoda pair of G and and C is the |F|-cyclotomic class modulo n=[H:K] 
## (w.r.t. the generator gq of H/K)
## args = [epi,gq]
##
InstallMethod( PrimitiveIdempotentsNilpotent,
    "for pairs of subgroups, one cyclotomic class, mapping and group element", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList], 
    0,
function(FG,H,K,C,args)

local G,				# underlying group of FG 
		F,					# underlying field of FG  
		N,					# Normalizer of K in G
		epi,				# N -->N/K
		gq,					# generator of H/K
		QNK,				# N/K
		C1,					# Cyclotomic class of q module n in H/K 
		QEK,				# E/K
		St,					# Stabilizer of C in H/K
		E_,					# Set of representatives of St by epi
		T_E,				# Right transversal of E_ in G
		e,					# e( G, H, K )
		eps,				# epsilon( G, H, K ) 
 		a,					# representative of generator of H/K
		E_even,			# 2-part of E
		E_odd,			# 2'-part of E
		H_even,			# 2-part of H
		H_odd,			# 2'-part of H
		QEK_even,		# E_even/K 
		QEK_odd,		# E_odd/K
		QHK_even,		# H_even/K
		QHK_odd,		# H_odd/K
		a_even,			# representative of generator of H_even/K
		a_odd,			# representative of generator of H_odd/K
		b_odd,			# representative of generator of cyclic complement of <a_odd>/K in E_odd/K
		T_odd,			# conjugating set which yields the result	
		T_even,			# conjugating set which yields the result	
		beta,				# element in FG which yields the result	
		Compl,			# ComplementClasses of H_even/K in E_even/K
		b_even,			# representative of a generator of order 2^k of E_even/K
		c_even,			# representative of a generator of order 2 of E_even/K
		n,					# Logarithm in base 2 of order of a_even mod K
		k,					# Logarithm in base 2 of order of b_even mod K
		QEH_even,		# (E_even/K) / (H_even/K)
		epi2,				# E_even/K --> (E_even/K) / (H_even/K)
		I,					# generators of (E_even/K) / (H_even/K) if metacyclic
		M_even,			# M_even/K is complement of H_even/K in E_even/K
		counter,		# counter
		counter1,		# counter
		counter2,		# counter
		i,					# parameter in presentation of E_even/K
		s,					# parameter in presentation of E_even/K
		r,					# parameter in presentation of E_even/K
		sum,				# auxiliary variable
		S,					# list [x,y] such that x^2 + y^2 = -1, y<>0
		x,					# x in S
		y,					# y in S
		L;					# output

### test requirements

G := UnderlyingMagma(FG);
if Intersection(PrimeDivisors(C[1]),PrimeDivisors(Index(H,K))) <> []
	then Error("Wedderga: The input cyclotomic class C does not contain faithfull characters","\n");
fi;
if not IsNilpotent(G) 
	then Error("Wedderga: The input group G needs to be nilpotent","\n");
fi;
if not IsSemisimpleFiniteGroupAlgebra(FG) 
	then Error("Wedderga: The input needs to be finite semisimple group algebra","\n"); 
fi;
if not IsStrongShodaPair(G,H,K) 
	then Error("Wedderga: The input (H,K) has to be a SSP of G","\n"); 
fi;

### setup all needed objects

F := LeftActingDomain(FG);
N := Normalizer( G, K );
epi := args[1];
gq := args[2];
QNK := Image( epi, N );
C1 := Set( List( C, ii -> gq^ii ) );
St := Stabilizer( QNK, C1, OnSets );
E_ := PreImage( epi, St );
QEK := Image( epi, E_ );
T_E := RightTransversal( G, E_);

eps := IdempotentBySubgroups( FG, H, K, C, [epi,gq] );
e := Sum( List( T_E, g -> eps^g ) );

a := gq;
E_even := SylowSubgroup(E_,2);
E_odd := [];
for counter in Set(Filtered(Factors(Size(E_)),IsOddInt)) do
	E_odd:=Concatenation(E_odd,AsList(SylowSubgroup(E_,counter)));
od;
E_odd:=Subgroup(G,E_odd);
H_even := SylowSubgroup(H,2);
H_odd := [];
for counter in Set(Filtered(Factors(Size(H)),IsOddInt)) do
	H_odd := Concatenation(H_odd,AsList(SylowSubgroup(H,counter)));
od;
H_odd := Subgroup(G,H_odd);

QEK_even := Image(epi,E_even);
QEK_odd := Image(epi,E_odd);
QHK_even := Image(epi,H_even);
QHK_odd := Image(epi,H_odd);

if IsEmpty(PreImage(epi,MinimalGeneratingSet(QHK_even))) then
	a_even := Random(G)^0;
else 
	a_even := PreImage(epi,MinimalGeneratingSet(QHK_even))[1];
fi;
if IsEmpty(PreImage(epi,MinimalGeneratingSet(QHK_odd))) then
	a_odd := Random(G)^0;
else 
	a_odd := PreImage(epi,MinimalGeneratingSet(QHK_odd))[1];
fi;
if a_odd = Random(G)^0 then
	if IsEmpty(PreImage(epi,MinimalGeneratingSet(QEK_odd))) then
		b_odd := Random(G)^0;
	else
		b_odd := PreImage(epi,MinimalGeneratingSet(QEK_odd))[1];
	fi;
else 
	if Index(QEK_odd,Subgroup(QNK,[Image(epi,a_odd)]))=1 then
		b_odd := Random(G)^0;
	else
		b_odd := PreImage(epi,MinimalGeneratingSet(ComplementClassesRepresentatives(
													QEK_odd,Subgroup(QNK,[Image(epi,a_odd)]))[1]))[1];
	fi;
fi;

T_odd := [];
for counter in [0..Index(E_odd,H_odd)-1] do
	Add(T_odd,a_odd^counter);
od;

Compl := ComplementClassesRepresentatives(QEK_even,QHK_even);
if IsEmpty(Compl) then
	n := Log(Order(Image(epi,a_even)),2);
	k := Log(Size(QEK_even)/(2^n),2)-1;

 epi2 := NaturalHomomorphismByNormalSubgroup(QEK_even, 
  											Subgroup(QEK_even,[Image(epi,a_even)]) );
	QEH_even := Image(epi2,QEK_even);
	
	# we search for generators b_even and c_even of QEH_even
	if IsCyclic(QEH_even) then
		b_even := a_even^0;
		c_even := PreImagesRepresentativeNC(epi,PreImagesRepresentativeNC(epi2,MinimalGeneratingSet(QEH_even)[1]));				
		counter := 1;
		while Image(epi,c_even)^2 <> Image(epi,a_even)^(2^(n-1)*counter) do
			counter := counter+2;
		od;
		a_even := a_even^counter;
	else # we know that it is now metacyclic
		I := IndependentGeneratorsOfAbelianGroup(QEH_even);
		if Order(I[1]) = 2
			then 
				c_even := PreImagesRepresentativeNC(epi,PreImagesRepresentativeNC(epi2,I[1]));
				b_even := PreImagesRepresentativeNC(epi,PreImagesRepresentativeNC(epi2,I[2]));
		else
			c_even := PreImagesRepresentativeNC(epi,PreImagesRepresentativeNC(epi2,I[2]));
			b_even := PreImagesRepresentativeNC(epi,PreImagesRepresentativeNC(epi2,I[1]));
		fi;

		# replace b_even and c_even if needed in order to obtain the desired presentation

		s := 0;
		while Image(epi,b_even)^(2^k) <> Image(epi,a_even)^(s) do
			s := s+1;
		od;
		i := 0;
		while Image(epi,c_even)*Image(epi,b_even) <> Image(epi,a_even)^(i)*
												Image(epi,b_even)*Image(epi,c_even) do
			i := i+1;
		od;	
		r := 0;
		while Image(epi,b_even)^(-1)*Image(epi,a_even)*Image(epi,b_even) <> 
												Image(epi,a_even)^(r) do
			r := r+1;
		od;		
		
		sum := 0;
		for counter in [0 .. 2^k-1] do
			sum := sum + r^counter;
		od;
		counter1 := 0;
		while (counter1 * sum + s) mod 2^n = 0 do
			counter1 := counter1+1;
		od;		
		b_even := PreImagesRepresentativeNC(epi,Image(epi,b_even)*Image(epi,a_even)^counter1);

		counter2 := 0;
		while counter2*(r-1)+i*r mod 2^n = 0 do 
			counter2 := counter2+1;
		od;	
		c_even := PreImagesRepresentativeNC(epi,Image(epi,a_even)^counter2*Image(epi,c_even));
		
	fi;

	S := SolveEquation@(F);
	x := S[1];
	y := S[2];
	beta := AverageSum(FG,[b_even])*1/2*ElementOfMagmaRing(FamilyObj(Zero(FG)), 
												Zero(F), 
												[One(F),x,y], 
												[a_even^0,a_even^(2^(n-2)),a_even^(2^(n-2))*c_even]); 
	T_even := [];
	for i in [0..2^k-1] do
		Add(T_even,a_even^i);
		Add(T_even,c_even*a_even^i);
	od;	
else 
	if Index(QEK_even,QHK_even) = 1 then
		M_even := Group(Random(G)^0);
	else
		M_even := PreImage(epi,Compl[1]);
	fi;
	beta := AverageSum(FG,M_even);

	if IsCyclic(Image(epi,M_even)) then
		n := Log(Order(Image(epi,a_even)),2);
		k := Log(Size(Image(epi,M_even)),2);
	else
		n := Log(Order(Image(epi,a_even)),2);
		k := Log(Size(Image(epi,M_even)),2)-1;
	fi;

	T_even := [];
	if (n<=1 or IsCentral(QEK_even,Group(Image(epi,a_even^(2^(n-2)))))) and 
								IsCyclic(Image(epi,M_even))
		then
			for i in [0..2^k-1] do
				Add(T_even,a_even^i);
			od;
	else 	
		for i in [0..Index(E_even,H_even)/2-1] do
			Add(T_even,a_even^i);
			Add(T_even,a_even^(2^(n-2)+i));
		od;
	fi;
fi;

### Construct the idempotents 

L := List( Product3Lists([T_odd,T_even,T_E]) , i -> (AverageSum(FG,[b_odd])*beta*eps)^i);

return L; 
end);



#############################################################################
##
#M PrimitiveIdempotentsTrivialTwisting( FG,H,K,C,args )
##
## 
## The function PrimitiveIdempotentsTrivialTwisting computes a complete set of 
## orthogonal primitive idempotents of FGe where F is a finite field, G is a 
## finite strongly monomial group such that the twisting of the simple component 
## FGe is trivial, FG is semisimple and e=e_C(G,H,K) is a primitive central 
## idempotent of FG (i.e. (H,K) is a strong Shoda pair of G and and C is the 
## |F|-cyclotomic class modulo n=[H:K] (w.r.t. the generator gq of H/K)
## args = [epi,gq]
##
InstallMethod( PrimitiveIdempotentsTrivialTwisting,
    "for pairs of subgroups, one cyclotomic class, mapping and group element", 
    true, 
    [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList], 
    0,
function(FG,H,K,C,args)

local G,					# underlying group of FG
			F,					# underlying field of FG
			N,					# Normalizer of K in G
			epi,				# N -->N/K
			gq,					# generator of H/K
			QNK,				# N/K
			C1,					# Cyclotomic class of q module n in H/K 
			St,					# Stabilizer of C in H/K
			E_,					# Set of representatives of St by epi
			epi2,				# E --> E/H
			QEH,				# E/H
			T_2,				# Right transversal of E_ in G
			T_1,				# Right transversal of H in E_
			e,					# e( G, H, K )
			eps,				# epsilon( G, H, K ) 
			A,					# permutation matrix
			P,					# some matrix	
			xi,					# primitive [H:K]-th root of unity in F
			F_1,				# extension of F with xi
			F_2,				# subfield of F_1 of degree [E:H] 
			B,					# normal basis of F_1/F_2
			FH,					# Group ring FH
			FHeps,			# simple component FHeps
			m,					# isomorphism of fields from F_1 to FHeps, sending xi^C[1] to gq*eps
			twist,			# twisting of crossed product F_1 * E/H
			act,				# action of crossed product F_1 * E/H
			D,					# crossed product F_1 * E/H
			B1,					# F_1-basis of crossed product F_1 * E/H
			B2,					# F_2-basis of crossed product F_1 * E/H (F_2 is the center of F_1 * E/H)
			b,					# variable
			b1,					# variable
			c,					# variable
			coef,				# variable
			B3,					# F_2 basis of matrix algebra (of size [E:H]) over F_2
			x_e,				# psi^(-1)(PAP^(-1))
			i,					# counter	
			j, 					# counter	
			L;					# output list

G := UnderlyingMagma(FG);
F := LeftActingDomain(FG);

### test requirements

if Intersection(PrimeDivisors(C[1]),PrimeDivisors(Index(H,K)))<>[]
	then Error("Wedderga: The input cyclotomic class C does not contain faithfull characters","\n");
fi;
if not IsTwistingTrivial(G,H,K) 
	then Error("Wedderga: The associated twisting is not trivial","\n"); 
fi;
if not IsSemisimpleFiniteGroupAlgebra(FG) 
	then Print("Wedderga: The input needs to be finite semisimple group algebra","\n"); 
fi;
if not IsStrongShodaPair(G,H,K) 
	then Error("Wedderga: (H,K) has to be a SSP of G","\n"); 
fi;
if Size(F)^First([1..1000],m->RemInt(Size(F)^m,Index(H,K))=1)>2^16 
	then Print("Warning: this method might be not applicable since GAP has problems with handeling large finite fields");
fi;

### setup all needed objects

N := Normalizer( G, K );
epi := args[1];
gq := args[2];
QNK := Image( epi, N );
C1 := Set( List( C, ii -> gq^ii ) );
St := Stabilizer( QNK, C1, OnSets );
E_ := PreImage( epi, St );
epi2 := NaturalHomomorphismByNormalSubgroup(E_,H);
QEH := Image(epi2,E_);
T_2 := RightTransversal( G, E_);
T_1 := RightTransversal( E_, H);

eps := IdempotentBySubgroups( FG, H, K, C, [epi,gq] );
e := Sum( List( T_2, g -> eps^g ) );

# define field extensions
xi := PrimRootOfUnity(F,Index(H,K));  
F_1 := GF(F,MinimalPolynomial(F,xi));
F_2 := Filtered(Subfields(F_1),x->Size(x) = Characteristic(F_1)^(Log(Size(F_1),Characteristic(F_1))/Index(E_,H)))[1]; 
B := Basis(AsField(F_2,F_1),NormalBase(AsField(F_2,F_1)));

# define FHeps being the image of F_1
FH := Subalgebra(FG,AsList(Image(Embedding(G,FG),H)));
FHeps := Subalgebra(FH,Set(Basis(FH),b->b*eps));    
m := AlgebraHomomorphismByImages(F_1,FHeps,[xi^C[1]],[PreImagesRepresentativeNC(epi,gq)*eps]);

# define the needed matrices A and P
A := [];
for i in [1..Index(E_,H)] do
	Add(A,List([1..Index(E_,H)],i->Zero(F_2)));
od;
A[1][Index(E_,H)] := One(F_2);
for i in [1..Index(E_,H)-1] do
	A[i+1][i] := One(F_2);
od;

P := [];
for i in [1..Index(E_,H)] do
	Add(P,List([1..Index(E_,H)],i->Zero(F_2)));
od;
for i in [1..Index(E_,H)] do
	P[i][i] := -One(F_2);
	P[1][i] := One(F_2);
	P[i][1] := One(F_2);
od;

# define the crossed product F_1 * E/H
twist := function(RG,g,h) return One(LeftActingDomain(RG)); end;
act := function(RG,g) return ReturnGalElement(PreImagesRepresentativeNC(epi2,g),E_,H,K,F_1,xi); end;
D := CrossedProduct(F_1,QEH,act,twist);
B1 := Basis(D); # basis over F_1
B2 := []; # basis over F_2
for b in B do 
	for b1 in B1 do
		c:=CoefficientsAndMagmaElements(b1);
		Add(B2,ElementOfCrossedProduct(FamilyObj(Zero(D)),Zero(F_2),[b*c[2]],[c[1]]));
	od;
od;

# Map the F_2-basis of F_1 * E/H to a F_2-basis of the matrix algebra over F_2
B3 := List(B2, b-> MakeMatrixByBasis(CompositionMapping(
						LeftMultiplicationBy(CoefficientsAndMagmaElements(b)[2],F_1),
						ReturnGalElement(PreImagesRepresentativeNC(epi2,CoefficientsAndMagmaElements(b)[1]),E_,H,K,F_1,xi)),B));
B3 := Basis(MatrixAlgebra(F_2,Index(E_,H)),B3); 


### Construct the idempotents 

# The idempotents are conjugates of AverageSum(FG,T_1)*eps by T_2 and powers of x_e,
# where x_e is the preimage of P*A*P^(-1) under the mapping used in the definition of B3
# x_e can be constructed as follows: 
# 1. Compute the coefficients of P*A*P^(-1) in terms of the basis B3
# 2. The element with exactly these coefficients and corresponding basis elements of B2 is the element of F_1 * E/H mapping to P*A*P^(-1)
# 3. use the knowledge to write the element of F_1 * E/H as an element in FEeps = FHeps * E/H
x_e := Zero(FG);
c := Coefficients(B3,P*A*P^(-1));
for i in [1..Size(c)] do
	coef := CoefficientsAndMagmaElements(B2[i]); # we assume that the ordering of B_2 is the same as B_3
	x_e := x_e + MakeLinearCombination(FG, [PreImagesRepresentativeNC(epi2,coef[1])], [Image(m,coef[2])*Image(m,c[i])]); 
od;

L := Flat(List( T_2, i -> List(List([0..Index(E_,H)-1],i->x_e^i*(AverageSum(FG,T_1)*eps)*x_e^(Index(E_,H)-i)),j->j^i)));

return L; 
end);


#############################################################################
##
#E
##
