#####################################################
# Given a simple algebra output by WedderburnDecompositionInfo
# or SimpleAlgebraByCharacterInfo from wedderga,
# this program determines its actual matrix degree and division
# algebra part in terms of local indices at all primes.
#####################################################
#############################################
# Necessary arithmetic functions for subroutines
#############################################
InstallGlobalFunction( PPartOfN, function(n,p)
local i,a,b;

b:=n;
a:=0;
while b/p in PositiveIntegers do
b:=b/p;
a:=a+1;
od;

return p^a;
end);

###############################
InstallGlobalFunction( PDashPartOfN, function(n,p)
local m;

m:=n;
while m/p in PositiveIntegers do
m:=m/p;
od;

return m;
end);

#########################################
# Cyclotomic reciprocity functions for the extension
# F(E(n))/F at the prime p.  These calculate the
# splitting degree g(F(E(n))/F,p),
# residue degree f(F(E(n))/F,p), and
# ramification index e(F(E(n))/F,p).  Using
# suitable quotients of these, one can obtain
# the e, f, and g for any extension K/F of
# abelian number fields.
#########################################
InstallGlobalFunction( PSplitSubextension, function(F,n,p)
local a,y1,L,i,n1,n0,f,L1,b,F1;

a:=PrimitiveElement(F);
L:=[];
for i in [1..n] do
  if Gcd(i,n)=1 and GaloisCyc(a,i)=a then
    Add(L,i);
  fi;
od;

n1:=PDashPartOfN(n,p);
f:=1;
if n1>1 then
while not(PowerMod(p,f,n1)= 1) do
f:=f+1;
od;
fi;

n0:=PPartOfN(n,p);
L1:=[];
for b in L do
if GaloisCyc(E(n1),b)=E(n1) then
AddSet(L1,b);
else
for i in [1..f] do
if b mod n1 = PowerMod(p,i,n1) then
AddSet(L1,b);
fi;
od;
fi;
od;

F1:=NF(n,L1);
######bugfix-not returning extension of F-04/04/2020#######
y1:=PrimitiveElement(F1);
F1:=Field([a,y1]);
#########################################################

return F1;
end);

###################################
InstallGlobalFunction( SplittingDegreeAtP, function(F,n,p)
local K,g;

K:=PSplitSubextension(F,n,p);
g:=Trace(K,F,1);

return g;
end);

###################################
InstallGlobalFunction( ResidueDegreeAtP, function(F,n,p)
local K,a,n1,L,f;

K:=PSplitSubextension(F,n,p);
a:=PrimitiveElement(K);
n1:=PDashPartOfN(n,p);
L:=Field([a,E(n1)]);
f:=Trace(L,K,1);

return f;
end);

#################################
InstallGlobalFunction( RamificationIndexAtP, function(F,n,p)
local n0,n1,a,U,i,U1,e;

a:=PrimitiveElement(F);
n0:=Conductor([a,E(n)]);
U:=[];
for i in [1..n0] do
  if Gcd(i,n0)=1 and GaloisCyc(a,i)=a then
   Add(U,i);
  fi;
od;
n1:=PDashPartOfN(n0,p);
U1:=[];
for i in U do
  if GaloisCyc(E(n1),i)=E(n1) then
   Add(U1,i);
  fi;
od;
e:=Size(U1);

return e;
end);
#################################
# (27/03/2020) Character Descent Functions
# - added by Allen Herman to optimize Clifford theory reductions
# needed for the Local Index functions to give correct output,
# this fixes the bug created by the new WedderburnDecompositionInfo
# command outputting crossed product algebras that are wildly
# ramified at odd primes, which the Local Index functions were
# not designed to handle.
##################################
################################################
# 2) Add a new Global Splitting function that
#    reduces algebras whose factor set has too
#    many zeroes or is globally trivial.
################################################
InstallGlobalFunction( GlobalSplittingOfCyclotomicAlgebra, function(A)
local A1,m,F,a1,m1,a,b,c,n,m2,g,g1,g2,g3,b1,c1,a2,b2,c2,F1,f,F2,t,d,d1,b11,i,j,cont;

A1:=A;
if Length(A)=5 then

  F:=A[2];
  if ForAll(A[4],x->x[3]=0) and ForAll(A[5],x->ForAll(x,y->y=0)) then
      return [A[1]*Product(A[4],x->x[1]),F];
  fi;

  A:=A1;

###########################################################################
# (01/04/2020) PUTTING ZEROES IN THE FIFTH ENTRY OF A CYCLOTOMIC ALGEBRA
# WITH ARBITRARY NUMBER OF GENERATORS FOR THE GALOIS GROUP
###########################################################################

  if Length(A)=5 then
    A1:=KillingCocycle(A);
  fi;


  while A1<>fail do
    A:=A1;
    A1:=ReducingCyclotomicAlgebra(A);
  od;
fi;



################## OLD CODE ###################


#   if Length(A)=5 and Length(A[4])=2 then
#
#     A1:=A;
#     if not(A[5][1][1]=0) then
#       m:=A[1];
#       a1:=A[4][1][1];
#       b1:=A[4][1][2];
#       c1:=A[4][1][3];
#       a2:=A[4][2][1];
#       b2:=A[4][2][2];
#       c2:=A[4][2][3];
#       d:=A[5][1][1];
#       F:=A[2];
#       a:=PrimitiveElement(F);
#       m1:=A[3];
#       F2:=Field([a,E(m1)]);
#
#       t:=0;
#       for i in [1..m1-1] do
#         if ((((1-b2)*i+d) mod m1) = (0 mod m1)) then
#           d1:=1;
#           for j in [1..(a1-1)] do
#             d1:=1+b1*d1;
#           od;
#           A1:=[m,F,m1,[[a1,b1,((c1+i*d1) mod m1)],[a2,b2,c2]],[[0]]];
#           t:=1;
#           break;
#         fi;
#       od;
#
#       if t=0 then
#         b11:=PowerMod(b1,-1,m1);
#         for i in [1..m1-1] do
#           if ((((b11-1)*i+d) mod m1) = (0 mod m1)) then
#             d1:=1;
#             for j in [1..(a2-1)] do
#               d1:=1+b2*d1;
#             od;
#             A1:=[m,F,m1,[[a1,b1,c1],[a2,b2,((c2+i*d1) mod m1)]],[[0]]];
#             t:=1;
#             break;
#           fi;
#         od;
#       fi;
#
#     fi;

#    A:=A1;

#    if A[4][2][3]=0 and A[5][1][1]=0 then
#      m:=A[1];
#      a1:=A[4][1][1];
#      b1:=A[4][1][2];
#      c1:=A[4][1][3];
#      a2:=A[4][2][1];
#      b2:=A[4][2][2];
#      c2:=A[4][2][3];
#      F:=A[2];
#      a:=PrimitiveElement(F);
#      m1:=A[3];
#      F2:=Field([a,E(m1)]);
#
#      for m2 in [2..m1-1] do
#        if GaloisCyc(E(m1)^m2,b2)=E(m1)^m2 then
#          g1:=E(m1)^m2;
#          break;
#        fi;
#      od;

#      if OrderMod(b1,Order(g1))= a1 then
#        f:=Order(g1);
#        for m2 in [2..a1] do
#          if g1^m2 = E(m1)^c1 then
#            c1:=m2;
#            break;
#          fi;
#        od;
#        A1:=[m*a2, F, f, [a1,b1 mod f,c1]];
#      fi;

#    fi;

#    if A[4][1][3]=0 and A[5][1][1]=0 then
#      m:=A[1];
#      a1:=A[4][1][1];
#      b1:=A[4][1][2];
#      c1:=A[4][1][3];
#      a2:=A[4][2][1];
#      b2:=A[4][2][2];
#      c2:=A[4][2][3];
#      F:=A[2];
#      a:=PrimitiveElement(F);
#      m1:=A[3];
#      F2:=Field([a,E(m1)]);

#      for m2 in [2..m1-1] do
#        if GaloisCyc(E(m1)^m2,b1)=E(m1)^m2 then
#          g1:=E(m1)^m2;
#          break;
#        fi;
#      od;

#      if OrderMod(b2,Order(g1))= a2 then
#        f:=Order(g1);
#        for m2 in [2..a2] do
#          if g1^m2 = E(m1)^c2 then
#            c2:=m2;
#            break;
#          fi;
#        od;
#        A1:=[m*a1, F, f, [a2,b2 mod f,c2]];
#      fi;

#    fi;

#  fi;
#fi;

A1:=A;

if Length(A) = 4 then

  F:=A[2];
  a1:=PrimitiveElement(F);
  m1:=A[3];
  a:=A[4][1];
  b:=A[4][2];
  c:=A[4][3];

  n:=Conductor(F);
  if IsOddInt(n) then 
    n:=2*n; 
  fi;
  for m2 in [1..n] do
    if E(n)^m2 in F then
       break;
    fi;
  od;
  g:=Order((E(n)^m2)^a);

  g1:=E(m1);
  for m2 in [1..(a-1)] do
    g1:=E(m1)*g1^b;
  od;
  g1:=Order(g1);
  g2:=Order(E(m1)^c);
  g3:=Lcm(g,g1);
  if (g3/g2 in Integers) then
   A1:=[A[1]*a,F];
  fi;
fi;

return A1;
end);
#######################################################
# Finds group over which cyclotomic algebra of length 4 or 5
# is faithfully represented.
#######################################################
InstallGlobalFunction( DefiningGroupAndCharacterOfCyclotAlg, function(A)
#local l,f,a,b,c,d,g,I,g1,S,m,n,i,chi,F,u,V,U,F1;;
local l,f,a,b,c,d,g,I,g1,S,m,n,i,chi,F,u,V,U,F1,k,gen,ord,hs,rs,ss,cs,relact,relpow,relcom,rel;

l:=Length(A);
# g1:="fail";
if l=2 then
  return fail;
fi;

########## DEALING WITH ARBITRARY NUMBER OF GENERATORS FOR THE GALOIS GROUP #####################

if l=5 then
  k:=Length(A[4]);
  f:=FreeGroup(k+1);
  gen := GeneratorsOfGroup(f);
  ord:=A[3];
  hs:=List([1..k],x->A[4][x][1]);
  rs:=List([1..k],x->A[4][x][2]);
  ss:=List([1..k],x->A[4][x][3]);
  cs := A[5];
  relact := List([1..k], i -> (gen[1]^gen[i+1])*gen[1]^-rs[i]);
  relpow := List([1..k], i -> (gen[i+1]^hs[i])*gen[1]^-ss[i]);
  relcom := Concatenation(List([1..k-1], i -> List([1..k-i], j -> gen[i+j+1]^-1*gen[i+1]^-1*gen[i+j+1]*gen[i+1]*gen[1]^-cs[i][j])));
  rel := Concatenation([gen[1]^ord],relact,relpow,relcom);
  g := f/rel;
fi;

################ OLD CODE ######################

# if (l=5 and Length(A[4])=2) then
#   f:=FreeGroup("a","b","c");
#   a:=f.1;
#   b:=f.2;
#   c:=f.3;
#   g:=f/[a^A[3],b^A[4][1][1]*a^(-A[4][1][3]),c^A[4][2][1]*a^(-A[4][2][3]),b^(-1)*a*b*a^(-A[4][1][2]),
#        c^(-1)*a*c*a^(-A[4][2][2]),c^(-1)*b^(-1)*c*b*a^(-A[5][1][1])];
# fi;

# if (l=5 and Length(A[4])=3) then
#   f:=FreeGroup("a","b","c","d");
#   a:=f.1;
#   b:=f.2;
#   c:=f.3;
#   d:=f.4;
#   g:=f/[a^A[3],b^A[4][1][1]*a^(-A[4][1][3]),c^A[4][2][1]*a^(-A[4][2][3]),d^A[4][3][1]*a^(-A[4][3][3]),
#     b^(-1)*a*b*a^(-A[4][1][2]),c^(-1)*a*c*a^(-A[4][2][2]), d^(-1)*a*d*a^(-A[4][3][2]),    c^(-1)*b^(-1)*c*b*a^(-A[5][1][1]),
#     d^(-1)*b^(-1)*d*b*a^(-A[5][1][2]),d^(-1)*c^(-1)*d*c*a^(-A[5][2][1])];
# fi;

########## END OF DEALING WITH ARBITRARY NUMBER OF GENERATORS FOR THE GALOIS GROUP #####################
if (l=4) then
  f:=FreeGroup("a","b");
  a:=f.1;
  b:=f.2;
  g:=f/[a^A[3],b^A[4][1]*a^(-A[4][3]),b^(-1)*a*b*a^(-A[4][2])];
fi;

I:=IsomorphismSpecialPcGroup(g);
g1:=Image(I);

S:=[];
S[1]:=g1;

if Length(A)=2 then d:=1; fi;
if Length(A)=4 then d:=A[4][1]; F1:=NF(A[3],[A[4][2]]); fi;
if (Length(A)=5 and Length(A[4])=2) then
   d:=A[4][1][1]*A[4][2][1];
   F1:=NF(A[3],[A[4][1][2],A[4][2][2]]);
fi;
if (Length(A)=5 and Length(A[4])=3) then
   d:=A[4][1][1]*A[4][2][1]*A[4][3][1];
   F1:=NF(A[3],[A[4][1][2],A[4][2][2],A[4][3][2]]);
fi;

n:=Size(Irr(g1)) ;
m:=Trace(F1,Rationals,1);
U:=[];
for i in [1..n] do
chi:=Irr(g1)[n-i+1];
V:=ValuesOfClassFunction(chi);
F:=FieldByGenerators(V);
if V[1]/d in PositiveIntegers then
if Size(KernelOfCharacter(chi))=1 then
if FieldByGenerators(V)=F1 then
	Add(U,n-i+1);
fi;
fi;
fi;
od;
if Size(U)=m then
u:=U[1];
chi:=Irr(g1)[u];
else
chi:=U;
fi;

S[2]:=chi;

return S;
end);

#######################################################
InstallGlobalFunction( DefiningGroupOfCyclotomicAlgebra, function(A)
local l,f,a,b,c,d,g,I,g1,k,gen,ord,hs,rs,ss,cs,relact,relpow,relcom,rel;

l:=Length(A);
if l=2 then g1:=SmallGroup(1,1);
else
g1:="fail";

########## DEALING WITH ARBITRARY NUMBER OF GENERATORS FOR THE GALOIS GROUP #####################

if l=5 then
  k:=Length(A[4]);
  f:=FreeGroup(k+1);
  gen := GeneratorsOfGroup(f);
  ord:=A[3];
  hs:=List([1..k],x->A[4][x][1]);
  rs:=List([1..k],x->A[4][x][2]);
  ss:=List([1..k],x->A[4][x][3]);
  cs := A[5];
  relact := List([1..k], i -> (gen[1]^gen[i+1])*gen[1]^-rs[i]);
  relpow := List([1..k], i -> (gen[i+1]^hs[i])*gen[1]^-ss[i]);
  relcom := Concatenation(List([1..k-1], i -> List([1..k-i], j -> gen[i+j+1]^-1*gen[i+1]^-1*gen[i+j+1]*gen[i+1]*gen[1]^-cs[i][j])));
  rel := Concatenation([gen[1]^ord],relact,relpow,relcom);
  g := f/rel;
fi;

################ OLD CODE ######################

# if (l=5 and Length(A[4])=2) then
# f:=FreeGroup("a","b","c");
# a:=f.1;
# b:=f.2;
# c:=f.3;
# g:=f/[a^A[3],b^A[4][1][1]*a^(-A[4][1][3]),c^A[4][2][1]*a^(-A[4][2][3]),b^(-1)*a*b*a^(-A[4][1][2]),
#  c^(-1)*a*c*a^(-A[4][2][2]),c^(-1)*b^(-1)*c*b*a^(-A[5][1][1])];
# fi;

# if (l=5 and Length(A[4])=3) then
# f:=FreeGroup("a","b","c","d");
# a:=f.1;
# b:=f.2;
# c:=f.3;
# d:=f.4;
# g:=f/[a^A[3],b^A[4][1][1]*a^(-A[4][1][3]),c^A[4][2][1]*a^(-A[4][2][3]),d^A[4][3][1]*a^(-A[4][3][3]),
#  b^(-1)*a*b*a^(-A[4][1][2]),c^(-1)*a*c*a^(-A[4][2][2]), d^(-1)*a*d*a^(-A[4][3][2]),    c^(-1)*b^(-1)*c*b*a^(-A[5][1][1]),
# d^(-1)*b^(-1)*d*b*a^(-A[5][1][2]),d^(-1)*c^(-1)*d*c*a^(-A[5][2][1])];
# fi;

########## END OF DEALING WITH ARBITRARY NUMBER OF GENERATORS FOR THE GALOIS GROUP #####################

if (l=4) then
f:=FreeGroup("a","b");
a:=f.1;
b:=f.2;
g:=f/[a^A[3],b^A[4][1]*a^(-A[4][3]),b^(-1)*a*b*a^(-A[4][2])];
fi;

I:=IsomorphismSpecialPcGroup(g);
g1:=Image(I);

fi;

return g1;
end);

#################################################
InstallGlobalFunction( DefiningCharacterOfCyclotomicAlgebra, function(A)
local g1,d,m,n,i,chi,F,u,V,U,F1;

if Length(A)=2 then u:=1; else
if Length(A)>2 then
g1:=DefiningGroupOfCyclotomicAlgebra(A);
if Length(A)=4 then d:=A[4][1]; F1:=NF(A[3],[A[4][2]]); fi;
###Adjust to arbitrary number of generators ###
#if (Length(A)=5 and Length(A[4])=2) then
#   d:=A[4][1][1]*A[4][2][1];
#   F1:=NF(A[3],[A[4][1][2],A[4][2][2]]);
#fi;
#if (Length(A)=5 and Length(A[4])=3) then
#   d:=A[4][1][1]*A[4][2][1]*A[4][3][1];
#   F1:=NF(A[3],[A[4][1][2],A[4][2][2],A[4][3][2]]);
#fi;
if Length(A)=5 then
  V:=[];
  d:=1;
  for i in [1..Length(A[4])] do Add(V,A[4][i][2]); od;
  for i in [1..Length(A[4])] do d:=d*A[4][i][1]; od;
fi;

n:=Size(Irr(g1)) ;
m:=Trace(F1,Rationals,1);
U:=[];
for i in [1..n] do
chi:=Irr(g1)[n-i+1];
V:=ValuesOfClassFunction(chi);
F:=FieldByGenerators(V);
if V[1]/d in PositiveIntegers then
if Size(KernelOfCharacter(chi))=1 then
if FieldByGenerators(V)=F1 then
	Add(U,n-i+1);
fi;
fi;
fi;
od;
if Size(U)=m then
u:=U[1];
else
u:=U;
fi;

fi;
fi;

return u;
end);

##########################################
#  The next function was created to replace SimpleAlgebraByCharacterInfo
#  before it was fixed to work over larger fields.
##########################################
InstallGlobalFunction( SimpleComponentOfGroupRingByCharacter, function(F,G,n)
local R,chi,B;

R:=GroupRing(F,G);
if IsPosInt(n) then
  if HasOrdinaryCharacterTable(G) then
    chi:=Irr(G)[n];
  else
    Error("The group has no ordinary character table yet. To avoid randomisation errors, you should compute it first\n");
  fi;
elif IsCharacter(n) then
  chi:=n;
else
  Error("The third argument must be a character or its number\n");
fi;
B:=SimpleAlgebraByCharacterInfo(R,chi);

return B;
end);

#######################################################
#  Global Character Descent functions - this
#  adds as much Clifford theory as is possible over
#  the global field before local methods need to be
#  used.  This is needed because wedderga's local index
#  functions are designed for crossed products that
#  have been reduced in this way.
#######################################################
#######################################################
# These functions are intended to provide an enhancement
# to wedderga to extend its capabilities from small groups
# to some medium-sized groups. Given chi Irr(G)[n] and a
# cyclotomic field K, you want to determine the simple
# component of KG corresponding to chi.  The algorithm
# initially searches the irreducible characters of
# maximal subgroups M for a constituent phi for which
# (chi_M, phi) is coprime to chi(1) and K(phi)=K, and
# when it finds such a pair (M,phi) it replaces (G,chi)
# with (M,phi) and repeats.  This works because this
# condition implies the simple component of KM corresponding
# to phi is Brauer equivalent to the simple component of KG
# corresponding to chi. Once this process terminates it
# reverts to wedderga's functions for expressing the
# simple component.
#
# The main computational barrier to its effectiveness is
# for some larger groups have too many maximal subgroups.
# Since it must store them all one can get memory crashes.
# Another issue is that there are groups with characters
# have no global reduction with arbitrarily large size, or
# the algorithm terminates too quickly to a subgroup that
# is still too large for wedderga to handle.  In these
# unlucky situations this algorithm will be a waste of time.
#
# For calculating the simple component of KG corresponding
# to chi = Irr(G)[n], the command is
#
# SimpleComponentByCharacterDescent(K,G,n);
#
# For the Wedderburn Decomposition of KG the command is
#
# WedderburnDecompositionByCharacterDescent(K,G);
#
# (note that this differs from other wedderga functions since
# the input here is not the group ring). This will perform
# the calculation one character at a time, its performance
# is only saved because the maximal subgroup lattice won't be
# recalculated each time.
################################################
InstallGlobalFunction( CharacterDescent, function(F,G,n,e,H)
local chi,M,m,r,y,U,F1,D,y1,i,psi,C,j,m1,F2,y2,chi1,k,n1;

chi:=Irr(G)[n];
m:=chi[1];
r:=1;
y:=PrimitiveElement(F);
F1:=Field(Concatenation(chi,[y]));
D:=[r,F1,G,n];
y1:=PrimitiveElement(F1);
psi:=RestrictedClassFunction(chi,H);
C:=ConstituentsOfCharacter(psi);
for j in [1..Size(C)] do
  m1:=ScalarProduct(psi,C[j]);
  if Gcd(m1,e)=1 then
    F2:=Field(C[j]);
    y2:=PrimitiveElement(F2);
    if y2 in F1 then
      for k in [1..Size(Irr(H))] do
        if Irr(H)[k]=C[j] then
          n1:=k;
          r:=m/(ValuesOfClassFunction(C[j])[1]);
          D:=[r,F1,H,k];
          break;
        fi;
      od;
    fi;
  fi;
od;

return D;
end);
########################
InstallGlobalFunction( GlobalCharacterDescent, function(F,G,n)
local e,r,t,y,V,F1,R,M,m,n1,s,i,H,D;

e:=Size(G);
r:=1;
t:=0;
y:=PrimitiveElement(F);
V:=Concatenation(ValuesOfClassFunction(Irr(G)[n]),[y]);
F1:=Field(V);
R:=[r,F1,G,n];

n1:=Irr(G)[n][1];
if n1>1 then
while t=0 do
M:=ConjugacyClassesMaximalSubgroups(R[3]);
m:=Size(M);
i:=m;
for i in [1..m] do
  H:=Representative(M[m-i+1]);
  D:=CharacterDescent(F1,R[3],R[4],e,H);
  if Size(D[3])<Size(R[3]) then
    s:=n1/(Irr(D[3])[D[4]][1]);
    R:=[r*s,F1,D[3],D[4]];
    t:=1;
    break;
  fi;
od;
if i=m and t=0 then t:=1; fi;
od;
fi;

return R;
end);
########################################
InstallGlobalFunction( GaloisRepsOfCharacters, function(F1,G)
local U,U1,V,y,y1,T,n,k,i,j,i1,t,F,m;

U:=[1];
U1:=[];
y:=PrimitiveElement(F1);
T:=Irr(G);
n:=Size(Irr(G));
i:=1;
for i in [2..n] do
 V:=ValuesOfClassFunction(T[i]);
 F:=Field(V);
 y1:=PrimitiveElement(F);
 if y1 in F1 then
  Add(U,i);
 else
  Add(U1,i);
 fi;
od;
if not(U1=[]) then
Add(U,U1[1]);
for i in [2..Length(U1)] do
 F:=Field(Union(ValuesOfClassFunction(T[U1[i]]),[y]));
 k:=Conductor(F);
 t:=1;
 for j in [2..k-1] do
   if Gcd(j,k)=1 then
   if (GaloisCyc(y,j)=y) then
   for i1 in [1..i-1] do
    if GaloisCyc(ValuesOfClassFunction(T[U1[i1]]),j)=ValuesOfClassFunction(T[U1[i]]) then
    t:=0;
    break;
    fi;
   od;
   fi;
   fi;
 od;
 if t=1 then Add(U,U1[i]); fi;
od;
fi;

Sort(U);

return U;
end);
###########################
InstallGlobalFunction( SimpleComponentByCharacterDescent, function(F,G,n)
local R,n1,n2,r,S,S0,T;

n1:=Size(G);
R:=GlobalCharacterDescent(F,G,n);
n2:=Size(R[3]);
r:=R[1];
while n2<n1 do
  n2:=Size(R[3]);
  R:=GlobalCharacterDescent(R[2],R[3],R[4]);
  n1:=Size(R[3]);
  r:=r*R[1];
od;
T:=SimpleComponentOfGroupRingByCharacter(R[2],R[3],R[4]);
T[1]:=r*T[1];

T:=GlobalSplittingOfCyclotomicAlgebra(T);

return T;
end);
###########################
InstallGlobalFunction( WedderburnDecompositionByCharacterDescent, function(F,G)
local R,S,y,n,U,F1,T;

R:=[];
S:=GaloisRepsOfCharacters(F,G);
y:=PrimitiveElement(F);
for n in S do
  U:=Union(ValuesOfClassFunction(Irr(G)[n]),[y]);
  F1:=Field(U);
  T:=SimpleComponentByCharacterDescent(F1,G,n);
  Add(R,T);
od;

return R;
end);
######################

################################
# Given a simple component of a rational group algebra whose
# "WedderburnDecompositionInfo" output has 4 terms, the next
# three functions compute its indices at odd primes, infinity,
# and 2.
################################
InstallGlobalFunction( LocalIndexAtOddP, function(A,q)
local m,F,n,a,b,c,n1,e,f,h,f1,e1,k;

m:=1;
F:=A[2];
a:=A[4][1];
b:=A[4][2];
c:=A[4][3];
n:=Lcm(Conductor(F),A[3]);
n1:=PDashPartOfN(n,q);
##########################
e:=RamificationIndexAtP(F,n,q);
if e>1 and c>0 and A[3]/q in PositiveIntegers then
f:=ResidueDegreeAtP(Rationals,n,q);
h:=ResidueDegreeAtP(F,n,q);
f1:=f/h;
e1:=Gcd(q^f1-1,e);
k:=(q^f1-1)/e1;
while not(k/(Order(E(A[3])^(c*m))) in PositiveIntegers) do
  m:=m+1;
od;
fi;

return m;
end);




####################################
# Cyclotomic reciprocity calculation:
#U:=[1];
#if n1 > 1 then
#i:=1;
#while not(PowerMod(q,i,n1)=1) do
#Add(U,PowerMod(q,i,n1));
#i:=i+1;
#od;
#fi;
#U1:=[];
#n2:=PDashPartOfN(A[3],q);
#for u in U do
#Add(U1,u mod n2);
#od;
#g:=1;
#while not(PowerMod(b,g,n2) in U1) do
#g:=g+1;
#od;
#h:=OrderMod(b^g,n2);
#e:=OrderMod(b^g,A[3])/h;
#####################################
# Now g, h, and e are the splitting,
# residue degree and ramification of
# Q(E(n))/F at q.
#####################################
#if e>1 and c>0 and A[3]/q in PositiveIntegers then
#####################################
# Compute residue degree f of CF(n) at q
#f:=1;
#while not(n1=1 or PowerMod(q,f,n1)=1) do
#f:=f+1;
#od;
#########################
# Now F_q contains E(q^(f/h)-1).  We find the least
# power m of E(A[3])^c that lies in the group generated
# by E(q^(f/h)-1)^e.
#########################
#f1:=f/h;
#e1:=Gcd(q^f1-1,e);
#k:=(q^f1-1)/e1;
#while not(k/(Order(E(A[3])^(c*m))) in PositiveIntegers) do
#  m:=m+1;
#od;
#fi;

#return m;
#end);


###############################
# For the computation of the index of a cyclic
# cyclotomic algebra at infinity, we determine the
# nature of the algebra as a quadratic algebra
# over the reals.
###############################
InstallGlobalFunction( LocalIndexAtInfty, function(A)
local m,n,s,n1;

m:=1;
n:=A[3];
if n>2 then
  s:=ANFAutomorphism(A[2],-1);
  if s=ANFAutomorphism(A[2],1) then
     if E(n)^A[4][3]=-1 then
      m:=2;
    fi;
  fi;
fi;

return m;
end);

###############################
# For the local index at 2, we detect if the cyclic
# cyclotomic algebra will be of nonsplit quaternion type
# over the 2-adics
###############################
InstallGlobalFunction( LocalIndexAtTwo, function(A)
local n,m,a,K,b,c,f1,f,e1,e,h,g,n2,n3,n4,i,u,U,U1;

m:=1;
a:=A[4][1];
if A[3]/4 in PositiveIntegers and IsEvenInt(a) then
n:=Lcm(Conductor(A[2]),A[3]);
f1:=ResidueDegreeAtP(A[2],n,2);
f:=ResidueDegreeAtP(Rationals,n,2);
e1:=RamificationIndexAtP(A[2],n,2);
e:=RamificationIndexAtP(Rationals,n,2);
if IsOddInt(f/f1) and IsOddInt(e/e1) then
#K:=PSplitSubextension(A[2],n,2);
#if IsOddInt(Trace(K,A[2],1)) then
b:=A[4][2];
c:=A[4][3];
n2:=PPartOfN(n,2);
n4:=PPartOfN(A[3],2);
  if E(n2)^b=E(n2)^(-1) and E(n4)^c=-1 then
   m:=2;
  fi;
fi;
fi;

return m;
end);

#################################
#n1:=PDashPartOfN(n,2);
#f:=OrderMod(2,n1);
#n2:=PPartOfN(n,2);
#e:=Phi(n2);

#U:=[1];
#i:=1;
#while not(PowerMod(2,i,n1) in U) do
#Add(U,PowerMod(2,i,n1));
#i:=i+1;
#od;
#U1:=[];
#n3:=PDashPartOfN(A[3],2);
#n4:=PPartOfN(A[3],2);
#for u in U do
#Add(U1,u mod n3);
#od;

#g:=1;
#while not(PowerMod(b,g,n3) in U1) do
#g:=g+1;
#od;
#h:=OrderMod(b^g,n3);
#e1:=OrderMod(b^g,A[3])/OrderMod(b^g,n3);
# if e1>1 and IsOddInt(e*f/e1*h) then
# n2:=PPartOfN(n,2);
#  if E(n2)^(b^g)=E(n2)^(-1) and E(n4)^c=-1 then
#   m:=2;
#  fi;
# fi;
#fi;
#
#return m;
#end);
##############################
# Given a group G and a simple component A whose
# WedderburnDecompositionInfo in wedderga has length 4,
# this program gives the list of local indices at
# all primes relevant to the rational Schur index
###############################
InstallGlobalFunction( LocalIndicesOfCyclicCyclotomicAlgebra, function(A)
local n,S,s,i,L,l,q,L1;

L:=[];
if A[4][3]=0 then
  L1:=[];
else
S:=AsSet(FactorsInt(A[3]));
s:=Size(S);
for i in [1..s] do
  if S[i]=2 then
  l:=LocalIndexAtTwo(A);
  L[i]:=[];
  L[i][1]:=2;
  L[i][2]:=l;
  else
  q:=S[i];
  l:=LocalIndexAtOddP(A,q);
  L[i]:=[];
  L[i][1]:=q;
  L[i][2]:=l;
  fi;
od;

l:=LocalIndexAtInfty(A);
L[s+1]:=[];
L[s+1][1]:= infinity;
L[s+1][2]:=l;

L1:=[];

s:=Size(L);
for i in [1..s] do
if L[i][2]>1 then
Add(L1,L[i]);
fi;
od;
fi;

return L1;
end);




##########################################
InstallGlobalFunction( IsDyadicSchurGroup, function(G)
local d,j,P,P1,V0,c,g,g1,z,V1,v1,V2,q,s,n1,r,i,y,p1,x,U,V,P2,L;

d:=false;
j:=0;
P:=SylowSubgroup(G,2);
q:=Size(G)/Size(P);
if IsPrimeInt(q) then
U:=SylowSubgroup(G,q);
V:=Centralizer(P,U);

# First look for G of type (Q_8,q)
if IdSmallGroup(V)=[8,4] then
 V1:=Centralizer(P,V);
 L:=UnionSet(Elements(V),Elements(V1));
 P2:=GroupByGenerators(L);
 if P=P2 then
  d:=true;
  j:=1;
 fi;
else
# Checks if P is of type (QD,q)
P1:=DerivedSubgroup(P);
 s:=LogInt(Size(P)/Size(P1),2)-1;
 if s=LogInt(PPartOfN(OrderMod(2,q),2),2) then
 if not(IsAbelian(V)) then
 if Size(P)/Size(V)=2^s then
 d:=true;
 fi;
 fi;
 fi;
fi;

fi;

return d;
end);

##########################################
InstallGlobalFunction( LocalIndexAtInftyByCharacter, function(F,G,n)
local m,T,v2,a,pos;

if IsPosInt(n) then
  if HasOrdinaryCharacterTable(G) then
    pos:=n;
  else
    Error("The group has no ordinary character table yet. To avoid randomisation errors, you should compute it first\n");
  fi;
elif IsCharacter(n) then
  pos := Position( Irr(G), n );
else
  Error("The fourth argument must be a character or its number\n");
fi;

m:=1;
T:=CharacterTable(G);
v2:=Indicator(T,2)[pos];
a:=PrimitiveElement(F);
if GaloisCyc(a,-1)=a then
if v2=-1 then
m:=2;
fi;
fi;

return m;
end);

###########################################
InstallGlobalFunction( FinFieldExt, function(F,G,p,n,n1)
local chi,V,Y,h,a,m1,d1,L,i,z,l,m,K,B,d,M,C,D,b,j,F1,M1,M2,psi,U,k,F2,t;

if IsPosInt(n) then
  if HasOrdinaryCharacterTable(G) then
    chi:=Irr(G)[n];
  else
    Error("The group has no ordinary character table yet. To avoid randomisation errors, you should compute it first\n");
  fi;
elif IsCharacter(n) then
  chi:=n;
else
  Error("The fourth argument must be a character or its number\n");
fi;

if IsPosInt(n1) then
  if HasOrdinaryCharacterTable(G) and IsBound( ComputedBrauerTables( CharacterTable(G) )[p] ) then
    psi:=Irr( BrauerTable(G,p) )[n1];
  else
    Error("The group has no Brauer character table for p=", p, " yet. To avoid randomisation errors, you should compute it first\n");
  fi;
elif IsCharacter(n1) and n1 in Irr( BrauerTable(G,p) ) then
  psi:=n1;
else
  Error("The fifth argument must be a Brauer character at ", p, " or its number\n");
fi;

V:=ValuesOfClassFunction(chi);
Y:=OrdersClassRepresentatives( CharacterTable(G) );
h:=Size(Y);
a:=PrimitiveElement(F);
m1:=PDashPartOfN(Conductor(a),p);
for i in [1..m1] do if (m1=1 or PowerMod(p,i,m1)=1) then d1:=i; break; fi; od;

L:=[];
for i in [1..h] do if Gcd(Y[i],p) = 1 then Add(L,V[i]); fi; od;
l:=Size(L);
m:=Conductor(L);
K:=CF(m);
B:=Basis(K);
for i in [1..m] do if (m=1 or PowerMod(p,i,m)=1) then d:=i; break; fi; od;
z:=Z(p^d)^((p^d-1)/m);
M:=[];
D:=[];
for i in [1..Size(B)] do
for j in [1..m] do
if B[i]=E(m)^j then
D[i]:=j;
fi;
od;
od;

for i in [1..l] do
  C:=Coefficients(B,L[i]);
  b:=0;
for j in [1..Size(B)] do
  b:=b+C[j]*z^(D[j]);
od;
  M[i]:=b;
od;
M1:=UnionSet(M,[Z(p^d1)]);
F1:=FieldByGenerators(M1);

U:=ValuesOfClassFunction(psi);
m:=Conductor(U);
K:=CF(m);
B:=Basis(K);
for i in [1..m] do if (m=1 or PowerMod(p,i,m)=1) then d:=i; break; fi; od;
z:=Z(p^d)^((p^d-1)/m);
M:=[];
D:=[];
for i in [1..Size(B)] do
for j in [1..m] do
if B[i]=E(m)^j then
D[i]:=j;
break;
fi;
od;
od;

for i in [1..l] do
  C:=Coefficients(B,U[i]);
  b:=0;
for j in [1..Size(B)] do
  b:=b+C[j]*z^(D[j]);
od;
  M[i]:=b;
od;

M2:=UnionSet(M,M1);
F2:=FieldByGenerators(M2);
t:=LogInt(Size(F2),Size(F1));

return t;
end);

##############################################
# Oct 2014 New Defect Group Functions
##############################################
InstallGlobalFunction( DefectGroupOfConjugacyClassAtP, function(G,c,p)
local C,g1,H,D;

C:=ConjugacyClasses(G);
g1:=Representative(C[c]);
H:=Centralizer(G,g1);
D:=SylowSubgroup(H,p);

return D;
end);
#############################
InstallGlobalFunction( DefectGroupsOfPBlock, function(G,n,p)
local D1,D2,r,U,U1,T,chi,C,c,i,m,h1,a1,a2,b1,A1,D;

if IsPosInt(n) then
  if HasOrdinaryCharacterTable(G) then
    chi:=Irr(G)[n];
  else
    Error("The group has no ordinary character table yet. To avoid randomisation errors, you should compute it first\n");
  fi;
elif IsCharacter(n) then
  chi:=n;
else
  Error("The third argument must be a character or its number\n");
fi;

T:=CharacterTable(G);
C:=ConjugacyClasses(G);
c:=Size(C);
U:=[];
U1:=[];
for i in [1..c] do
  m:=OrdersClassRepresentatives(T)[i];
  if not(m mod p = 0 mod p) then
    AddSet(U,i);
    AddSet(U1,i);
  fi;
od;

for i in U1 do
  h1:=Size(C[i]);
  a1:=chi[i];
  b1:=chi[1];
  A1:=(h1*a1)/b1;
  r:=PDashPartOfN(Size(G),p);
  a2:=Norm(r*A1);
  if not(a2 in Integers) or (a2/p in Integers) then
     RemoveSet(U,i);
  fi;
od;

D2:=[];
for i in U do
  D:=DefectGroupOfConjugacyClassAtP(G,i,p);
  AddSet(D2,D);
od;

if Length(D2)>1 then
  D1:=D2[1];
for i in [2..Size(D2)] do
  if Size(D2[i])<Size(D1) then
    D1:=D2[i];
  fi;
od;
else
D1:=D2[1];
fi;

D:=ConjugacyClassSubgroups(G,D1);

return D;
end);
####################
InstallGlobalFunction( DefectOfCharacterAtP, function(G,n,p)
local D1,D,q,d;

D1:=DefectGroupsOfPBlock(G,n,p);
D:=Representative(D1);
q:=Size(D);
d:=LogInt(q,p);

return d;
end);

###################################################
#InstallGlobalFunction( PossibleDefectGroups, function(G,n,p)
#local S,U,U0,U1,Q,Q1,Q2,i,j,I,T,b,k,d,H,a;
#
#T:=CharacterTable(G);
#S:=T mod p;
#b:=BlocksInfo(S);
#for j in [1..Size(b)] do
#if n in b[j].ordchars then
#  k:=b[j].modchars[1];
#  d:=b[j].defect;
#  break;
#fi;
#od;
#Q:=SylowSubgroup(G,p);
#U:=[];
#U0:=[];
#U1:=[];
#if Size(Q)>p^d then
#H:=ConjugacyClasses(G);
#for j in [2..Size(H)] do
#  if Gcd(OrdersClassRepresentatives(T)[j],p)=1 then
#     a:=Elements(H[j])[1];
#     Q1:=Intersection(Q,Q^a);
#     if Size(Q1)=p^d then
#        AddSet(U0,ConjugacyClassSubgroups(G,Q1));
#     fi;
#     Q2:=SylowSubgroup(Centralizer(G,a),p);
#     if Size(Q2)=p^d then
#        AddSet(U1,ConjugacyClassSubgroups(G,Q2));
#     fi;
#  fi;
#od;
#U:=Intersection(U0,U1);
#else
#AddSet(U,ConjugacyClassSubgroups(G,Q));
#fi;
#
#return U;
#end);

##########################################
InstallGlobalFunction( LocalIndexAtPByBrauerCharacter, function(F,G,n,p)
local chi,n1,V,a,V1,C,m1,b,j,k,u,t,T,S,U,f,m2,n0,K0,d0,F1,K1,d1;

if IsPosInt(n) then
  if HasOrdinaryCharacterTable(G) then
    chi:=Irr(G)[n];
    n1:=n;
  else
    Error("The group has no ordinary character table yet. To avoid randomisation errors, you should compute it first\n");
  fi;
elif IsCharacter(n) then
  chi:=n;
  n1:=Position(Irr(G),chi);
else
  Error("The third argument must be a character or its number\n");
fi;

V:=ValuesOfClassFunction(chi);
a:=PrimitiveElement(F);
V1:=Union(V,[a]);
F1:=FieldByGenerators(V1);
C:=FieldByGenerators(V);
m1:=[];
m1[1]:=1;
m1[2]:="DGisCyclic";
T:=CharacterTable(G);
S:=T mod p;
b:=BlocksInfo(S);
for j in [1..Size(b)] do
if n1 in b[j].ordchars then
  k:=b[j].modchars[1];
  break;
fi;
od;
#####################
# Adapted to new defect group function
#####################
#U:=PossibleDefectGroups(G,n,p);
#f:=0;
#for u in [1..Size(U)] do
# if not(IsCyclic(Elements(U[u])[1])) then
#   f:=f+1;
# fi;
#od;
#if not(f=0) then
#  if f<Size(U) then
#   m1[2]:="DGmaybeCyclic";
#  else
#   m1[2]:="DGnotCyclic";
#  fi;
#fi;
U:=DefectGroupsOfPBlock(G,n,p);
if not(IsCyclic(Representative(U))) then
  m1[2]:="DGnotCyclic";
fi;

####################################

t:=FinFieldExt(C,G,p,n,k);
if t>1 then
m1[1]:=t;
fi;

m2:=m1;
if m2[2]="DGisCyclic" then
m1:=m2[1];
a:=PrimitiveElement(F);
V1:=Union(V,[a]);
n0:=Conductor(V1);
K0:=PSplitSubextension(C,n0,p);
d0:=Trace(CF(n0),K0,1);
F1:=FieldByGenerators(V1);
K1:=PSplitSubextension(F1,n0,p);
d1:=Trace(CF(n0),K1,1);
m1:=m1/Gcd(m1,d0/d1);
fi;

return m1;
end);

###########################################
InstallGlobalFunction( LocalIndexAtOddPByCharacter, function(F,G,n,p)
local m,B,K,B1,g,n1;

m:=1;
B:=SimpleComponentOfGroupRingByCharacter(F,G,n);

if Length(B)=2 then
m:=1;
fi;
if Length(B)=4 then
m:=LocalIndexAtOddP(B,p);
fi;

if Length(B)=5 then
K:=PSplitSubextension(F,B[3],p);
B1:=SimpleComponentOfGroupRingByCharacter(K,G,n);
g:=DefiningGroupAndCharacterOfCyclotAlg(B1);
m:=LocalIndexAtPByBrauerCharacter(K,g[1],g[2],p);
#g:=DefiningGroupOfCyclotomicAlgebra(B1);
#n1:=DefiningCharacterOfCyclotomicAlgebra(B1);
#m:=LocalIndexAtPByBrauerCharacter(K,g,n1,p);
fi;

return m;
end);

###########################################
InstallGlobalFunction( LocalIndexAtTwoByCharacter, function(F,G,n)
local m,chi,g,g1,B1,W,i,a,B,K,V,V1,a1,F0,F1,n0,n1,n01,n02,n11,n12,f,f0,f1,m2;

if IsPosInt(n) then
  if HasOrdinaryCharacterTable(G) then
    chi:=Irr(G)[n];
  else
    Error("The group has no ordinary character table yet. To avoid randomisation errors, you should compute it first\n");
  fi;
elif IsCharacter(n) then
  chi:=n;
else
  Error("The third argument must be a character or its number\n");
fi;

m2:=1;
m:=0;
B:=SimpleComponentOfGroupRingByCharacter(F,G,n);

if Length(B)=2 then
m2:=1;
fi;
if Length(B)=4 then
m2:=LocalIndexAtTwo(B);
fi;

if Length(B)=5 then
K:=PSplitSubextension(F,B[3],2);
B1:=SimpleComponentOfGroupRingByCharacter(K,G,n);
g:=DefiningGroupAndCharacterOfCyclotAlg(B1);
m2:=LocalIndexAtPByBrauerCharacter(K,g[1],g[2],2);
#g:=DefiningGroupOfCyclotomicAlgebra(B1);
#n1:=DefiningCharacterOfCyclotomicAlgebra(B1);
#m2:=LocalIndexAtPByBrauerCharacter(F,g,n1,2);
if not(m2 in Integers) then
m:=1;
  if IsDyadicSchurGroup(g[1]) then
  m:=2;
  V:=ValuesOfClassFunction(chi);
  F0:=FieldByGenerators(V);
  F1:=B1[2];
    if not(F0=F1) then
      if E(4) in F1 then
        m:=1;
      else
        n0:=Conductor(F0);
        n02:=PPartOfN(n0,2);
        n1:=Conductor(F1);
        n12:=PPartOfN(n1,2);
          if not(n02=n12) then
            m:=1;
          else
            n11:=PDashPartOfN(n1,2);
            f1:=OrderMod(2,n1);
            n01:=PDashPartOfN(n0,2);
            f0:=OrderMod(2,n0);
            f:=f1/f0;
              if (f/2 in PositiveIntegers) then
                m:=1;
              fi;
          fi;
       fi;
    fi;
  fi;
fi;
fi;

if m>0 then m2:=m; fi;

return m2;
end);

#############################################
InstallGlobalFunction( LocalIndicesOfCyclotomicAlgebra, function(A)
local L,F,l,d,G,n,m0,m2,m,P,p,l1,l2,i,L1;

##################
# bugfix lines (20/03/2020)
##################
A:=GlobalSplittingOfCyclotomicAlgebra(A);
l:=Length(A);
if l>4 then
  l2:=Length(A[4]);
  l1:=l2-1;
  l:=Length(A);
  while l>4 and l1<l2 do
    l2:=Length(A[4]);
    G:=DefiningGroupOfCyclotomicAlgebra(A);
    n:=DefiningCharacterOfCyclotomicAlgebra(A);
    A:=SimpleComponentByCharacterDescent(A[2],G,n);
    l:=Length(A);
    if l>4 then l1:=Length(A[4]); fi;
  od;
fi;
##################

L:=[];
L1:=[];
F:=A[2];
l:=Length(A);

if l=5 then
  d:=DefiningGroupAndCharacterOfCyclotAlg(A);
  G:=d[1];
  n:=d[2];
  #G:=DefiningGroupOfCyclotomicAlgebra(A);
  #n:=DefiningCharacterOfCyclotomicAlgebra(A);
  m0:=LocalIndexAtInftyByCharacter(F,G,n);
  Add(L1,[infinity,m0]);

  P:=AsSet(Factors(Size(G)));
  if P[1]=2 then
    m2:=LocalIndexAtTwoByCharacter(F,G,n);
    Add(L1,[2,m2]);
  fi;

  P:=Difference(P,[2]);
  if Size(P)>0 then
    for i in [1..Size(P)] do
      p:=P[i];
      m:=LocalIndexAtOddPByCharacter(F,G,n,p);
      Add(L1,[p,m]);
    od;
  fi;

  l1:=Size(L1);
  for i in [1..l1] do
   if not(L1[i][2]=1) then
     Add(L,L1[i]);
   fi;
  od;
fi;

if (l=4 and not(A[4][3]=0)) then
  L:=LocalIndicesOfCyclicCyclotomicAlgebra(A);
fi;

return L;
end);

############################################
InstallGlobalFunction( RootOfDimensionOfCyclotomicAlgebra, function(A)
local d,i;

if Length(A)<4 then
d:=A[1];
fi;

if Length(A)=4 then
d:=A[1]*A[4][1];
fi;

if Length(A)=5 then
d:=A[1];
for i in [1..Length(A[4])] do
d:=d*A[4][i][1];
od;
fi;

return d;
end);

###########################################
# Calculates the least common multiple of the list of
# local indices
###########################################

InstallGlobalFunction( GlobalSchurIndexFromLocalIndices, function(L)
local l,m,i;

l:=Length(L);
m:=1;
if l>0 then
  m:=L[1][2];
fi;
if l>1 then
for i in [2..l] do
  m:=Lcm(m,L[i][2]);
od;
fi;

return m;
end);

###########################################
InstallGlobalFunction( CyclotomicAlgebraWithDivAlgPart, function(A)
local L,m,d,B,D;

L:=LocalIndicesOfCyclotomicAlgebra(A);
m:=GlobalSchurIndexFromLocalIndices(L);
d:=RootOfDimensionOfCyclotomicAlgebra(A);
if m>1 then
D:=rec(DivAlg:=true, Center:=A[2], SchurIndex:=m, LocalIndices:=L);
B:=[d/m,D];
else
B:=[d,A[2]];
fi;

return B;
end);

###############################################
# Main function for obtaining the Wedderburn decomposition
# for R = GroupRing(F,G) with division algebra parts identified
# in terms of local indices
###############################################
InstallGlobalFunction( WedderburnDecompositionWithDivAlgParts, function(R)
local W,w,i,W1;

W:=WedderburnDecompositionInfo(R);
w:=Size(W);
W1:=[];
for i in [1..w] do
if Length(W[i]) < 4 then
  Add(W1,W[i]);
else
W1[i]:=CyclotomicAlgebraWithDivAlgPart(W[i]);
fi;
od;

return W1;
end);

#############################
# Given a Schur algebra output from "wedderga" with 5 terms
# that decomposes as the tensor product of two generalized
# quaternion algebras, the first function determines this
# tensor decomposition.
# #############################
InstallGlobalFunction( DecomposeCyclotomicAlgebra, function(A)
local B,B1,m1,n,m,d,c,z,r,s,t,u,v,u1,i,j,b,F,w;

B:=[];
B[1]:=[];
B[2]:=[];
n:=A[3];

if (Length(A)=5 and A[4][1][1]=2 and A[4][2][1]=2) then

if not(A[5][1][1]=0) then
  d:=A[5][1][1];
  z:=E(n)^d;

if z=-1 and not(E(4) in A[2]) then
  if E(4)^A[4][2][2]=E(4) then
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[E(n)^A[4][1][3]] ];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[E(4)^2*E(n)^A[4][2][3]]];
  else
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[E(4)*GaloisCyc(E(4),A[4][1][2])*E(n)^A[4][1][3]]];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[E(n)^A[4][2][3]] ];
  fi;
fi;

if z=-E(4) and not(E(4) in A[2]) then
  if E(4)^A[4][2][2]=E(4) then
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[E(n)^A[4][1][3]] ];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[(1-E(4))^2*E(n)^A[4][2][3]]];
  else
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[(1+E(4))*GaloisCyc((1+E(4)),A[4][1][2])*E(n)^A[4][1][3]]];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[ E(n)^A[4][2][3]] ];
  fi;
fi;

if z=E(4) and not(E(4) in A[2]) then
  if E(4)^A[4][2][2]=E(4) then
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[E(n)^A[4][1][3]] ];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[(1+E(4))^2*E(n)^A[4][2][3]]];
  else
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[(1-E(4))*GaloisCyc((1-E(4)),A[4][1][2])*E(n)^A[4][1][3]]];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[ E(n)^A[4][2][3]] ];
  fi;
fi;

if not(z^4=1) and not(z in A[2]) then
m:=Order(z);
u:=E(m);

if E(m)^A[4][2][2]=E(m) then
r:=A[4][1][2]^-1 mod n;
w:=z^-1;
fi;

if E(m)^A[4][1][2]=E(m) then
r:=A[4][2][2]^-1 mod n;
w:=GaloisCyc(z,A[4][2][2])^-1;
fi;
t:=DescriptionOfRootOfUnity(w)[2];
v:=w*u;

for i in [1..n] do
v:=v^r;
if v=E(m) then
   break;
else
u:=u+v;
v:=E(m)^t*v;
fi;
od;

  if E(m)^A[4][2][2]=E(m) then
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[E(n)^A[4][1][3]] ];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[u*GaloisCyc(u,A[4][2][2])*E(n)^A[4][2][3]]];
  else
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[u*GaloisCyc(u,A[4][1][2])*E(n)^A[4][1][3]]];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[ E(n)^A[4][2][3]] ];
  fi;
fi;

for i in [1,2] do
  if not(IsSubset(B[i][2],B[i][1])) then
  c:=PrimitiveElement(B[i][1]);
  d:=PrimitiveElement(B[i][2]);
  B[i][2]:=FieldByGenerators([c,d]);
  fi;
od;

#  c:=[];
#  for j in [1,2] do
#  d:=Conductor(B[i][j]);
#  F:=CF(d);
#  if not(IsCyclotomicField(B[i][j])) then
#    c[j]:=Trace(F,B[i][1],E(d));
#  else
#    c[j]:=E(d);
#  fi;
#  od;
#  B[i][2]:=FieldByGenerators(c);

  B1:=B;

else
  B[1][1]:=A[2];
  B[1][2]:=NF(n,[1,A[4][2][2]]);
  B[1][3]:=[E(n)^A[4][1][3]];
  B[2][1]:=A[2];
  B[2][2]:=NF(n,[1,A[4][1][2]]);
  B[2][3]:=[E(n)^A[4][2][3]];
  B1:=B;
fi;
fi;

if Length(A)=5 and (A[4][1][1]>2 or A[4][2][1]>2) then

  if not(A[5][1][1]=0) then
  d:=A[5][1][1];
  z:=E(n)^d;
  m:=Order(z);
  u:=E(m);

if E(m)^A[4][2][2]=E(m) then
r:=A[4][1][2]^-1 mod n;
w:=z^-1;
else
r:=A[4][2][2]^-1 mod n;
w:=GaloisCyc(z,A[4][2][2])^-1;
fi;
t:=DescriptionOfRootOfUnity(w)[2];
v:=w*u;

for i in [1..n] do
v:=v^r;
if v=E(m) then
   break;
else
u:=u+v;
v:=E(m)^t*v;
fi;
od;

if E(m)^A[4][2][2]=E(m) then
u1:=u;
for i in [1..(A[4][2][1]-1)] do
  u1:=u1*GaloisCyc(u1,A[4][2][2]);
od;
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[E(n)^A[4][1][3]] ];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[u1*E(n)^A[4][2][3]]];
  if not(B[2][3][1] in A[2]) then B[2][4]:="fail"; fi;
else
u1:=u;
for i in [1..(A[4][1][1]-1)] do
  u1:=u1*GaloisCyc(u1,A[4][1][2]);
od;
  B[1]:=[A[2],NF(n,[1,A[4][2][2]]),[u1*E(n)^A[4][1][3]]];
  B[2]:=[A[2],NF(n,[1,A[4][1][2]]),[ E(n)^A[4][2][3]] ];
  if not(B[1][3][1] in A[2]) then B[1][4]:="fail"; fi;
fi;

for i in [1,2] do
  if not(IsSubset(B[i][2],B[i][1])) then
  c:=PrimitiveElement(B[i][1]);
  d:=PrimitiveElement(B[i][2]);
  B[i][2]:=FieldByGenerators([c,d]);
  fi;
od;

B1:=B;

else
  B[1][1]:=A[2];
  B[1][2]:=NF(n,[1,A[4][2][2]]);
  B[1][3]:=[E(n)^A[4][1][3]];
  B[2][1]:=A[2];
  B[2][2]:=NF(n,[1,A[4][1][2]]);
  B[2][3]:=[E(n)^A[4][2][3]];
  B1:=B;
fi;
fi;

return B1;
end);

################################################
# The next few functions allow conversions between
# cyclic algebras and quaternion algebras.
########################################################
InstallGlobalFunction( ConvertQuadraticAlgToQuaternionAlg, function(A)
local d,t,n,i,B;

n:=Conductor(A[2]);
i:=0;
t:=Trace(A[2],A[1],1);

if t=2 then
for d in [1..n] do
  if Sqrt(d) in A[2] and not(Sqrt(d) in A[1]) then
     i:=d;
     break;
  fi;
  if Sqrt(-d) in A[2] and not(Sqrt(-d) in A[1]) then
     i:=-d;
     break;
  fi;
od;
fi;

if not(i=0) then
  B:=QuaternionAlgebra(A[1],i,A[3][1]);
else
  B:="fail";
fi;

return B;
end);

#####################################################
InstallGlobalFunction( ConvertCyclicCyclotomicAlgToCyclicAlg, function(A)
local n,a,K,B;

if Length(A)=4 then
n:=A[3];
a:=PrimitiveElement(A[2]);
K:=FieldByGenerators([a,E(n)]);
B:=[A[2],K,[E(n)^A[4][3]]];
else
B:="fail";
fi;
return [A[1],B];
end);

###############################################
InstallGlobalFunction( ConvertCyclicAlgToCyclicCyclotomicAlg, function(A)
local F,K,n,i,j,M,k,l,m,B;

F:=A[1];
K:=A[2];
B:="fails";
if IsCyclotomicField(K) then
n:=Conductor(K);
if IsOddInt(n) then n:=2*n; fi;
if A[3][1]^n=1 then

k:=0;
for i in [1..n-1] do
if Gcd(i,n) = 1 then
m:=OrderMod(n,i);
M:=[];
for j in [1..m] do
AddSet(M,i^j mod n);
od;
if F=NF(n,M) then k:=i; break; fi;
fi;
od;

if k>0 then
m:=Order(ANFAutomorphism(K,k));
for i in [0..n] do
if E(n)^i=A[3][1] then l:=i; break; fi;
od;
B:=[1,F,n,[m,k,l]];
fi;

fi;
fi;

return B;
end);

#####################################################

#################################################
InstallGlobalFunction( ConvertQuaternionAlgToQuadraticAlg, function(A)
local F,K,B,b,d,d1,i,a;

d:=[];
b:=Elements(Basis(A));
for i in [1..4] do
  if b[i]=Identity(A) then
    d[i]:=0;
  else
   d[i]:=Sum(Coefficients(Basis(A),b[i]^2));
  fi;
od;

d1:=[];
for i in [1..4] do
  if not(d[i]=0) then
    Add(d1,d[i]);
  fi;
od;

Sort(d1);
F:=LeftActingDomain(A);
a:=PrimitiveElement(F);
if not(d1[1]+d1[3]<0) then
K:=FieldByGenerators([a,Sqrt(d1[1])]);
B:=[F,K,[d1[2]]];
else
if d1[3]<0 then
K:=FieldByGenerators([a,Sqrt(d1[2])]);
B:=[F,K,[d1[3]]];
else
K:=FieldByGenerators([a,Sqrt(d1[3])]);
B:=[F,K,[d1[2]]];
fi;
fi;

return B;
end);

##########################################
# The next few functions allow one to compute the
# local indices of rational quaternion algebras,
# and determine if it is a division algebra.
# The first one computes the local index of the
# symbol algebra (p,q) over Q when p and q are -1
# or a prime.  Warning: It will not work when p or
# q are other integers, and it does not check this fact.
##########################################
InstallGlobalFunction( LocalIndicesOfRationalSymbolAlgebra, function(a,b)
local p,q,L;

p:=a;
q:=b;
L:=[];
if p < q then
  p:=b;
  q:=a;
fi;

if p=-1 then L:=[[infinity,2],[2,2]]; fi;
if p=2 then L:=[]; fi;
if p>2 then
  if q=-1 then
    if Legendre(q,p)=-1 then L:=[[2,2],[p,2]]; else L:=[]; fi;
  fi;
  if q=2 then
    if Legendre(2,p)=-1 then L:=[[2,2],[p,2]]; else L:=[]; fi;
  fi;
  if p>q and q>2 then
    if Legendre(q,p)=-1 then
       if Legendre(p,q)=-1 then L:=[[q,2],[p,2]]; else L:=[[2,2],[p,2]]; fi;
    else
       if Legendre(p,q)=-1 then L:=[[2,2],[q,2]]; fi;
    fi;
  fi;
  if p=q then
    if Legendre(-1,p)=-1 then L:=[[2,2],[p,2]]; fi;
  fi;
fi;
return L;
end);

##################################################
InstallGlobalFunction( LocalIndicesOfTensorProductOfQuadraticAlgs, function(L,M)
local i,j,m,S,L1;

S:=[];
L1:=[];

if L=[] then L1:=M; else if M=[] then L1:=L;
else

for i in [1..Length(L)] do
AddSet(S,L[i][1]);
od;
for j in [1..Length(M)] do
AddSet(S,M[j][1]);
od;

for i in [1..Length(S)] do
m:=1;
for j in [1..Length(L)] do
if L[j][1]=S[i] then
 m:=(m+2) mod 4;
fi;
od;
for j in [1..Length(M)] do
if M[j][1]=S[i] then
 m:=(m+2) mod 4;
fi;
od;
if m>1 then
Add(L1,[S[i],2]);
fi;
od;

fi;
fi;
return L1;
end);

##########################################
# The next function computes local indices for
# quaternion algebras over the rationals.  For
# quaternion algebras over larger number fields,
# we convert to quadratic algebras and use the
# cyclotomic algebra functions.
############################################
InstallGlobalFunction( LocalIndicesOfRationalQuaternionAlgebra, function(A)
local b,D1,D2,p,i,j,M,F,F1,L;

L:=fail;
if LeftActingDomain(A)=Rationals then
D1:=[];
D2:=[];
b:=Elements(Basis(A));
p:=Sum(Coefficients(Basis(A),b[3]^2));
F:=Factors(p);
for i in [1..Size(F)] do
if (p/(F[i]^2) in Integers) then
  p:=p/F[i]^2;
fi;
od;
F:=Factors(p);
for i in [1..Size(F)] do
if F[i]<0 then
  AddSet(D1,-1);
  AddSet(D1,-F[i]);
else
  AddSet(D1,F[i]);
fi;
od;

p:=Sum(Coefficients(Basis(A),b[2]^2));
F:=Factors(p);
for i in [1..Size(F)] do
if (p/(F[i]^2) in Integers) then
  p:=p/F[i]^2;
fi;
od;
F:=Factors(p);
for i in [1..Size(F)] do
if F[i]<0 then
  AddSet(D2,-1);
  AddSet(D2,-F[i]);
else
  AddSet(D2,F[i]);
fi;
od;

L:=[];
for i in [1..Size(D1)] do
for j in [1..Size(D2)] do
  M:=LocalIndicesOfRationalSymbolAlgebra(D1[i],D2[j]);
  L:=LocalIndicesOfTensorProductOfQuadraticAlgs(L,M);
od;
od;
fi;

return L;
end);

##############################################
# The next function checks if a Rational Quaternion Algebra
# is a division ring.
##############################################
InstallGlobalFunction( IsRationalQuaternionAlgebraADivisionRing, function(A)
local L,V;

L:=LocalIndicesOfRationalQuaternionAlgebra(A);
if L=[] then
V:=false;
else
V:=true;
fi;

return V;
end);

#################################################
InstallGlobalFunction( SchurIndex, function(A)
 local m,i,l,L,B,C,D;

m:="fail: Unrecognized Algebra";
if IsAlgebra(A) then
  if IsQuaternionCollection(Basis(A)) then
    if LeftActingDomain(A)=Rationals then
      L:=LocalIndicesOfRationalQuaternionAlgebra(A);
      l:=Length(L);
      m:=1;
      if l>0 then 
        m:=L[1][2]; 
      fi;
      if l>1 then 
        for i in [2..l] do 
          m:=Lcm(m,L[i][2]); 
        od; 
      fi;
    fi;
  else
    m:="fail: Quaternion Algebra Over NonRational Field, use another method.";
  fi;
fi;

if IsRecord(A) then 
  m:=A.SchurIndex; 
fi;

if IsList(A) then
  l:=Length(A);
    if Length(A)=2 and IsField(A[2]) then 
      m:=1; 
    fi;
    if Length(A)=2 and IsRecord(A[2]) then 
      m:=A[2].SchurIndex; 
    fi;
  if Length(A)=3 and IsField(A[1]) and IsField(A[2]) then 
    m:="fail: Cyclic Algebra, use another method.";
  fi;
  if Length(A)=4 then
    L:=LocalIndicesOfCyclicCyclotomicAlgebra(A);
    m:=GlobalSchurIndexFromLocalIndices(L);
  fi;
  if Length(A)=5 then
    L:=LocalIndicesOfCyclotomicAlgebra(A);
 m:=GlobalSchurIndexFromLocalIndices(L);
 fi;
fi;

return m;
end);

############################################
InstallGlobalFunction( SchurIndexByCharacter, function(F,G,n)
local m,A,B,n1;

#A:=SimpleComponentOfGroupRingByCharacter(F,G,n);
B:=Irr(G);
if IsOrdinaryCharacter(n) then
  n1:=Position(Irr(G),n);
else
  n1:=n;
fi;
A:=SimpleComponentByCharacterDescent(F,G,n1);
m:=SchurIndex(A);

return m;
end);
#############################################
InstallGlobalFunction( SimpleComponentByCharacterAsSCAlgebra, function(F,G,n)
local chi,F0,y0,y,F1,I,g,a,A;

if IsPosInt(n) then
  if HasOrdinaryCharacterTable(G) then
    chi:=Irr(G)[n];
  else
    Error("The group has no ordinary character table yet. To avoid randomisation errors, you should compute it first\n");
  fi;
elif IsCharacter(n) then
  chi:=n;
else
  Error("The third argument must be a character or its number\n");
fi;

F0:=Field(chi);
y0:=PrimitiveElement(F0);
y:=PrimitiveElement(F);
F1:=Field([y,y0]);
I:=IrreducibleRepresentationsDixon(G,chi);
g:=Image(I);
a:=Algebra(F1,GeneratorsOfGroup(g));
A:=Image(IsomorphismSCAlgebra(a));

return A;
end);
############################
InstallGlobalFunction( CyclotomicAlgebraAsSCAlgebra, function(A)
local g,m,F,a;

g:=DefiningGroupAndCharacterOfCyclotAlg(A);
#g:=DefiningGroupOfCyclotomicAlgebra(A);
#m:=DefiningCharacterOfCyclotomicAlgebra(A);
F:=A[2];
a:=SimpleComponentByCharacterAsSCAlgebra(F,g[1],g[2]);
#a:=SimpleComponentByCharacterAsSCAlgebra(F,g,m);

return a;
end);
###########################
InstallGlobalFunction( WedderburnDecompositionAsSCAlgebras, function(R)
local W,l,W1,A,i;

W:=WedderburnDecompositionInfo(R);
l:=Size(W);
W1:=[];
for i in [1..l] do
  if Size(W[i])=2 then
  if W[i][1]=1 then
     W1[i]:=W[i][2];
  else
     W1[i]:=MatrixAlgebra(W[i][2],W[i][1]);
  fi;
  fi;
  if Size(W[i])>2 then
  if W[i][1]=1 then
    W1[i]:=CyclotomicAlgebraAsSCAlgebra(W[i]);
  else
    A:=CyclotomicAlgebraAsSCAlgebra(W[i]);
    W1[i]:=MatrixAlgebra(A,W[i][1]);
  fi;
  fi;
od;

return W1;
end);

##########################

################################################
# AntiSymMatUpMat is a technical function which outputs an
# antisymmetric with input matrix from its upper
# triangular part.
# The input is given by a list of list of decreasing length
################################################

InstallGlobalFunction( AntiSymMatUpMat, function(x)
#AntiSymMatUpMat := function(x)
local k,y,i;
k := Length(x)+1;
y:=List([1..k],i->[]);

for i in [1..k-1] do
  y[i] := Concatenation(List([1..i-1], j -> -x[j,i-j]),[0],x[i]);
od;
y[k] := Concatenation(List([1..k-1], j -> -x[j,k-j]),[0]);

return y;
end);

################################################
# KillingCocycle outputs the numerical information
# describing a cyclotomic algebra, equivalent to the input
# which is also the numerical information of a cyclotomic
# algebra, trying to put as zeroes in the fifth entry as possible.
################################################

InstallGlobalFunction( KillingCocycle, function(A)
#KillingCocycle := function(A)

local n,F,m,hrs,k,c,d,e,i,h,r,s,md,x,a,as,j,r1;

if Length(A) < 5 then
  return A;
fi;

n := A[1];
F := A[2];
m := A[3];
hrs := List(A[4],x->List(x,y->y));
k := Length(hrs);
c := List(A[5],x->List(x,y->y));
e := AntiSymMatUpMat(c);

for i in [1..k] do
  h:=hrs[i][1];
  r:=hrs[i][2];
  s:=hrs[i][3];
  if r mod m <> 1 then
    x:=Filtered([1..k],j-> e[i,j] mod Gcd(m,hrs[j][2]-1) = 0);
    if Length(x) > 0 then
      as:=[];
      for j in x do
        r1 := hrs[j][2];
        d := Gcd(m,r1-1);
        md := m/d;
        a := Int(ZmodnZObj(e[i,j]/d,md)*ZmodnZObj((r1-1)/d,md)^-1 );
        Add(as,List([0..d-1],y->a+y*md mod m));
      od;
      as:=Intersection(as);
      if Size(as)>1 then
        a:=as[1];
        hrs[i][3]:=(s+a*(r^h-1)/(r-1)) mod m;
        for j in Difference(x,[i]) do
          if j>i then
            c[i][j-i]:=0;
          else
            c[j][i-j]:=0;
          fi;
        od;
        e := AntiSymMatUpMat(c);
      fi;
    fi;
  fi;
od;

return [n,F,m,hrs,c];

end);

################################################
# IsCyclotomicExtension checks whether the input are two number fields 
# and the first is a cyclotomic extension of the second
################################################

InstallGlobalFunction( IsCyclotomicExtension, function(K,F)
#IsCyclotomicExtension := function(K,F)

local c,pr,i;

if not IsNumberField(K) or not IsNumberField(F) or not IsSubset(K,F) then 
  return fail;
fi;

c:=Conductor(K);
pr := PrimitiveElement(F);

for i in [1..c] do
  if Field([pr,E(i)])=K then 
    return true;
  fi;
od;

return false;

end);

################################################
# ReducingCyclotomicAlgebra TRIES TO REDUCE THE NUMBER OF GENERATORS
################################################

InstallGlobalFunction( ReducingCyclotomicAlgebra, function(A)
#ReducingCyclotomicAlgebra := function(A)

local m,F,pF,K,con,k,e,i,x,y,F1,F2,c,g,h,n,act,coc1,coc2,l,pos,ex,a,j,A1,split,v,w,d1,cls,gcls,lc,A2,s1,s2;

if Length(A) < 5 then
  return fail;
fi;

m:=A[3];
F:=A[2];
pF := PrimitiveElement(F);
K := Field([pF,E(m)]);
con := Lcm(2,Conductor(K));
e := AntiSymMatUpMat(A[5]);
k := Length(A[4]);


for x in [1..k] do
  F2:=NF(con,[A[4][x][2]]);
  if ForAll([1..k],j->e[x][j]=0) and IsCyclotomicExtension(F2,F) then
    y := Difference([1..k],[x]);
    F1:=NF(con,List(y,j->A[4][j][2]));
    c := Lcm(2,Conductor(F1));
    g := E(c);
    while not g in F1 do
      g:=g*E(c);
    od;
    h := E(con)^A[4][x][3];
    n := Order(g);
    split := false;
    a := 1;
    for i in [0..n-1] do
      split := h=Norm(F1,F,a);
      if split then
        break;
      fi;
      a:=a*g;
    od;
## If not we check whether it is a cyclotomic and compute the Schur index to check whether it is 1
    if not split and IsCyclotomicExtension(F1,F) then
      a := 1;
      ex := 0;
      for i in [0..n-1] do
        if a = h then
          break;
        else
          ex:=ex+1;
          a:=a*g;
        fi;
      od;
      A1 := [1,F,c,[A[4][x][1],A[4][x][2] mod c,ex]];
      split := SchurIndex(A1)=1;  
    fi;
    if split then
      c := Lcm(2,Conductor(F2));
      g:=E(c);
      act := [];
      l := Length(y);
      pos := 1;
      if l>1 then 
        coc1 := [];    
      fi;
      for i in y do
        a:=1;
        h:=E(con)^A[4][i][3];
        for ex in [0..c] do
          if a = h then
            break;
          fi;
          a:=a*g;
        od;
        if ex=c then 
          Print("\n The algebra is not a genuine cyclotomic algebra \n");
          return fail;
        fi;          
        Add(act,[A[4][i][1],A[4][i][2] mod c,ex]); 
        if pos < l then
          coc2 := [];
          for j in [pos+1..l] do
            a:=1;
            h:=E(con)^A[5][i][y[j]-i];
            for ex in [0..c] do
              if a = h then
                break;
              fi;
              a:=a*g;
            od;
            if ex=c then 
             Print("\n The algebra is not a genuine cyclotomic algebra \n");
              return fail;
            fi;          
            Add(coc2,ex); 
          od;
          Add(coc1,coc2);
        fi;
        pos := pos+1;
      od;
      if l=1 then 
        return [A[1]*A[4][x][1],F,c,act[1]];
      else 
        return [A[1]*A[4][x][1],F,c,act,coc1];
      fi;
    fi;
  fi;
od;

v := [1..k];
d1 := List(v,i->Filtered(v,j->i=j or e[i,j]<>0));
cls := [];

w:=v;

while w <> [] do
  i:=w[1];
  x:=[];
  y:=[i];
  while x<>y do
    x:=y;
    y:=Union(x,Concatenation(List(x,j->d1[j])));
    w:=Difference(w,y);
  od;
  Add(cls,x);
od;

lc := Length(cls);

if lc = 1 then 
  return fail;
fi;

### WE CALCULATE ALL THE UNIONS OF CLASSES WITH AT LEAST TWO GENERATORS AND AT MOST HALF OF THE NUMBER OF GENERATORS

gcls := Filtered(SSortedList(Arrangements(cls),Union),x->Size(x)>1 and 2*Size(x) <= k);
SortBy(gcls,Size);

for x in gcls do
  y := Difference(v,x);
  F1:=NF(con,List(y,j->A[4][j][2]));
  F2:=NF(con,List(x,j->A[4][j][2]));
  if IsCyclotomicExtension(F1,F) and IsCyclotomicExtension(F2,F) then    

# Construction of the first factor
    c := Lcm(2,Conductor(F1));
    g:=E(c);
    act := [];
    l := Length(x);
    pos := 1;
    if l>1 then 
      coc1 := [];
    fi;
    for i in x do
      h:=E(con)^A[4][i][3];
      a := 1;
      ex := 0;
      for j in [0..c-1] do
        if a = h then
          break;
        else
          ex:=ex+1;
          a:=a*g;
        fi;
      od;
      if ex=c then 
        Print("\n The algebra is not a genuine cyclotomic algebra \n");
        return fail;
      fi;
      Add(act,[A[4][i][1],A[4][i][2] mod c,ex]);                  
      if pos < l then
        coc2 := [];
        for j in [pos+1..l] do
          a:=1;
          h:=E(con)^A[5][i][x[j]-i];
          for ex in [0..c] do
            if a = h then
              break;
            fi;
            a:=a*g;
          od;
          if ex=c then 
            Print("\n The algebra is not a genuine cyclotomic algebra \n");
            return fail;
          fi;          
          Add(coc2,ex); 
        od;
        Add(coc1,coc2);
      fi;
      pos := pos+1;
    od;
    if l=1 then 
      A1 := [A[1],F,c,act[1]];
    else 
      A1 := [A[1],F,c,act,coc1];
    fi;

# Construction of the second factor
    c := Lcm(2,Conductor(F2));
    g:=E(c);
    act := [];
    l := Length(y);
    pos := 1;
    if l>1 then 
      coc1 := [];    
    fi;
    for i in y do
      h:=E(con)^A[4][i][3];
      a := 1;
      ex := 0;
      for j in [0..c-1] do
        if a = h then
          break;
        else
          ex:=ex+1;
          a:=a*g;
        fi;
      od;
      if ex=c then 
        Print("\n The algebra is not a genuine cyclotomic algebra \n");
        return fail;
      fi;
      Add(act,[A[4][i][1],A[4][i][2] mod c,ex]);
      
      if pos < l then
        coc2 := [];
        for j in [pos+1..l] do
          a:=1;
          h:=E(con)^A[5][i][y[j]-i];
          for ex in [0..c] do
            if a = h then
              break;
            fi;
            a:=a*g;
          od;
          if ex=c then 
            Print("\n The algebra is not a genuine cyclotomic algebra \n");
            return fail;
          fi;          
          Add(coc2,ex); 
        od;
        Add(coc1,coc2);
      fi;
      pos := pos+1;
    od;
    
    if l=1 then 
      A2 := [A[1],F,c,act[1]];
    else 
      A2 := [A[1],F,c,act,coc1];
    fi;
#    Print("\n", [A1,A2]);
    s1 := SchurIndex(A1);
    s2 := SchurIndex(A2);
    if s1=1 then
      if s2 =1 then
        return [A[1],F];
      else
        A2[1] := A2[1]*Product(A1[4],x->x[1]);
        return A2;
      fi;
    elif s2=1 then
      A1[1] := A1[1]*Product(A2[4],x->x[1]); 
      return A1;
    fi;
  fi;
od;

return fail;

end);



################################################
# FactoringCycAlg TRIES TO WRITE A CYCLOTOMIC ALGEBRA AS
# A TENSOR PRODUCT OF SMALLER ALGEBRAS
################################################

FactoringCycAlg := function(A)

local m,F,pF,K,con,c,e,k,d1,i,v,w,cls,x,y,lc,gcls,F1,F2,g,act,coc1,coc2,l,pos,ex,a,j,A1,A2,h,s1,s2;



if Length(A) < 5 then
  return A;
fi;

m:=A[3];
F:=A[2];
pF := PrimitiveElement(F);
K := Field([pF,E(m)]);
con := Lcm(2,Conductor(K));

### HERE THE GENERATORS OF THE GALOIS GROUP ARE CLASSIFIED INTO CLASSES SO THAT THE GENERATORS IN ONE CLASS COMMUTE
### WITH THE GENERATORS OF THE OTHER CLASSES

c:=A[5];
e := AntiSymMatUpMat(c);
k := Length(A[4]);
v := [1..k];
d1 := List(v,i->Filtered(v,j->i=j or e[i,j]<>0));
cls := [];

w:=v;

while w <> [] do
  i:=w[1];
  x:=[];
  y:=[i];
  while x<>y do
    x:=y;
    y:=Union(x,Concatenation(List(x,j->d1[j])));
    w:=Difference(w,y);
  od;
  Add(cls,x);
od;

lc := Length(cls);

if lc = 1 then 
  return fail;
fi;

### WE CALCULATE ALL THE NON-EMPTY UNIONS OF CLASSES WITH AT MOST HALF OF THE NUMBER OF GENERATORS

gcls := Filtered(SSortedList(Arrangements(cls),Union),x->Size(x)>1 and 2*Size(x) <= k);
SortBy(gcls,Size);

for x in gcls do
  y := Difference(v,x);
  F1:=NF(con,List(y,j->A[4][j][2]));
  F2:=NF(con,List(x,j->A[4][j][2]));
  if IsCyclotomicExtension(F2,F) and IsCyclotomicExtension(F2,F) then    

# Construction of the first factor
    c := Lcm(2,Conductor(F1));
    g:=E(c);
    act := [];
    l := Length(x);
    pos := 1;
    if l>1 then 
      coc1 := [];
    fi;
    for i in x do
      h:=E(con)^A[4][i][3];
      a := 1;
      ex := 0;
      for j in [0..c-1] do
        if a = h then
          break;
        else
          ex:=ex+1;
          a:=a*g;
        fi;
      od;
      if ex=c then 
        Print("\n The algebra is not a genuine cyclotomic algebra \n");
        return fail;
      fi;
      Add(act,[A[4][i][1],A[4][i][2] mod c,ex]);                  
      if pos < l then
        coc2 := [];
        for j in [pos+1..l] do
          a:=1;
          h:=E(con)^A[5][i][x[j]-i];
          for ex in [0..c] do
            if a = h then
              break;
            fi;
            a:=a*g;
          od;
          if ex=c then 
            Print("\n The algebra is not a genuine cyclotomic algebra \n");
            return fail;
          fi;          
          Add(coc2,ex); 
        od;
        Add(coc1,coc2);
      fi;
      pos := pos+1;
    od;
    if l=1 then 
      A1 := [A[1],F,c,act[1]];
    else 
      A1 := [A[1],F,c,act,coc1];
    fi;

# Construction of the second factor
    c := Lcm(2,Conductor(F2));
    g:=E(c);
    act := [];
    l := Length(y);
    pos := 1;
    if l>1 then 
      coc1 := [];    
    fi;
    for i in y do
      h:=E(con)^A[4][i][3];
      a := 1;
      ex := 0;
      for j in [0..c-1] do
        if a = h then
          break;
        else
          ex:=ex+1;
          a:=a*g;
        fi;
      od;
      if ex=c then 
        Print("\n The algebra is not a genuine cyclotomic algebra \n");
        return fail;
      fi;
      Add(act,[A[4][i][1],A[4][i][2] mod c,ex]);
      
      if pos < l then
        coc2 := [];
        for j in [pos+1..l] do
          a:=1;
          h:=E(con)^A[5][i][y[j]-i];
          for ex in [0..c] do
            if a = h then
              break;
            fi;
            a:=a*g;
          od;
          if ex=c then 
            Print("\n The algebra is not a genuine cyclotomic algebra \n");
            return fail;
          fi;          
          Add(coc2,ex); 
        od;
        Add(coc1,coc2);
      fi;
      pos := pos+1;
    od;
    
    if l=1 then 
      A2 := [A[1],F,c,act[1]];
    else 
      A2 := [A[1],F,c,act,coc1];
    fi;
#    Print("\n", [A1,A2]);
    s1 := SchurIndex(A1);
    s2 := SchurIndex(A2);
    if s1=1 then
      if s2 =1 then
        return [A[1],F];
      else
        A2[1] := A2[1]*Product(A1[4],x->x[1]);
        return A2;
      fi;
    elif s2=1 then
      A1[1] := A1[1]*Product(A2[4],x->x[1]); 
      return A1;
    fi;
  fi;
od;

return fail;

end;
