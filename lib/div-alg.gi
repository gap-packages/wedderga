
####################################################################################
# The following functions are needed for "CrossedProductWithDivisionAlgebraPart",  #
# which expresses the division algebra part of a simple component output by        # 
# "WedderburnDecompositionInfo" when the output has 4 terms (cyclic cyclotomic     #
# algebras) or in special cases when the output has 5 terms (a cyclotomic crossed  #
# product of Y type) that occurs in a group of order up to 200.  In the latter     #
# case, the program searches a small library to determine the local indices, so    # 
# this has limited functionality.  For general Y-type crossed products, the        #
# function "DecomposeYAlgebra" can be used to decompose the algebra into the direct# 
# product of two cyclic algebras, whose Schur indices can be computed by solving   #
# appropriate norm equations.                                                      #
####################################################################################

#############################
# Given a Schur algebra output from "wedderga" with 5 terms
# that decomposes as the tensor product of two generalized 
# quaternion algebras, the first function determines this 
# tensor decomposition.   
# #############################

DecomposeYAlgebra:=function(A)
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
  c:=[];
  for j in [1,2] do
  d:=Conductor(B[i][j]);
  F:=CF(d);
  if not(IsCyclotomicField(B[i][j])) then 
    c[j]:=Trace(F,B[i][1],E(d));
  else
    c[j]:=E(d); 
  fi;
  od; 
  B[i][2]:=FieldByGenerators(c);
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
  c:=[];
  for j in [1,2] do
  d:=Conductor(B[i][j]);
  F:=CF(d);
  if not(IsCyclotomicField(B[i][j])) then 
    c[j]:=Trace(F,B[i][1],E(d));
  else
    c[j]:=E(d); 
  fi;
  od; 
  B[i][2]:=FieldByGenerators(c);
  fi;
od; 
  
else
  B[1][1]:=A[2]; 
  B[1][2]:=NF(n,[1,A[4][2][2]]); 
  B[1][3]:=[E(n)^A[4][1][3]];
  B[2][1]:=A[2];
  B[2][2]:=NF(n,[1,A[4][1][2]]);
  B[2][3]:=[E(n)^A[4][2][3]];
fi;
fi;
B1:=B;
return B1; 

end; 
################################################
# Local Schur indices of cyclic algebras (K/F,a) and (K/F,b) 
# will match when a^-1*b is a norm in the extension K/F. 
# Determining this is difficult in general, but it is at 
# least possible to dispense of easy cases by trying to 
# match a few known norms. 
################################################
NormMatch:=function(K,F,c)
local a,n,m,k,d,i,N;

N:=[];
n:=Conductor(K);
if IsOddInt(n) then n:=2*n; fi;
d:=DivisorsInt(n);
a:=PrimitiveElement(K);
m:=Degree(MinimalPolynomial(F,a));
if m=2 then 
for i in d do 
 if Sqrt(i) in K and not(Sqrt(i) in F) then 
  a:=Sqrt(i);
  break;
 else
 if Sqrt(-i) in K and not(Sqrt(-i) in F) then 
  a:=Sqrt(-i);
  break;
 fi; 
 fi;
od; 
fi;

for d in [-n..n] do 
for k in [0..(m-1)] do
  AddSet(N,Norm(K,F,d+a^k));
od;
od;
m:=Size(N);
for i in [1..m] do 
if not(N[i]=0) then
for d in [1..m] do 
  AddSet(N,N[i]^(-1)*N[d]);
  AddSet(N,N[i]*N[d]);
od; 
fi;
od;
if c in N then 
  return true;
else
  return false;
fi;
end;

######################
##############
# The next three functions are used to convert 
# between the wedderga lists for algebras and 
# actual quaternion algebras.
###############################################
ConvertCyclicAlgToCyclicCyclotomicAlg:=function(A)
local F,K,n,i,j,M,k,l,m,B;

F:=A[1];
K:=A[2];
if IsCyclotomicField(K) then 
n:=Conductor(K);
if IsOddInt(n) then n:=2*n; fi;
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

m:=Order(ANFAutomorphism(K,k));
l:="fails";
for i in [0..n] do 
if NormMatch(K,F,E(n)^(-i)*A[3][1]) then l:=i; break; fi; 
od; 
B:=[1,F,n,[m,k,l]];
else
B:="fails";
fi;
return B;
end; 
#####################################################

ConvertCyclicCyclotomicAlgtoQuadraticAlg:=function(A)
local n,B;

if Length(A)=4 and A[4][1]=2 then 
n:=A[3];
B:=[A[2],CF(n),[E(n)^A[4][3]]];
else
B:="fails";
fi;
return [A[1],B];
end;
########################################################
ConvertQuadraticAlgToQuaternionAlg:=function(A)
local d,n,i,B;

n:=Conductor(A[2]);
i:=0;
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
  if not(i=0) then 
     B:=QuaternionAlgebra(A[1],i,A[3][1]);
  else
     B:="fails";
  fi; 
return B; 
end;
#################################################
ConvertQuaternionAlgToCyclicAlg:=function(A)
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
end;

##########################################
# The next function is needed for local indices 
# rational quaternion algebras.  It computes the
# local index of the symbol algebra (p,q) over Q 
# when p and q are -1 or a prime.  Warning: It 
# will not work when p or q are other integers,
# and it does not check this fact. 
##########################################
LocalIndicesOfRationalSymbolAlgebra:=function(a,b) 
local p,q,L;

p:=a;
q:=b;
L:=[];
if p < q then 
  p:=b;
  q:=a;
fi;

if p=-1 then L:=[["infty",2],[2,2]]; fi;
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
end;

##################################################
LocalIndicesOfTensorProductOfQuadraticAlgs:=function(L,M)
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
end;
##########################################
# The next function computes local indices for 
# quaternion algebras over the rationals.  For 
# quaternion algebras over larger number fields, 
# we convert to quadratic algebra and make use 
# of "DivAlgPartOfCyclicAlg", which searches a 
# local index library. 
############################################

LocalIndicesOfRationalQuaternionAlgebra:=function(A)
local b,D1,D2,p,i,j,M,F,F1,L;

L:="fail";
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
end;

##############################################
# The next function checks if a Rational Quaternion Algebra
# is a division ring.
##############################################

IsRationalQuaternionAlgebraADivisionRing:=function(A)
local L,V;

L:=LocalIndicesOfRationalQuaternionAlgebra(A);
if L=[] then 
V:=false;
else
V:=true;
fi;

return V;
end;

##########################################
# Local index library for generalized quaternion 
# and quadratic algebras that occur for small groups, 
# including all of those needed for groups of order 
# up to 200.  
###########################################

DivAlgPartOfCyclicAlg:=function(A)
local L,D,d;

D:=A;
L:=[1,D];
d:=Trace(A[2],A[1],1);
if NormMatch(A[2],A[1],A[3][1]) then D:=A[1]; L:=[d,D]; fi; 

if A[1]=Rationals and A[2]=CF(4) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then 
   D:=["DivAlg",Rationals,"LocInds=",[["infty",2],[2,2]]]; L:=[1,D]; fi;

if A[1]=Rationals and A[2]=CF(3) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then  
   D:=["DivAlg",Rationals,"LocInds=",[["infty",2],[3,2]]]; L:=[1,D]; fi;

if A[1]=Rationals and A[2]=CF(3) and NormMatch(A[2],A[1],(-2)^(-1)*A[3][1]) then 
   D:=["DivAlg",Rationals,"LocInds=",[["infty",2],[2,2]]]; L:=[1,D]; fi;

if A[1]=Rationals and A[2]=CF(3) and NormMatch(A[2],A[1],(2)^(-1)*A[3][1]) then 
   D:=["DivAlg",Rationals,"LocInds=",[[2,2],[3,2]]]; L:=[1,D]; fi;

if A[1]=Rationals and A[2]=NF(12,[1,11]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",Rationals,"LocInds=",[[2,2],[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(5,[1,4]) and A[2]=CF(5) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(5,[1,4]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(5,[1,4]) and A[2]=CF(5) and NormMatch(A[2],A[1],(2)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(5,[1,4]),"LocInds=",[[2,2],[5,2]]]; L:=[1,D]; fi;

if A[1]=NF(5,[1,4]) and A[2]=CF(5) and NormMatch(A[2],A[1],(-2)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(5,[1,4]),"LocInds=",[["infty",2],[2,2],[5,2]]]; L:=[1,D]; fi;

if A[1]=NF(5,[1,4]) and A[2]=NF(20,[1,9]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(5,[1,4]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(5,[1,4]) and A[2]=NF(15,[1,4]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(5,[1,4]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(8,[1,7]) and A[2]=NF(24,[1,7]) and NormMatch(A[2],A[1],(2-Sqrt(2))^(-1)*A[3][1]) then D:=["DivAlg",NF(8,[1,7]),"LocInds=",[[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(8,[1,7]) and A[2]=NF(24,[1,7]) and NormMatch(A[2],A[1],(2+Sqrt(2))^(-1)*A[3][1]) then D:=["DivAlg",NF(8,[1,7]),"LocInds=",[[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(8,[1,7]) and A[2]=NF(24,[1,7]) and NormMatch(A[2],A[1],(-2+Sqrt(2))^(-1)*A[3][1]) then D:=["DivAlg",NF(8,[1,7]),"LocInds=",[["infty",2],[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(8,[1,7]) and A[2]=NF(24,[1,7]) and NormMatch(A[2],A[1],(-2-Sqrt(2))^(-1)*A[3][1]) then D:=["DivAlg",NF(8,[1,7]),"LocInds=",[["infty",2],[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(8,[1,7]) and A[2]=NF(24,[1,7]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(8,[1,7]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(8,[1,7]) and A[2]=CF(8) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then 
   D:=["DivAlg",NF(8,[1,7]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(8,[1,3]) and A[2]=NF(24,[1,11]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=NF(8,[1,3]); L:=[2,D]; fi;

if A[1]=NF(24,[1,5,7,11]) and A[2]=NF(24,[1,5]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=NF(24,[1,5,7,11]); L:=[2,D]; fi;

if A[1]=NF(24,[1,5,19,23]) and A[2]=NF(24,[1,19]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(24,[1,5,19,23]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(24,[1,5,19,23]) and A[2]=NF(24,[1,5]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(24,[1,5,19,23]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(7,[1,6]) and A[2]=CF(7) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(7,[1,6]),"LocInds=",[["infty",2],[7,2]]]; L:=[1,D]; fi;

if A[1]=NF(7,[1,6]) and A[2]=CF(7) and NormMatch(A[2],A[1],(2)^(-1)*A[3][1]) then
   D:=NF(7,[1,6]); L:=[2,D]; fi;

if A[1]=NF(7,[1,6]) and A[2]=CF(7) and NormMatch(A[2],A[1],(-2)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(7,[1,6]),"LocInds=",[["infty",2],[7,2]]]; L:=[1,D]; fi;

if A[1]=NF(7,[1,6]) and A[2]=NF(28,[1,13]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(7,[1,6]),"LocInds=",[["infty",2],[2,2]]]; L:=[1,D]; fi;

if A[1]=NF(7,[1,6]) and A[2]=NF(21,[1,13]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(7,[1,6]),"LocInds=",[["infty",2],[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(9,[1,8]) and A[2]=CF(9) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(9,[1,8]),"LocInds=",[["infty",2],[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(9,[1,8]) and A[2]=CF(9) and NormMatch(A[2],A[1],(-2)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(9,[1,8]),"LocInds=",[["infty",2],[2,2]]]; L:=[1,D]; fi;

if A[1]=NF(9,[1,8]) and A[2]=CF(9) and NormMatch(A[2],A[1],(2)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(9,[1,8]),"LocInds=",[[2,2],[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(9,[1,8]) and A[2]=NF(36,[1,17]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then 
   D:=["DivAlg",NF(9,[1,8]),"LocInds=",[["infty",2],[2,2]]]; L:=[1,D]; fi;

if A[1]=NF(40,[1,9,31,39]) and A[2]=NF(40,[1,31]) and NormMatch(A[2],A[1],(2-E(8)+E(8)^3)^(-1)*A[3][1]) then D:=["DivAlg",NF(40,[1,9,31,39]),"LocInds=",[[5,2]]]; L:=[1,D]; fi;

if A[1]=NF(40,[1,9,31,39]) and A[2]=NF(40,[1,31]) and NormMatch(A[2],A[1],(2+E(8)-E(8)^3)^(-1)*A[3][1]) then D:=["DivAlg",NF(40,[1,9,31,39]),"LocInds=",[[5,2]]]; L:=[1,D]; fi;

if A[1]=NF(40,[1,9,31,39]) and A[2]=NF(40,[1,31]) and NormMatch(A[2],A[1],(-2+E(8)-E(8)^3)^(-1)*A[3][1]) then 
   D:=["DivAlg",NF(40,[1,9,31,39]),"LocInds=",[[5,2],["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(40,[1,9,31,39]) and A[2]=NF(40,[1,31]) and NormMatch(A[2],A[1],(-2-E(8)+E(8)^3)^(-1)*A[3][1]) then  
   D:=["DivAlg",NF(40,[1,9,31,39]),"LocInds=",[[5,2],["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(40,[1,9,31,39]) and A[2]=NF(40,[1,31]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(40,[1,9,31,39]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(40,[1,11,29,39]) and A[2]=NF(40,[1,11]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(40,[1,11,29,39]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(40,[1,11,29,39]) and A[2]=NF(40,[1,29]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(40,[1,11,29,39]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;


if A[1]=Rationals and A[2]=NF(20,[1,3,7,9]) and NormMatch(A[2],A[1],(2)^(-1)*A[3][1]) then
   D:=["DivAlg",Rationals,"LocInds=",[[2,2],[5,2]]]; L:=[1,D]; fi;

if A[1]=Rationals and A[2]=NF(20,[1,3,7,9]) and NormMatch(A[2],A[1],(-2)^(-1)*A[3][1]) then
   D:=["DivAlg",Rationals,"LocInds=",[[5,2],["infty",2]]]; L:=[1,D]; fi;

if A[1]=Rationals and A[2]=NF(20,[1,3,7,9]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",Rationals,"LocInds=",[["infty",2],[2,2]]]; L:=[1,D]; fi;

if A[1]=Rationals and A[2]=CF(5) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then 
   D:=["DivAlg",Rationals,"LocInds=",[[5,2],["infty",2]]]; L:=[2,D]; fi;

if A[1]=NF(40,[1,9,11,19]) and A[2]=NF(40,[1,9]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=NF(40,[1,9,11,19]); L:=[2,D]; fi;

if A[1]=NF(40,[1,9,11,19]) and A[2]=NF(40,[1,11]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=NF(40,[1,9,11,19]); L:=[2,D]; fi;

if A[1]=NF(40,[1,19,29,31]) and A[2]=NF(40,[1,19]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=NF(40,[1,19,29,31]); L:=[2,D]; fi;

if A[1]=NF(40,[1,19,29,31]) and A[2]=NF(40,[1,29]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=NF(40,[1,19,29,31]); L:=[2,D]; fi;

if A[1]=NF(11,[1,10]) and A[2]=CF(11) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(11,[1,10]),"LocInds=",[["infty",2],[11,2]]]; L:=[1,D]; fi;

if A[1]=NF(11,[1,10]) and A[2]=NF(44,[1,21]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(11,[1,10]),"LocInds=",[["infty",2],[2,2]]]; L:=[1,D]; fi;

if A[1]=NF(11,[1,10]) and A[2]=CF(11) and NormMatch(A[2],A[1],(-2)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(11,[1,10]),"LocInds=",[[2,2],["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(11,[1,10]) and A[2]=CF(11) and NormMatch(A[2],A[1],(2)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(11,[1,10]),"LocInds=",[[2,2],[11,2]]]; L:=[1,D]; fi;

if A[1]=NF(16,[1,15]) and A[2]=CF(16) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(16,[1,15]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(16,[1,15]) and A[2]=NF(48,[1,31]) and NormMatch(A[2],A[1],(2+E(16)^3-E(16)^5)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(16,[1,15]),"LocInds=",[[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(16,[1,15]) and A[2]=NF(48,[1,31]) and NormMatch(A[2],A[1],(2-E(16)^3+E(16)^5)^(-1)*A[3][1]) then 
   D:=["DivAlg",NF(16,[1,15]),"LocInds=",[[3,2]]]; L:=[1,D]; fi;

if A[1]=NF(16,[1,15]) and A[2]=NF(48,[1,31]) and NormMatch(A[2],A[1],(-2+E(16)^3-E(16)^5)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(16,[1,15]),"LocInds=",[[3,2],["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(16,[1,15]) and A[2]=NF(48,[1,31]) and NormMatch(A[2],A[1],(-2+E(16)-E(16)^7)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(16,[1,15]),"LocInds=",[[3,2],["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(16,[1,15]) and A[2]=NF(48,[1,31]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then
   D:=["DivAlg",NF(16,[1,15]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(48,[1,23,31,41]) and A[2]=NF(48,[1,31]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=NF(48,[1,23,31,41]); L:=[2,D]; fi;

if A[1]=NF(48,[1,23,31,41]) and A[2]=NF(48,[1,23]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=NF(48,[1,23,31,41]); L:=[2,D]; fi;

if A[1]=NF(48,[1,7,41,47]) and A[2]=NF(48,[1,41]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(48,[1,7,41,47]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(48,[1,7,41,47]) and A[2]=NF(48,[1,7]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then D:=["DivAlg",NF(48,[1,7,41,47]),"LocInds=",[["infty",2]]]; L:=[1,D]; fi;

if A[1]=NF(16,[1,7]) and A[2]=CF(16) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then 
   D:=NF(16,[1,7]); L:=[2,D]; fi;

if A[1]=NF(16,[1,7]) and A[2]=NF(48,[1,7]) and NormMatch(A[2],A[1],(-1)^(-1)*A[3][1]) then  
   D:=NF(16,[1,7]); L:=[2,D]; fi;

return L;
end;

############################################
# This Function returns the local indices of 
# the central simple algebra that is the tensor   
# product of two simple algebras presented as 
# matrix algebras over division algebras by 
# wedderga's div-alg functions.
############################################

DivAlgPartOfTensorProductOfCyclicAlgs:=function(A,B)
local S,p,i,j,l,m,L,D,L1,r1,r2,r,m1,m2;

r1:=A[1]; 
r2:=B[1];

if IsField(A[2]) then 
if IsField(B[2]) then 
if A[2]=B[2] then r:=r1*r2; D:=[r,A[2]]; else D:="TensorProductUndefined"; fi;
else 
if A[2]=B[2][2] then r:=r1*r2; D:=[r,B[2]]; else D:="TensorProductUndefined"; fi;
fi;
fi;

if not(IsField(A[2])) and IsField(B[2]) then 
if A[2][2]=B[2] then r:=r1*r2; D:=[r,A[2]]; else D:="TensorProductUndefined"; fi;  
fi;

if not(IsField(A[2])) and not(IsField(B[2])) then 
r:=r1*r2;
S:=[];
L:=[];
L1:=[];
m1:=1;
m2:=1;
for p in [1..Length(A[2][4])] do 
AddSet(S,A[2][4][p][1]);
m1:=Lcm(m1,A[2][4][p][2]);
od;
for p in [1..Length(B[2][4])] do 
AddSet(S,B[2][4][p][1]);
m2:=Lcm(m2,A[2][4][p][2]);
od; 
for i in [1..Size(S)] do 
L[i]:=[];
L[i][1]:=S[i];
m:=0;
for p in [1..Length(A[2][4])] do 
if A[2][4][p][1]=S[i] then 
  m:=A[2][4][p][2]; 
fi; 
od; 
for p in [1..Length(B[2][4])] do 
if B[2][4][p][1]=S[i] then
if not(B[2][4][p][2]=m) then 
m:=Lcm(B[2][4][p][2],m);
else
if (m=2 and B[2][4][p][2]=2) then 
m:=1;
else
m:=["divides",m];
fi; 
fi;
fi;
od; 
L[i][2]:=m;
od;

l:=Length(L);
if l>0 then 
for i in [1..l] do 
if L[i][2]>1 then 
   Add(L1,L[i]);
fi;
od; 
fi;

l:=Length(L1);
m:=1;
if l>0 then
for j in [1..l] do
if (L1[j][2] in Integers) then  
m:=Lcm(m,L1[j][2]);
fi;
od;
fi;
r:=r*m1*m2/m;
if L1=[] then 
  D:=[r,A[2][2]]; 
else 
  D:=[r,["DivAlg",A[2][2],"LocInds=",L1]];
fi;

fi;

return D; 
end;

#############################################
# Necessary arithmetic functions for subroutines
#############################################
PPartOfN:=function(n,p)
local i,a,b;

b:=n;
a:=0;
while b/p in PositiveIntegers do 
b:=b/p;
a:=a+1;
od;

return p^a;
end;

###############################
PDashPartOfN:=function(n,p)
local m;

m:=n; 
while m/p in PositiveIntegers do 
m:=m/p;
od; 

return m; 
end;
######################################
# Cyclotomic reciprocity functions, not needed anymore.
######################################
#ResidueDegreeAtQ:=function(n,q)
#local f,m;
#m:=PDashPartOfN(n,q); 
#f:=1; 
#while not((q^f-1)/m in Integers) do 
#f:=f+1;
#od;
#return f; 
#end; 
################################
#SplittingDegreeAtQ:=function(n,q)
#local m,f,g;
#f:=ResidueDegreeAtQ(n,q);
#m:=PDashPartOfN(n,q);
#g:=Phi(m)/f; 
#return g; 
#end;
################################

################################
# Given a simple component of a rational group algebra whose 
# "WedderburnDecompositionInfo" output has 4 terms, the next 
# three functions compute its indices at odd primes, infinity, 
# and 2.    
################################

LocalIndexAtOddP:=function(A,q)
local m,n,b,k,f,g,e,e1,k1,c,F,K,m1,n1;

m:=1; 
b:=A[4][3];
if b>0 then 
n:=Lcm(Conductor(A[2]),A[3]);
c:=Trace(CF(n),A[2],E(n));
F:=FieldByGenerators([c,E(q-1)]);
n1:=Conductor([E(n),E(q-1)]);
K:=CF(n1);
  if A[3]/q in PositiveIntegers then 
  e1:=OrderMod(A[4][2],q);
   if e1>1 then 
##################
#  Replaced cyclotomic reciprocity functions
#    k:=A[3]/Gcd(A[3],b);
#   f:=ResidueDegreeAtQ(n,q);
#    g:=SplittingDegreeAtQ(n,q);
#    e:=Phi(n)/f*g;
#####################
m1:=PDashPartOfN(n,q); 
f:=1; 
while not((q^f-1)/m1 in Integers) do 
f:=f+1;
od;
####################
    k:=Trace(K,F,1);
    k1:=(q^f-1)/k;
    while not(k1/(Order(E(A[3])^(b*m))) in PositiveIntegers) do 
      m:=m+1;
    od; 
    fi; 
  fi;
fi;

return m;
end;
###############################
# For the computation of the index of a cyclic 
# cyclotomic algebra at infinity, we cannot 
# use the Frobenius-Schur indicator because we 
# have no way of knowing which group character 
# is involved.  Instead we determine the 
# nature of the algebra as a quadratic algebra
# over the reals. 
###############################
LocalIndexAtInfty:=function(A)
local m,n,s,n1;

m:=1;
n:=Lcm(Conductor(A[2]),A[3]);
if n>2 then 
  s:=ANFAutomorphism(A[2],-1);  
  if s=ANFAutomorphism(A[2],1) then 
    n1:=PPartOfN(A[3],2);
    if E(n1)^A[4][3]=-1 then 
      m:=2;
    fi;
  fi;
fi;

return m;
end;
###############################
# For the local index at 2, we detect if the cyclic 
# cyclotomic algebra will be of nonsplit quaternion type  
# over the 2-adics
###############################
LocalIndexAtTwo:=function(A)
local n,m,n1,n2,b,f; 

n:=Lcm(Conductor(A[2]),A[3]);
m:=1;
n1:=PDashPartOfN(n,2); 
#f:=ResidueDegreeAtQ(n1,2);
f:=1; 
while not((2^f-1)/n1 in Integers) do 
f:=f+1;
od;

if IsOddInt(f) then 
n2:=PPartOfN(A[3],2);
b:=PPartOfN(A[4][1],2);
if Phi(n2)/b=1 then
if E(n2)^A[4][3]=-1 then 
m:=2;
fi;
fi;
fi;

return m;
end;

##############################
# Given a group G and a simple component A whose 
# WedderburnDecompositionInfo in wedderga has length 4,
# this program gives the list of local indices at 
# all primes relevant to the rational Schur index
###############################
LocalIndicesOfCyclicCyclotomicAlgebra:=function(A)
local n,S,s,i,L,l,q,L1;

L:=[];
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
L[s+1][1]:="infty";
L[s+1][2]:=l;

L1:=[];

s:=Size(L);
for i in [1..s] do 
if L[i][2]>1 then 
Add(L1,L[i]);
fi;
od;

return L1; 
end;
###########################################
# Program that obtains rational Schur index 
# of a cyclic cyclotomic algebra from its local 
# indices.  Algebras given by a list of of length 3 
# have to be in the form of a quadratic algebra - 
# this can be used to calculate Schur indices of 
# quaternion algebras after they are converted. 
# Algebras given by a list of length 
# 4 or 5 need to be in the form of an algebra output 
# by "CrossedProductWithDivisionAlgebraPart". 
############################################

RationalSchurIndex:=function(A)
local m,i,l,L,B,C,D;

if IsAlgebra(A) and IsQuaternionCollection(Basis(A)) then 
if LeftActingDomain(A)=Rationals then 
L:=LocalIndicesOfRationalQuaternionAlgebra(A);
l:=Length(L);
m:=1;
if l>0 then m:=L[1][2]; fi;
if l>1 then for i in [2..l] do m:=Lcm(m,L[i][2]); od; fi;
else
B:=ConvertQuaternionAlgToCyclicAlg(A);
D:=DivAlgPartOfCyclicAlg(B);
L:=D[2][4];
l:=Length(L);
m:=1;
if l>0 then m:=L[1][2]; fi;
if l>1 then for i in [2..l] do m:=Lcm(m,L[i][2]); od; fi;
fi;
fi;

if IsList(A) then 
if Length(A)<3 then 
L:=[];
m:=1; 
fi; 

if Length(A)=3 then 
B:=DivAlgPartOfCyclicAlg(A);
if Length(B[2])=4 then L:=B[2][4]; fi;
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

if Length(A)=4 then 
L:=LocalIndicesOfCyclicCyclotomicAlgebra(A);
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

if Length(A)=5 then 
  D:=DecomposeYAlgebra(A); 
  C:=[DivAlgPartOfCyclicAlg(D[1]),DivAlgPartOfCyclicAlg(D[2])];
  B:=DivAlgPartOfTensorProductOfCyclicAlgs(C[1],C[2]);
  L:=B[2][4];
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
fi;

return m; 
end;

############################### 
# This function will convert cyclic
# cyclotomic algebras output by wedderga into matrix algebras 
# over Schur division algebras with indicated center and specified  
# local indices at rational primes. 
###############################

DivAlgPartOfCyclicCyclotomicAlgebra:=function(A)
local L,m,r,D,B;

L:=LocalIndicesOfCyclicCyclotomicAlgebra(A);
m:=RationalSchurIndex(A);
r:=A[1]*A[4][1]/m;
if m=1 then 
B:=[r,A[2]];
fi;
if m>1 then
D:=[];
D[1]:="DivAlg";
D[2]:=A[2];
D[3]:="LocInds=";
D[4]:=LocalIndicesOfCyclicCyclotomicAlgebra(A);
B:=[r,D];
fi;

return B;
end; 

#####################################################
# Given a simple algebra output by WedderburnDecompositionInfo
# or SimpleAlgebraByCharacterInfo from wedderga, 
# this program determines its actual matrix degree and division 
# algebra part in terms of local indices at all primes.
# If it cannot find the division algebra part, it leaves the 
# simple algebra as it was.
#####################################################

CrossedProductWithDivisionAlgebraPart:=function(A)
local l,B,D,C;

l:=Length(A); 
B:=A;  
if l=4 then B:=DivAlgPartOfCyclicCyclotomicAlgebra(A); fi;
if l=5 then 
  D:=DecomposeYAlgebra(A); 
  C:=[DivAlgPartOfCyclicAlg(D[1]),DivAlgPartOfCyclicAlg(D[2])];
  B:=DivAlgPartOfTensorProductOfCyclicAlgs(C[1],C[2]);
fi;

return B;
end;

#####################################
# These functions support the main function 
# "DivisionAlgebraPartByBrauerCharacters" which 
# is able to accurately compute the local Schur index 
# of an irreducible character of a finite group G 
# at infty (using the Frobenius-Schur indicator)
# and at any prime p for which the defect group at p 
# of the character is cyclic using Benard's theorem.
# WARNING: Knowledge of Brauer characters in GAP is 
# limited, and only available for certain small 
# p-solvable groups.  Since GAP does not provide the 
# defect group, only its order, we can be certain the 
# answer given by this program is correct only when 
# all p-subgroups of that order are cyclic.  
# In cases where the defect group is not cyclic, we know 
# the answers produced by this function are often incorrect.
######################################

FinFieldExt:=function(G,p,n,n1)
local T,chi,V,Y,h,L,i,z,l,m,K,B,d,M,C,D,b,j,F1,T1,psi,U,k,F2,t;

T:=CharacterTable(G);
chi:=Irr(G)[n];
V:=ValuesOfClassFunction(chi);
Y:=OrdersClassRepresentatives(T);
h:=Size(Y);
L:=[];
for i in [1..h] do if Gcd(Y[i],p) = 1 then Add(L,V[i]); fi; od; 
l:=Size(L);
m:=Conductor(L);
K:=CF(m);
B:=Basis(K);
for i in [1..m] do if (p^i-1)/m in Integers then d:=i; break; fi; od; 
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
#Print(M," ");
F1:=FieldByGenerators(M);

T1:=T mod p; 
psi:=IBr(G,p)[n1];
U:=ValuesOfClassFunction(psi);
m:=Lcm(m,Conductor(U));
K:=CF(m);
B:=Basis(K);
for i in [1..m] do if (p^i-1)/m in Integers then d:=i; break; fi; od; 
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

F2:=FieldByGenerators(M);
t:=LogMod(0,Size(F1),Size(F2));


return t; 
end;
##############################################
PossibleDefectGroups:=function(G,n,p)
local S,U,Q,Q1,i,j,I,T,b,k,d,H,a;

T:=CharacterTable(G);
S:=T mod p; 
b:=BlocksInfo(S);
for j in [1..Size(b)] do 
if n in b[j].ordchars then 
  k:=b[j].modchars[1];
  d:=b[j].defect;
  break;
fi;
od;
Q:=SylowSubgroup(G,p);
U:=[];
if Size(Q)>p^d then 
H:=ConjugacyClasses(G);
for j in [2..Size(H)] do 
  if Gcd(OrdersClassRepresentatives(T)[j],p)=1 then 
     a:=Elements(H[j])[1];
     Q1:=Intersection(Q,Q^a);
     if Size(Q1)=p^d then 
        AddSet(U,ConjugacyClassSubgroups(G,Q1)); 
     fi; 
  fi;
od;
else
AddSet(U,ConjugacyClassSubgroups(G,Q));
fi;

return U;
end;

##############################################
# Determines the division algebra part of the simple component 
# of QG corresponding to the complex character 
# I:=Irr(G)[n] in terms of its local indices at all primes
# using the dimension of the field of values of a Brauer character in the same 
# block over the field of values of I.
# Warning: This works as long as the character lies in a block 
# with cyclic defect group.  
###############################################

DivisionAlgebraPartByBrauerCharacters:=function(F,G,n)
local W,T,I,V,C,L,P,Q,Q1,q,p,f,i,j,l,d,S,a,b,m,t,k,u,H,U,K,D,B,L1;

T:=CharacterTable(G);
I:=Irr(G)[n];
V:=ValuesOfClassFunction(I);
C:=FieldByGenerators(F,V);
L:=[];
W:=[];
f:=Indicator(T,2)[n]; 
if f=-1 then 
L[1]:=["infty",2];
fi;

P:=AsSet(Factors(Size(G)));
q:=Size(P);
if q=1 then 
B:="failsforPgroups";
else

for i in [1..q] do 
p:=P[i];
S:=T mod p; 
b:=BlocksInfo(S);
for j in [1..Size(b)] do 
if n in b[j].ordchars then 
  k:=b[j].modchars[1];
  break;
fi;
od;

U:=PossibleDefectGroups(G,n,p);
f:=0;
for u in [1..Size(U)] do 
 if not(IsCyclic(Elements(U[u])[1])) then 
   f:=f+1;
 fi;
od; 
if not(f=0) then
  if f<Size(U) then  
   AddSet(W,[p,"maybe"]);
  else 
   AddSet(W,[p,"fail"]);
  fi;
fi;
   
t:=FinFieldExt(G,p,n,k);
if t>1 then 
Add(L,[p,t]);
fi;
od; 

l:=Size(L); 
if l=0 then 
B:=[I[1],C];
else
m:=L[1][2];
if l>1 then 
for i in [2..l] do 
  m:=Lcm(m,L[i][2]);
od; 
fi;
D:=["DivAlg",C,"LocInds=",L];
B:=[V[1]/m,D,W];
fi;

fi;

return B;
end;

###################################################
# The next function is useful for converting crossed 
# algebras presented as a record by 
# "WedderburnDecompositionInfo" into actual algebras 
# defined by Structure Constants as in Section 62 of 
# the GAP manual.   
# Waiting for new GAP functions in "Algebras with 
# Finite Presentation" to be developed.
###################################################
#ConvertWDInfoToAlgebraBySC:=function(A)
#local i,j,m,B,M; 
#if Length(A)=1 then 
#M:=A[1];
#fi;
#if Length(A)=2 then 
#M:=MatrixAlgebra(A[2],A[1]);
#fi; 
#if Length(A)=4 then 
#if A[4][3]=0 or RationalSchurIndex(A)=1 then
#M:=MatrixAlgebra(A[2],A[1]*A[4][1]);
#else
#M:=A;
#fi;
#fi;
#if Length(A)=5 then
#if RationalSchurIndex(A)=1 then 
#M:=MatrixAlgebra(A[2],A[1]*A[4][1][1]*A[4][2][1]);
#else
#M:=A;
#fi; 
#fi;
#return M; 
#end;
####################################################

