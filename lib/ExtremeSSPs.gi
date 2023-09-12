###########################################################################
## IsMaximalAbelianFactorGroup:=function(N,A,D)
##
## The function IsMaximalAbelianFactorGroup checks whether A/D
## is a maximal abelian subgroup of N/D.
##
InstallGlobalFunction( IsMaximalAbelianFactorGroup, function(N,A,D) 

local  RTND,  # right transversal of D in N
       RTND0, # elements of RTND which are not in A 
       RTAD,  # right transversal of D in A
       RTAD0, # elements of RTAD which are not in D
       x,     # x in RTND0
       y;     # y in RTAD0

RTND:=RightTransversal(N,D);
RTND0:=Filtered(RTND,x->not x in A);
RTAD:=RightTransversal(A,D);
RTAD0:=Filtered(RTAD,x->not x in D);
for x in RTND0 do
  if ForAll(RTAD0,y->Comm(x,y) in D) then
    return false;
  fi;
od; 
return true;

end);


#############################################################################
## ExtSSPAndDim:=function(G)
##
## The attribute ExtSSPAndDim of G returns a record with components 
## ExtremelyStrongShodaPairs and SumDimension where
## 
## ExtremelyStrongShodaPairs = list of extremely strong Shoda pairs of G 
## that cover the complete set of primitive central idempotents of QG 
## realizable by extremely strong Shoda pairs, 
##
## SumDimension= Sum of Q-dimensions of the simple components corresponding
## to extremely strong Shoda pairs of G.
##
InstallMethod( ExtSSPAndDim,
    "for a finite group",
    true,
  [ IsGroup and IsFinite ],
  0,
function(G)
local ESSP,   # set of representatives of extremely strong Shoda pairs of G
      SumDim, # sum of the Q-dimensions of simple algebras associated elements of ESSP
      ND,     # proper normal subgroups of G in non-increasing order 
      DG,     # derived subgroup of G
      N,      # N in ND, a normal subgroup of G 
      AN,     # a normal subgroup of G, containing N such that AN/N is abelian normal 
              # subgroup of maximal order in G/N 
      ANs,    # list of AN 
      NA,     # Normal subgroups of AN
      NAC,    # normal subgroups D of AN, such that AN/D is cyclic 
      RNAC,   # representatives of G-conjugates in NAC 
      RNACs,  # list of RNAC 
      cont,   # loop controller
      i,j,    # counters  
      D,      # an element of NAC  
      NzD;    # Normalizer of D in G     
        
# INITIALIZATION 

ESSP:=[[G,G]];
SumDim:=1;
if SumDim=Size(G) then
  return rec(ExtremelyStrongShodaPairs:=ESSP, SumDimension:=SumDim);
fi;  
ND:=Difference(NormalSubgroups(G),[G]);
Sort(ND,function(a,b) return Size(a)>=Size(b); end);
DG:=DerivedSubgroup(G);
ANs:=[];
RNACs:=[];
# main loop running on the normal subgroups of G
for N in ND do
  # we first check if (G,N) is an extremely strong Shoda pair.
  if IsSubgroup(N,DG) then
    if IsCyclic(G/N) then
      Add(ESSP,[G,N]);
      SumDim:=SumDim+Phi(Size(G)/Size(N));
      if SumDim=Size(G) then
        SetIsNormallyMonomial( G , true );
        SetIsStronglyMonomial( G , true );
        # if G is normally monomial, then clearly G is strongly monomial
        return rec(ExtremelyStrongShodaPairs:=ESSP, SumDimension:=SumDim);
      fi;  
    fi;
  elif IsCyclic(Centre(G/N)) then 
  # the essps of the form (AN,N) are then collected  
    cont:=true;
    i:=0;
    while cont do
      i:=i+1;
      AN:=ND[i];
      if Size(AN)=Size(N) then
        cont:=false;
      elif IsSubgroup(AN,N) and IsAbelian(AN/N) then
        cont:=false;
        if IsCyclic(AN/N) then 
          if IsMaximalAbelianFactorGroup(G,AN,N) then
            Add(ESSP,[AN,N]); 
            SumDim:=SumDim+Phi(Size(AN)/Size(N))*(Size(G)/Size(AN));
            if SumDim=Size(G) then
              SetIsNormallyMonomial( G , true );
              SetIsStronglyMonomial( G , true );
              return rec(ExtremelyStrongShodaPairs:=ESSP, SumDimension:=SumDim);
            fi;    
          fi;  
        elif not AN in ANs then # no repetition of AN
          Add(ANs,AN); # AN for which AN/N is not cyclic are collected in ANs
          NA:=NormalSubgroups(AN); 
          # for AN in ANs, we now find possible subgroups D, so that (AN,D) is an 
          # essp and Core(D)=N. Moreover, for every AN in ANs, corresponding set of 
          # possible choices of D is saved in RNACs in the same position.  
          NAC:=Filtered(NA,x->IsCyclic(FactorGroupNC(AN,x))); # x is normal in AN
          RNAC:=[];
          while NAC<>[] do
            D:=NAC[1]; 
            NAC:=Difference(NAC,ConjugateSubgroups(G,D)); 
            # conjugate subgroups will yield equivalent essps 
            if not Core(G,D)=N then
              Add(RNAC,D);
            else  
              NzD:=Normalizer(G,D);
              if IsMaximalAbelianFactorGroup(NzD,AN,D) then
                Add(ESSP,[AN,D]); 
                SumDim:=SumDim+((Size(G)/Size(NzD))^2)*((Phi(Size(AN)/Size(D)))*(Size(NzD)/Size(AN)));
                if SumDim=Size(G) then
                  SetIsNormallyMonomial( G , true );
                  SetIsStronglyMonomial( G , true );
                  return rec(ExtremelyStrongShodaPairs:=ESSP, SumDimension:=SumDim);
                fi;  
              fi;
            fi;
          od;
          Add(RNACs,RNAC); 
        else 
          j:=Position(ANs,AN);
          RNAC:=RNACs[j];
          for D in RNAC do
            if Core(G,D)=N then
              Difference(RNACs[j],[D]);
              NzD:=Normalizer(G,D);
              if IsMaximalAbelianFactorGroup(NzD,AN,D) then
                Add(ESSP,[AN,D]); 
                SumDim:=SumDim+((Size(G)/Size(NzD))^2)*((Phi(Size(AN)/Size(D)))*(Size(NzD)/Size(AN)));
                if SumDim=Size(G) then
                  SetIsNormallyMonomial( G , true );
                  SetIsStronglyMonomial( G , true );
                  return rec(ExtremelyStrongShodaPairs:=ESSP, SumDimension:=SumDim);
                fi;  
              fi;
            fi;
          od; 
        fi;   
      fi;     
    od;
  fi;          
od;
SetIsNormallyMonomial( G , false );
return rec(ExtremelyStrongShodaPairs:=ESSP, SumDimension:=SumDim);

end);              


#############################################################################
## ExtremelyStrongShodaPairs:=function(G)
##
## The function ExtremelyStrongShodaPairs computes a non-redundant list of extremely 
## strong Shoda pairs of group G that covers the components of the rational group 
## algebra QG realizable by extremely strong Shoda pairs.
##
InstallGlobalFunction( ExtremelyStrongShodaPairs, function(G)
return ExtSSPAndDim(G).ExtremelyStrongShodaPairs;

end);


#############################################################################
## IsNormallyMonomial(G)
##
## The property IsNormallyMonomial checks whether a group is normally monomial or not.
##
InstallMethod( IsNormallyMonomial,
    "for finite groups",
    true,
  [ IsGroup ],
  0,
function( G )
if IsAbelian(G) then 
  return true;
elif IsAbelian(DerivedSubgroup(G)) then 
  return true;   # a finite metabelian group is normally monomial. 
elif ExtSSPAndDim(G).SumDimension=Size(G) then 
  return true;
else   
  return false;
fi;

end);

InstallTrueMethod( IsStronglyMonomial, IsNormallyMonomial );

#############################################################################
## PrimitiveCentralIdempotentsByExtSSP( QG )
##
## The function PrimitiveCentralIdempotentsByExtSSP computes the set of 
## primitive central idempotents of the group algebra QG, realizable by
## extremely strong Shoda pairs of G.
##
InstallMethod( PrimitiveCentralIdempotentsByExtSSP,
    "for a rational group algebra",
    true,
  [ IsSemisimpleRationalGroupAlgebra ],
  0,
function(QG)
local  G,    # underlying group of QG 
       PCIs, # the list of primitive central idempotents of QG
       ESSPD,# a complete and non-redundant list of extremely strong Shoda pairs of G
             # and Sum of the dimensions of simple components corresponding to them.
       P,    # an extremely strong Shoda pair of G from ESSPD
       IdP;  # the primitive central idempotents of QG associated to P

PCIs:=[];
G:=UnderlyingMagma(QG);
ESSPD:=ExtSSPAndDim(G);
for P in ESSPD.ExtremelyStrongShodaPairs do 
  IdP:=Idempotent_eGsum(QG,P[1],P[2])[2]; #idempotent of QG associated to P
  Add(PCIs,IdP);
od;

return PCIs;

end);




#############################################################################
## PrimitiveCentralIdempotentsByESSP( QG )
##
## The function PrimitiveCentralIdempotentsByESSP computes the set of 
## primitive central idempotents of the group algebra QG, realizable by
## extremely strong Shoda pairs of G.
##
InstallGlobalFunction( PrimitiveCentralIdempotentsByESSP, function(QG)
local  G;
G:=UnderlyingMagma(QG);
if not IsNormallyMonomial(G) then 
  Print("Wedderga: Warning!!!\nThe output is a NON-COMPLETE list of prim. central idemp.s of the input! \n");
fi;
return PrimitiveCentralIdempotentsByExtSSP(QG);
end);


#############################################################################
## SearchingNNKForSSP(QG,H)
##
## The function SearchingNNKForSSP searches a non-normal subgroup K such 
## that (K,H) is a strong Shoda Pair of G and returns [ [ K, H ], e( G, K, H ) ] 
## or returns fail, if such K doesn't exist.
##
InstallGlobalFunction( SearchingNNKForSSP, function(QG,H)
local  G,   # underlying group of QG
       NH,  # Normalizer of H in G
       Epi, # NH --> NH/H      
       NHH, # NH/H
       L,   # <NHH',Z(NHH)>
       Cen, # Centralizer of L in NHH
       K,   # The non-normal subgroup searched
       KH,  # K/H
       X;   # a subset of Cen

G:=UnderlyingMagma(QG);
NH:=Normalizer(G,H);
Epi:=NaturalHomomorphismByNormalSubgroup(NH,H);
NHH:=Image(Epi,NH);
L:=ClosureSubgroup(DerivedSubgroup(NHH),Centre(NHH));
if IsCyclic(L) then 
 Cen:=Centralizer(NHH,L);
  if IsAbelian(Cen) then
    if IsCyclic(Cen) and Centralizer(NHH,Cen)=Cen then
     K:=PreImagesNC(Epi,Cen);
      if IsNormal(G,K)=false then
        return Idempotent_eGsum(QG,K,H);
      fi;
    else
      return fail;
    fi;
  else
   X:=Difference(Elements(Cen),Elements(L));
    while X<>[] do
     KH:=ClosureGroup(L,[X[1]]);
      if IsCyclic(KH) and Centralizer(NHH,KH)=KH then
       K:=PreImagesNC(Epi,KH);
        if IsNormal(G,K)=false then
          return Idempotent_eGsum(QG,K,H);
        fi;
      fi;
       X:=Difference(X,KH);
    od;
  fi;
fi;
return fail;
end);


#############################################################################
## SSPNonESSPAndTheirIdempotents(QG)
##
## The attribute SSPNonESSPAndTheirIdempotents of the rational group algebra QG 
## returns a record of 2 components namely NonExtremelyStrongShodaPairs and 
## PCIsByNonESSPs, where
##
## NonExtremelyStrongShodaPairs = list of non-equivalent strong Shoda pairs of G
## which cover the simple components of QG, not covered by extremely strong Shoda 
## pairs of G.
##
## PCIsByNonESSPs = set of primitive central idempotents of QG realizable by strong
## Shoda pairs in NonExtremelyStrongShodaPairs.
##
InstallMethod( SSPNonESSPAndTheirIdempotents, 
    "for rational group algebra", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra ], 
    0,
function(QG)
local  G,       # underlying group of QG
       ESSPD,   # output of the function ExtSSPAndDim(G)
       SumDim,  # sum of the dimensions of simple algebras covered by strong Shoda pairs
       IdsESSP, # set of idempotents associated to extremely strong Shoda pairs of G
       IdsNESSP,# set of idempotents associated to other strong Shoda pairs of G
       NESSPs,  # set of representatives of strong Shoda pairs of G which are not 
                # extremely strong
       Con,     # conjugacy classes of G;
       Con1,    # conjugacy classes of G of Size greater than 1;
       C,       # C in Con1 
       H,       # representative of C
       HKId,    # strong Shoda pair (H,K) and its idempotent 
       SSP,     # SSP in NESSPs      
       IdSSP,   # the primitive central idempotents of QG associated toSSPP 
       Nz;      # normalizer of K in G 

#Initialization
G:=UnderlyingMagma(QG); 
ESSPD:=ExtSSPAndDim(G); 
SumDim:=ESSPD.SumDimension;

if SumDim=Size(G) then 
  SetIsStronglyMonomial( G , true );
  return rec(NonExtremelyStrongShodaPairs:=[], PCIsByNonESSPs:=[]);
  # if G is normally monomial, then every strong Shoda pair of G is an
  # extremely strong Shoda pair
fi;  

IdsESSP:=PrimitiveCentralIdempotentsByExtSSP(QG);
IdsNESSP:=[];
NESSPs:=[]; 
Con:=ConjugacyClassesSubgroups(G);
Con1:=Filtered(Con,x->Size(x)>1);
#here, we pick a non-abelian subgroup from different conjugacy class
for C in Con1 do
  H:=Representative(C);
  if IsCyclic(Centre(G/Core(G,H))) then
    HKId:=SearchingNNKForSSP(QG,H);
    if not HKId=fail then
      SSP:=HKId[1]; 
      IdSSP:=HKId[2];
      if not IdSSP in Union(IdsESSP,IdsNESSP) then 
        Add(IdsNESSP,IdSSP); 
        Add(NESSPs,SSP);
        Nz:=Normalizer(G,SSP[2]);
        SumDim:=SumDim+((Size(G)/Size(Nz))^2)*((Phi(Size(SSP[1])/Size(SSP[2])))*(Size(Nz)/Size(SSP[1])));
        if SumDim=Size(G) then 
          SetIsStronglyMonomial( G , true );
          return rec(NonExtremelyStrongShodaPairs:=NESSPs, PCIsByNonESSPs:=IdsNESSP);         
        fi;
      fi;  
    fi;
  fi; 
od;  
SetIsStronglyMonomial( G , false );      
return rec(NonExtremelyStrongShodaPairs:=NESSPs, PCIsByNonESSPs:=IdsNESSP); 
      
end);


#############################################################################
##
## PrimitiveCentralIdempotentsByStrongSP(QG)
##
## The function PrimitiveCentralIdempotentsByStrongSP computes the set of primitive
## central idempotents of the group algebra QG, realizable by strong Shoda pairs 
## of G.
##
InstallMethod( PrimitiveCentralIdempotentsByStrongSP,
    "for rational group algebra",
    true, 
    [ IsSemisimpleRationalGroupAlgebra ], 
    0,
    function(QG)
local  IdsESSP, IdsSSP, G; 

IdsESSP:=PrimitiveCentralIdempotentsByExtSSP(QG); 
IdsSSP:=SSPNonESSPAndTheirIdempotents(QG).PCIsByNonESSPs; 

G:=UnderlyingMagma(QG);
if not IsStronglyMonomial(G) then 
  Print("Wedderga: Warning!!!\nThe output is a NON-COMPLETE list of prim. central idemp.s of the input! \n");
fi;
return Concatenation(IdsESSP,IdsSSP);  

end);


#############################################################################
##
## IsExtremelyStrongShodaPair( G, H, K )
##
## The function IsExtremelyStrongShodaPair verifies if (H,K) is an extremely
## strong Shoda pair of G
##
InstallMethod( IsExtremelyStrongShodaPair,
    "for a group and two subgroups", 
    true,
    [ IsGroup, IsGroup, IsGroup ], 
    0,
function( G, H, K )
if IsNormal(G,H) then
  return IsStrongShodaPair( G, H, K );
else
  Info(InfoWedderga, 2, "Wedderga: IsSSP: <H> is not normal in <G>");
  return false;
fi;

end);

#############################################################################