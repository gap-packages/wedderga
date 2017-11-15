###########################################################################
## IsCyclicMaximalAbelianFactorGroup:=function(N,A,D)
##
## The function IsCyclicMaximalAbelianFactorGroup checks whether A/D is cyclic 
## maximal abelian subgroup of N/D.
##
IsCyclicMaximalAbelianFactorGroup:= function(N,A,D) 

local  ND,   # N/D       
	   AD,   # A/D
	   epi,  # N--->N/D
	   a,    # Generator of AD
       x;    # x in ND

epi := NaturalHomomorphismByNormalSubgroup(N,D);
ND := Image(epi);
AD := Image(epi,A);
if not IsCyclic(AD) then
  return false;
fi;
a := MinimalGeneratingSet(AD)[1];
for x in Difference(ND,AD) do
  if Comm(x,a) = One(ND) then
    return false;
  fi;
od; 
return true;
end;


#############################################################################
## ExtSSPAndDim:=function(G)
##
## The attribute ExtSSPAndDim of the G returns a record with components 
## ExtremelyStrongShodaPairs and SumDimension where
## 
##  ExtremelyStrongShodaPairs = list of extremely strong Shoda pairs of G 
##  that cover the complete set of primitive central idempotents of QG 
##  realizable by extremely strong Shoda pairs, 
##
##  SumDimension = Sum of the Q-SumDimensions of simple components corresponding
##  to extremely strong Shoda pairs of G.
##
ExtSSPAndDim:=function(G)
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
      RNACS,  # list of RNAC 
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
RANCs:=[];
# main loop running on the normal subgroups of G
for N in ND do
  # we first check if (G,N) is an extremely strong Shoda pair.
  if IsSubgroup(N,DG) then
    if IsCyclic(G/N) then
      Add(ESSP,[G,N]);
      SumDim:=SumDim+Phi(Size(G)/Size(N));
      if SumDim=Size(G) then
        return rec(ExtremelyStrongShodaPairs:=ESSP, SumDimension:=SumDim);
      fi;  
    fi;
  elif IsCyclic(Centre(G)/N) then 
    cont:=true;
    i:=0;
    while cont do
      i:=i+1;
      AN:=ND[i];
      if Size(AN)=Size(N) then
        cont:=false;
      elif IsSubgroup(AN,N) and IsAbelian(AN/N) then
        cont:=false;
        if not AN in ANs then
          Add(ANs,AN);
          NA:=NormalSubgroups(AN);
          NAC:=Filtered(NA,x->IsCyclic(AN/x));
          RNAC:=[];
          while NAC<>[] do
            D:=NAC[1]; 
            NAC:=Difference(NAC,ConjugateSubgroups(G,D));  
            if not Core(G,D)=N then
              Add(RNAC,D);
            else  
              NzD:=Normalizer(G,D);
              if IsCyclicMaximalAbelianFactorGroup(NzD,AN,D) then
                Add(ESSP,[AN,D]); 
                SumDim:=SumDim+((Size(G)/Size(NzD))^2)*((Phi(Size(AN)/Size(D)))*(Size(NzD)/Size(AN)));
                if SumDim=Size(G) then
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
              Remove(RNACs[j],D);
              NzD:=Normalizer(G,D);
              if IsCyclicMaximalAbelianFactorGroup(NzD,AN,D) then
                Add(ESSP,[AN,D]); 
                SumDim:=SumDim+((Size(G)/Size(NzD))^2)*((Phi(Size(AN)/Size(D)))*(Size(NzD)/Size(AN)));
                if SumDim=Size(G) then
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

return rec(ExtremelyStrongShodaPairs:=ESSP, SumDimension:=SumDim);

end;              


#############################################################################
## ExtStrongShodaPairs:=function(G)
##
## The function ExtStrongShodaPairs computes a list of extremely strong Shoda pairs
## of the group G that covers the complete set of primitive central idempotents of 
## the rational group algebra QG realizable by extremely strong Shoda pairs.
##
ExtStrongShodaPairs:=function(G)
return ExtSSPAndDim(G).ExtremelyStrongShodaPairs;
end;


#############################################################################
## IsNormallyMonomial(G)
##
## The function IsNormallyMonomial checks whether a group is normally monomial or not.
##
IsNormallyMonomial:=function(G)
if IsAbelian(G) then 
  return true;
elif IsAbelian(DerivedSubgroup(G)) then 
  return true;   # a finite metabelian group is normally monomial. 
elif ExtSSPAndDim(G).SumDimension=Size(G) then 
  return true;
else   
  return false;
fi;
end;


#############################################################################
## PrimitiveCentralIdempotentsByExtSSP( QG )
##
## The function PrimitiveCentralIdempotentsByExtSSP computes the set of 
## primitive central idempotents of the group algebra QG, realizable by
## extremely strong Shoda pairs of G.
##
PrimitiveCentralIdempotentsByExtSSP:=function(QG)
local  G,    # underlying group of QG 
       PCIs, # the list of primitive central idempotents of QG
       ESSP, # output of the function ExtSSPAndDim(G)
       P,    # P in ESSP, an extremely strong Shoda pair of G
       IdP;  # the primitive central idempotents of QG associated to P

PCIs:=[];
G:=UnderlyingMagma(QG);
ESSP:=ExtSSPAndDim(G);
for P in ESSP.ExtremelyStrongShodaPairs do 
 IdP:=Idempotent_eGsum(QG,P[1],P[2])[2]; #idempotent of QG associated to P
 Add(PCIs,IdP);
od;
if not ESSP(G).SumDimension = Size(G) then 
 Print("Warning! The output is not complete list of pcis of the input! \n");
fi;
return PCIs;
end;


#############################################################################
## SearchingNNKForSSP(QG,H)
##
## The function SearchingNNKForSSP searches a non-normal subgroup K such 
## that (K,H) is a strong Shoda Pair of G and returns [ [ K, H ], e( G, K, H ) ] 
## or returns fail, if such K doesn't exist.
##
SearchingNNKForSSP:=function(QG,H)
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
     K:=PreImages(Epi,Cen);
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
       K:=PreImages(Epi,KH);
        if IsNormal(G,K)=false then
          return Idempotent_eGsum(QG,K,H);
        fi;
      fi;
       X:=Difference(X,KH);
    od;
  fi;
fi;
return fail;
end;


#############################################################################
## StSP(QG)
##
## The attribute StSP of the rational group algebra QG returns a record of 
## extremely strong Shoda pairs of G, if G is normally monomial. Otherwise, 
## this attribute returns a record of 3 components namely SSP, SId and SumDim
## where 
##
## SSP = list of strong Shoda pairs of G that cover the complete set of primitive
## central idempotents of QG realizable by strong Shoda pairs, 
##
## SId = set of primitive central idempotents of QG realizable by strong Shoda pairs,
## and
## SumDim = Sum of the Q-dimensions of simple components corresponding to strong 
## Shoda pairs of G.
##
StSP:=function(QG)
local  G,      # underlying group of QG
       ESSP,   # output of the function ExtSSPAndDim(G)
       SSP,    # set of representatives of strong Shoda pairs of G    
       SumDim, # sum of the Q-dimensions of simple algebras associated elements of SSP
       SId,    # set of idempotents associated to strong Shoda pairs in SSP
       P,      # P in SSP, a strong Shoda pair of G
       IdP,    # the primitive central idempotents of QG associated to P 
       Con,    # conjugacy classes of G;
       Con1,   # conjugacy classes of G of Size greater than 1;
       c,      # c in Con1
       r,      # representative of c
       Cr,     # Core of r in G
       K,      # non-normal subgroup of G corresponding to r
       NZ;     # normalizer 

#Initialization
G:=UnderlyingMagma(QG); 
ESSP:=ExtSSPAndDim(G); 
SSP:=ESSP.ExtremelyStrongShodaPairs; 
SumDim:=ESSP.SumDimension;

if SumDim=Size(G) then 
  return [SSP];       #if G is normally monomial, then set of strong Shoda pairs of G
                      #is set of extremely strong Shoda pairs of G
else 
 SId:=[];
  for P in SSP do
   IdP:=Idempotent_eGsum(QG,P[1],P[2])[2]; 
   Add(SId,IdP); #pcis associated with extremely strong Shoda pairs are added to SId
  od;  
   Con:=ConjugacyClassesSubgroups(G);
   Con1:=Filtered(Con,x->Size(x)>1);
#here, we pick a non-abelian subgroup from different conjugacy class
     for c in Con1 do
      r:=Representative(c);
      Cr:=Core(G,r);
       if IsCyclic(Centre(G/Cr)) then
        K:=SearchingNNKForSSP(QG,r);
         if not K=fail then
          P:=K[1]; 
          IdP:=K[2];
           if not IdP in SId then 
            Add(SId,IdP); 
            Add(SSP,P);
            NZ:=Normalizer(G,P[2]);
            SumDim:=SumDim+
            ((Size(G)/Size(NZ))^2)*((Phi(Size(P[1])/Size(P[2])))*(Size(NZ)/Size(P[1])));
           fi;
         fi;
       fi;  
       if SumDim=Size(G) then 
        break;
       fi;
     od;
 return [SSP,SId,SumDim];
fi;
end;


#############################################################################
## StShodaPairs(G)
##
## The function StShodaPairs computes a list of strongly Shoda pairs of the
## group G that covers the complete set of primitive central idempotents of the 
## rational group algebra QG realizable by strong Shoda pairs.
##
StShodaPairs:=function(G)
return StSP(GroupRing(Rationals,G))[1];
end;


#############################################################################
####PrimitiveCentralIdempotentsByStSP(QG)
##
## The function PrimitiveCentralIdempotentsByStSP computes the set of primitive
## central idempotents of the group algebra QG, realizable by strong Shoda pairs 
## of G.
##
PrimitiveCentralIdempotentsByStSP:=function(QG)
local  G,    # underlying group of QG 
       S,  # Output of StSP(QG) 
       PCIs, # the set of primitive central idempotents of QG
       P,    # P in SSP[1], a strong Shoda pair of G
       IdP;  # the primitive central idempotents of QG associated to P
       
G:=UnderlyingMagma(QG);
S:=StSP(QG);
if Size(S)=1 then 
 PCIs:=[];
  for P in S[1] do 
   IdP:=Idempotent_eGsum(QG,P[1],P[2])[2]; 
   Add(PCIs,IdP);      
  od; 
  return PCIs; 
else 
 PCIs:=S[2];
  if not S[3] = Size(G) then 
   Print("Warning! The output is not complete list of pcis of the input! \n"); 
  fi;
  return PCIs; 
fi;  
end;

###########################################################################
## E
###########################################################################