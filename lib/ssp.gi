#############################################################################
##
#W  ssp.gi                The Wedderga package                Aurora Olivieri
#W                                                              Angel del Rio
#W                                                        Alexander Konovalov
##
#H  $Id$
##
#############################################################################


#############################################################################
## 
## eG:=function(QG,K,H)
##
## The following function computes e(G,K,H)    
## Note that actually the 1st argument is QG, not G
## It returns a list of the form [eGKH,[K,H]]
##
eG:=function(QG,K,H)
    local   
        G,      # underlying group of QG
        Eps,    # \varepsilon(K,H), 
        eGKH,   # is the final return that takes partial values 
                # with partial sums of conjugates of Eps  
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
        if NH<>G then
            NdK:=Normalizer(G,K);
            eGKH1:=Eps;
            if NH<>NdK then
                RTNH:=RightTransversal(NdK,NH); 
                nRTNH:=Length(RTNH);  
                for i in [2..nRTNH] do
                    g:=RTNH[i]^Embedding( G, QG );
                    eGKH1 := eGKH1 + Eps^g;    
                od;
            fi;
            eGKH:=eGKH1;        
            if NdK<>G then
                RTNdK:=RightTransversal(G,NdK); 
                nRTNdK:=Length(RTNdK);  
                for i in [ 2 .. nRTNdK ] do
                    g:=RTNdK[i]^Embedding( G, QG );
                    eGKH1g:=eGKH1^g;
                    if eGKH1*(eGKH1g) <> Zero(QG) then    
                        return  fail;
                    else
                        eGKH:= eGKH + eGKH1g;
                    fi;
                od;    
            fi; 
        else 
            return [Eps,[K,H]];
        fi;       
    return [eGKH,[K,H]];
    end;


#############################################################################
## 
## SearchingKForSSP(QG,H)
##
## The following function search an element K such that (K,H) is a SSP
## and returns eG(QG,K,H) in the form [eGKH,[K,H]] or returns fail, if
## such K doesn't exist
##
SearchingKForSSP:=function(QG,H)
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
    end;               


#############################################################################
## 
## The function PCIsFromSSP computes the primitive central idempotents of 
## the form e(G,K,H) for the SSPs(H,K) of a finite group G
##
#F PCIsFromSSP( QG )
##
InstallGlobalFunction( PCIsFromSSP, function(QG)
local   G,          # The group
        DG,         # Derived subgroup of G
        CCS,        # The conjugacy classes of subgroups
        LCCS,       # The lenght of CCS 
        eGKHs,      # The list of primitive central idempotents 
        SeGKHs,     # The sum of the elements of eGKHs
        H,          # Subgroup of G
        NH,         # Normalizer of H in G
        NHH,        # NH/H
        K,          # Subgroup of G 
        i, j,       # Counters
        eGKH,       # A primitive central idempotent of the form 
                    # e(G,K,H) for (K,H) a SSP
        KH;         # K/H

#We start checking if QG is a rational group algebra of a finite group

#Actually finiteness was not checked in the old code - ask Angel !!!

    if not IsRationalGroupAlgebra(QG) then
      Error("The input must be a rational group algebra !!!");
    fi;

#Initialization

    G := UnderlyingMagma(QG);
    CCS:=ConjugacyClassesSubgroups(G);
    LCCS:=Length(CCS);
    DG:=DerivedSubgroup(G); 
    eGKHs:=[];
    SeGKHs:=Zero(QG);

#Main loop

    if Size(G)=1 then
        return [One(QG)];
    fi;
    
    j:=LCCS;    

    while SeGKHs<>One(QG) and j>=1 do 
    
        H:=Representative(CCS[j]);        
        if IsSubset(H,DG) then
            if IsCyclic(FactorGroup(G,H)) then 
                eGKH:=eG(QG,G,H)[1]; 
                SeGKHs:= SeGKHs + eGKH;
                Add(eGKHs,eGKH);
            fi;
        else 
            eGKH:=SearchingKForSSP(QG,H);
            if eGKH <> fail and eGKH[1]*SeGKHs=Zero(QG) then 
                SeGKHs:= SeGKHs + eGKH[1];
                Add(eGKHs,eGKH[1]);
            fi;
        fi;
        j:=j-1;
    od;

#Output

    if SeGKHs <> One(QG) then 
        Print("Caution! It is not a complete set of primitive central idempotents \n");
    fi;
    return eGKHs;

end);


#############################################################################
##
#F StronglyShodaPairs( QG )
##
## The function StronglyShodaPairs computes a set of SSPs (H,K) 
## that covers the primitive central idempotents obtained from SSPs
##
InstallGlobalFunction( StronglyShodaPairs, function(QG)
local   G,          # The group
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
        eGKH,       # A primitive central idempotent of the form e(G,K,H) for (H,K) a SSP
        KH,         # K/H
        KHs;        # The list of SSP

#We start checking if QG is a rational group algebra of a finite group

    if not IsRationalGroupAlgebra(QG) then
        Error("The input must be a rational group algebra !!!");
    fi;

#Initialization

    G:=UnderlyingMagma(QG);
    CCS:=ConjugacyClassesSubgroups(G);
    LCCS:=Length(CCS);
    DG:=DerivedSubgroup(G); 
    KHs:=[];
    SeGKHs:=Zero(QG);
    
    if Size(G)=1 then
        return [[G,G]];
    fi;

    j:=LCCS;    
    while SeGKHs<>One(QG) and j>=1 do 
        H:=Representative(CCS[j]);        
        if IsSubset(H,DG) then
            if IsCyclic(FactorGroup(G,H)) then 
                eGKH:=eG(QG,G,H)[1]; 
                SeGKHs:= SeGKHs + eGKH;
                Add(KHs,[G,H]);
            fi;
        else 
            eGKH:=SearchingKForSSP(QG,H);
            if eGKH<>fail and eGKH[1]*SeGKHs=Zero(QG) then 
                SeGKHs:= SeGKHs + eGKH[1];
                Add(KHs,eGKH[2]);
            fi;
        fi;
    j:=j-1;
    od;

#Here finish the main loop

    if SeGKHs<>One(QG) then 
        Print("Caution! Some primitive central idempotents are not realizable by strongly Shoda pairs \n");
    fi;
    return KHs;
end);


#############################################################################
##
##  The function eGKHFromSSP verifies if (H,K) is a SSP and in that
##  case computes e(G,K,H).
##
#F eGKHFromSSP( QG, K, H )
##
InstallGlobalFunction(eGKHFromSSP, 
function(QG,K,H)
local   Verify, K1, H1, G, Emb, zero, NH, Eps,  NdK, eGKH1, RTNH, nRTNH, i, g, eGKH, RTNdK, nRTNdK, eGKH1g;

    if not IsRationalGroupAlgebra(QG) then
        Error("The first input must be a rational group algebra !!!");
    fi;
    
    if not IsSubgroup(UnderlyingMagma(QG),K) then
        Error("The group algebra does not correspond to the subgroups !!!");
    elif not(IsSubgroup(K,H) and IsNormal(K,H)) then
        Error("The second subgroup must be normal in the first one !!!");
    fi;

Verify:=function(K1,H1)
# verifies if H1 is a normal subgroup of K1 and
# K1/H1 is the maximal abelian subgroup in N1/H1,
# where N1 is the normalizer of H1 in G
local Epi, NHH, KH;

    if not(IsNormal(NH,K1)) then
        Print("The first subgroup must be normal in the normalizer of second one in G\n");
        return fail;
    fi;
    if K1<>NH then
        Epi:=NaturalHomomorphismByNormalSubgroup( NH, H1 ) ;
        NHH:=Image(Epi,NH); #It is isomorphic to the factor group NH/H.
        KH:=Image(Epi,K1); #It is isomorphic to the factor group K/H.
        if Centralizer(NHH,KH)<>KH then
            Print("The factor group (second input over third one) is not maximal abelian on the normalizer \n");
            return fail;
        fi;
    fi;
return true;
end;

G:=UnderlyingMagma(QG);
Emb:= Embedding( G, QG );
zero:=Zero(QG); 
NH:=Normalizer(G,H);
if Verify(K,H)=true then        
    Eps:=Epsilon(QG,K,H);
    if NH<>G and Eps in QG then
        NdK:=Normalizer(G,K);
        eGKH1:=Eps;
        if NH<>NdK then
            RTNH:=RightTransversal(NdK,NH); 
            nRTNH:=Length(RTNH);  
            for i in [2..nRTNH] do
                g:=RTNH[i]^Emb;
                eGKH1:=eGKH1+ Eps^g;    
            od;
        fi;
        eGKH:=eGKH1;        
        if NdK<>G then
            RTNdK:=RightTransversal(G,NdK); 
            nRTNdK:=Length(RTNdK);  
            i:=2;
            while i<=nRTNdK do
                g:=RTNdK[i]^Emb;
                eGKH1g:=eGKH1^g;
                if eGKH1*(eGKH1g)<>zero then    
                    Print("The conjugates of epsilon are not orthogonal \n");
                    return fail;
                else
                    eGKH:= eGKH + eGKH1g;
                fi;
            i:=i+1;
            od;    
        fi; 
    else 
    return Eps;
    fi;       
return eGKH;
else
    return fail;
fi;
end);


#############################################################################
##
##  The function eGKHsFromKHs computes the list of e(G,K,H)'s for (K,H) running
##  on a list of pairs of subgroups of G such that H is a normal subgroup of G
##  It is assumed to be used on a list of pairs (K,H) such that (H,K) is a strongly
##  SP
##
#F eGKHsFromKHs( QG, ListPairsOfSubg )
##
InstallGlobalFunction(eGKHsFromKHs, 
function(QG,list)
local   G, i;
    if not IsRationalGroupAlgebra(QG) then
        Error("The first input must be a rational group algebra !!!");
    elif not IsList(list) then
        Error("The second input must be a list (of pairs of subgroups)!!!");
    fi;
    G:=UnderlyingMagma(QG);
    for i in [1..Length(list)] do
        if not(IsSubgroup(G,list[i][1])) or not(IsSubgroup(G,list[i][2])) then
            Error("The ", Ordinal(i), 
            " element of the list is not a pair of subgroups of ", G, "!!!");
        elif ( not IsNormal( list[i][1], list[i][2]) ) or 
             ( not IsSubset( list[i][1], list[i][2]) ) then
            Error("The second entry of the ", Ordinal(i),
            " element of the list is not \n a normal subgroup of the ",
            " first entry of this element !!!");
        fi;
    od;    

return List([1..Length(list)], i->eGKH(QG,list[i][1],list[i][2]));
end);


#############################################################################
##
##  The function SimpleAlgebraFromSSP computes the data describing the simple
##  algebra QG e(G,K,H) for (H,K) a SSP of G
##
#F SimpleAlgebraFromSSP( G, K, H )
##
InstallGlobalFunction( SimpleAlgebraFromSSP, 
function(G,K,H)
local   N,          # Normalizer of H in G
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
                    # N/K
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
        x,i,j,ll;   # Controlers
        
    if G=H then 
        return [1,1,[],[]];
    else
    
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

# Output
        gen:=List([1..Length(Gen)],
                    i->PreImagesRepresentative(Epi2,Gen[i]));
        fi;
        return [Size(G)/Size(N),ok,
                List([1..Length(Gen)],
                i->[Order(Gen[i]),
                RemInt(Position(Potk,k^gen[i]),ok),
                RemInt(Position(Potk,gen[i]^Order(Gen[i])),ok)]),
                List([1..Length(Gen)-1],i->List([i+1..Length(Gen)],
                j->RemInt(Position(Potk,Comm(gen[j],gen[i])),ok)))];
fi;    
end);
    
 
#############################################################################
##
##  The function SimpleFactorsFromListOfSSP computes the data describing the simple
##  algebras QG e(G,K,H) for (H,K) running on a list of SSP of G
##
#F SimpleFactorsFromListOfSSP( G, ListPairsOfSubg )
##
InstallGlobalFunction(SimpleFactorsFromListOfSSP, 
function(G,list)
local   i;

    for i in [1..Length(list)] do
        if not(IsSubgroup(G,list[i][1])) or not(IsSubgroup(G,list[i][2])) then
            Error("The ", Ordinal(i), 
              " element of the list is not a pair of subgroups of ",G,"!!!");
        elif ( not IsNormal( list[i][1], list[i][2] ) ) or 
             ( not IsSubset( list[i][1], list[i][2] ) ) then
            Error("The first entry of the ", Ordinal(i), 
              " element of the list is not \n a normal subgroup ",
              " of the second entry of this element !!!");
        fi;
    od; 
return List([1..Length(list)], i->SimpleAlgebraFromSSP(G,list[i][1],list[i][2]));
end);