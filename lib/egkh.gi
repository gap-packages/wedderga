#############################################################################
##
#W  egkh.gi               The Wedderga package                Aurora Olivieri
#W                                                              Angel del Rio
#W                                                        Alexander Konovalov
##
#H  $Id$
##
#############################################################################


#############################################################################
##
##  The function EpsilonCyclic compute epsilon(G,K,H) for H and K subgroups of G
##  such that H is normal in K and K/H is cyclic.
##  This function is faster than Epsilon
##
#M EpsilonCyclic( QG, K, H )
##
InstallMethod( EpsilonCyclic,"for pairs of subgroups", true, 
[IsFreeMagmaRing,IsGroup, IsGroup ], 0,
function(QG,K,H)
local   Hat,    #Function(H)=|H|^-1 sum_{h\in H}h
        L,      #subgrup of G
        G,      #Group
        Emb,    #Embedding of G in QG
        zero,   #0 of QG
        one,    #1 of QG
        Epsilon,#central idempotent 
        Epi,    #K--->K/H
        KH,     #K/H
        n,      #Order of KH
        y,      #Generator of KH 
        x,      #Representative of preimage of y
        p,      #Set of divisors of n
        Lp,     #Length od p
        i;      #Counter

#begin of functions 
#Hat computes the idempotent of $\Q G$ defined 
#by a subgroup
    Hat:=function(L)
    local LnA, SumL; 
        
        LnA:=Image(Emb,L);
        SumL:=Sum(LnA);
        return 1/(Size(L))*SumL;
    end;
#end of functions

#PROGRAM

#We start checking if QG is a rational group algebra of a 
#finite group

    if not(IsFreeMagmaRing(QG)) or not(IsGroup(QG!.UnderlyingMagma)) 
        or not(QG!.LeftActingDomain!.Name="Rationals") then

        Print("The first input must be a rational group algebra \n");
        return fail;    
    fi;
#Then we check if K is subgroup of G, H is a normal
#subgroup of K and K/H is cyclic
    if not(IsSubgroup(QG!.UnderlyingMagma,K)) then
        Print("The group algebra does not correspond to the subgroups \n");
        return fail;
    elif not(IsSubgroup(K,H) or IsNormal(K,H)) then
        Print("The second subgroup must be a normal subgroup in the first one \n");
        return fail;
    elif not(IsCyclic(FactorGroup(K,H))) then
        Print("The factor group (second input over third one) is not a cyclic group \n");
        return fail;
    fi;

#Initialization
#Core of the program  
    G:=QG!.UnderlyingMagma;
    Emb:= Embedding( G, QG );
    Epsilon:=Hat(H);
    one:=One(QG);
    if K<>H then
        Epi:=NaturalHomomorphismByNormalSubgroup( K, H ) ;
        KH:=Image(Epi,K); 
        n:=Size(KH);
        y:=Product(IndependentGeneratorsOfAbelianGroup(KH)); 
        x:=PreImagesRepresentative(Epi,y); 
        p:= Set(FactorsInt(n));
        Lp:=Length(p);
        for i in [1..Lp] do
            Epsilon:=Epsilon*(one-Hat(Group(x^(n/p[i]))));
        od;
    fi;

#Output
return Epsilon; 
end);

#############################################################################
##
##  The function Epsilon computes epsilon(K,H) for H a normal subgroup of K
## 
#M Epsilon( QG, K, H )
##
InstallMethod( Epsilon,"for pairs of subgroups", true, 
[IsFreeMagmaRing,IsGroup, IsGroup ], 0,
function(QG,K,H)
local   Hat, L, G, Emb, hatH, Epsilon, Epi, KH, SNKH, LSNKH, i;

    if not(IsFreeMagmaRing(QG)) or not(IsGroup(QG!.UnderlyingMagma)) or 
        not(QG!.LeftActingDomain!.Name="Rationals") then

        Print("The first input must be a rational group algebra \n");
        return fail;    
    fi;
    if not(IsSubgroup(QG!.UnderlyingMagma,K)) then
        Print("The group algebra does not correspond to the subgroups \n");
        return fail;
    elif not(IsSubgroup(K,H))  then
        Print("The second subgroup must be contained in the first one \n");
        return fail;
    elif not(IsNormal(K,H)) then
        Print("The second subgroup must be normal in the first one \n");
        return fail;
    fi;

    Hat:=function(L)
    local LnA, SumL; 
        
        LnA:=Image(Emb,L);
        SumL:=Sum(LnA);
        return 1/(Size(L))*SumL;
    end;
    
    G:=QG!.UnderlyingMagma;
    Emb:= Embedding( G, QG );
    hatH:=Hat(H);
    Epsilon:=hatH;
    if K<>H then
        Epi:=NaturalHomomorphismByNormalSubgroup( K, H ) ;
        KH:=Image(Epi,K);
        SNKH:=NormalSubgroups(KH);
        LSNKH:=Length(SNKH); 
        for i in [1..LSNKH] do
            if IsPrime(Size(SNKH[i])) then
                L:=PreImage(Epi,SNKH[i]);
                Epsilon:=Epsilon*(hatH-Hat(L));
            fi;
        od;
    fi;

return Epsilon; 
end);

#############################################################################
##  The function CentralizerG computes the centralizer in G of an element of QG.
##  It is used to compute e(G,K,H) 
##
#M CentralizerG( QG, alpha )
##
#InstallMethod( CentralizerG,"for elements of rational group algebra", true, 
#[IsFreeMagmaRing,IsElementOfFreeMagmaRing ], 0,
#function(QG,alpha)
#local G, Emb, ElemG, OrderG, gen, i, g, g1;
#
#    if not(IsFreeMagmaRing(QG)) or not(IsGroup(QG!.UnderlyingMagma)) or 
#        not(QG!.LeftActingDomain!.Name="Rationals") then
#
#        Print("The first input must be a rational group algebra \n");
#        return fail;    
#    fi;
#    if not(alpha in  QG) then
#
#        Print("The second input must be an element of the rational group algebra \n");
#        return fail;    
#    fi;
#    
#    G:=QG!.UnderlyingMagma;
#    Emb:= Embedding( G, QG );    
#    #G1:=Image(Emb,G);
#    ElemG:=Elements(G);
#    OrderG:=Size(G);
#    gen:=[];
#    for i in [2..OrderG] do
#        g:=ElemG[i];
#        g1:=g^Emb;
#        if g1^(-1)*alpha*g1=alpha then
#            Add(gen,g);
#        fi;
#    od;
#return Subgroup(G,gen);
#end);

#############################################################################
##
##  The function eGKH computes e(G,K,H) for H and K subgroups of G such that H 
##  is normal in K
##
#M eGKH( QG, K, H )
##
InstallMethod( eGKH,"for pairs of subgroups", true, 
[IsFreeMagmaRing,IsGroup, IsGroup ], 0,
function(QG,K,H)
local   CentralizerG, alpha, G, Emb, zero, Eps, Cen, eGKH, RTCen, nRTCen, i, g;

    if not(IsFreeMagmaRing(QG)) or not(IsGroup(QG!.UnderlyingMagma)) or 
        not(QG!.LeftActingDomain!.Name="Rationals") then

        Print("The first input must be a rational group algebra \n");
        return fail;    
    fi;
    
    if not(IsSubgroup(QG!.UnderlyingMagma,K)) then
        Print("The group algebra does not correspond to the subgroups \n");
        return fail;
    elif not(IsSubgroup(K,H) or IsNormal(K,H)) then
        Print("The second subgroup must be normal in the first one \n");
        return fail;
    fi;

    CentralizerG:=function(alpha)
    local   ElemG,  # Elements of G
            OrderG, #|G|
            gen,    #generators of Centraliazer
            i,      #Counter
            g;      #g^Emb
    
        ElemG:=Elements(G);
        OrderG:=Size(G);
        gen:=[ ];
        for i in [2..OrderG] do
            g:=ElemG[i]^Emb;
            if g^(-1)*alpha*g=alpha then
                Add(gen,ElemG[i]);
            fi;
        od;
    return Subgroup(G,gen);
    end;

    G:=QG!.UnderlyingMagma;
    Emb:= Embedding( G, QG );
    Eps:=Epsilon(QG,K,H);
    if Eps in QG then
        Cen:=CentralizerG(Eps);
        eGKH:=Eps;
        if Cen<>G then
            RTCen:=RightTransversal(G,Cen); 
            nRTCen:=Length(RTCen);  
            for i in [2..nRTCen] do
                g:=RTCen[i]^Emb;
                eGKH:=eGKH+ Eps^g;    
            od;
        fi;
    else
        return fail;
    fi;

return eGKH;
end);

#############################################################################
##
##  The function VerifyShoda verifies if (H,K) is a SP
##
## 
#M VerifyShoda(G, K, H )
##
#InstallMethod( VerifyShoda,"for pairs of subgroups", true, 
#[IsGroup, IsGroup, IsGroup ], 0,
#function(G,K,H)
#    local   DGK, LDGK, ElemK, OrderK, KH, i, g, NoReady, j, k; 
#    
#    if not(IsSubgroup(G,H)) then 
#    Print("The second input must be a subgroup in the first one \n");
#        return fail;
#    elif not(IsSubgroup(K,H)) then
#        Print("The third subgroup must be contained in the second one \n");
#        return fail;
#    elif not(IsNormal(K,H)) then
#        Print("The third subgroup is not normal in the second one \n");
#        return false;        
#    elif not(IsCyclic(FactorGroup(K,H))) then
#        Print("The factor group (second input over first one) is not cyclic \n");
#        return false;
#    fi;
#
#DGK:=Difference(G,K);
#LDGK:=Length(DGK);
#ElemK:=Elements(K);
#OrderK:=Size(K);
#KH:=Difference(K,H);
#i:=1;
#while i<=LDGK do        
#    g:=DGK[i];
#    NoReady:=true;
#    j:=1;
#    while j<=OrderK and NoReady do
#    k:=ElemK[j];
#        if Comm(k,g) in KH then
#            NoReady:=false;
#        fi;
#    j:=j+1;
#    od;
#    if j=(OrderK+1) and NoReady then
#        return false;
#    fi;
#i:=i+1;
#od;
#return true;
#end);

#############################################################################
##
##  The function PCIFromSP checks if (H,K) is a SP of G and in that
##  case compute the primitive central idempotent of the form a e(G,K,H) where 
##  a is a rational number.
##
#M PCIFromSP( QG, K, H )
##
InstallMethod( PCIFromSP,"for pairs of subgroups", true, 
[IsFreeMagmaRing,IsGroup, IsGroup ], 0,
function(QG,K,H)
local   eG,             #Function for eGKH
        VerifyShoda,    #Function for Shoda
        K1, H1,         #Subgroups of G
        G,              #The group
        Emb,            # G --> QG     
        e1,             #eGKH
        e2,             #eGKH^2
        B,              #basis
        ce1,            #Coeficients of e1
        ce2,            #Coeficients of e1
        alpha;          #ce1[1]/ce2[1]

#begin of functions    

#The following function computes e(G,K,H)    

    eG:=function(K1,H1)
    local   CentralizerG,  #Function
            alpha,          #Element of QG   
            Eps,    # \varepsilon(K1,H1), 
            Cen,     # Centralizer of Eps in G
            RTCen,   # Right transversal of Cen in G
            g;      # element of G

#The following function computes the centralizer in G of an 
#element of QG.
        CentralizerG:=function(alpha)
        local   ElemG,  # Elements of G
                OrderG, #|G|
                gen,    #generators of Centraliazer
                i,      #Counter
                g;      #g^Emb
        
            ElemG:=Elements(G);
            OrderG:=Size(G);
            gen:=[ ];
            for i in [2..OrderG] do
                g:=ElemG[i]^Emb;
                if g^(-1)*alpha*g=alpha then
                    Add(gen,ElemG[i]);
                fi;
            od;
        return Subgroup(G,gen);
        end;

        
        Eps:=EpsilonCyclic(QG,K1,H1); #ver C.2.1
        Cen:=CentralizerG( Eps );
        RTCen:=RightTransversal(G,Cen);
    return Sum(List(RTCen,g->Eps^(g^Emb)));
    end;
   
#The function VerifyShoda verifies (S3) for a pair (H,K) 
#which satisfies (S1) and (S2) for K/H cyclic

    VerifyShoda:=function(K1,H1)
    local   DGK,        # G\K1 
            LDGK,       #|DGK|
            ElemK,      #Elements of K1
            OrderK,     #|K1|
            KH,         #K1\H1
            i,j,        #Counters
            k,          #Element of K1
            g,          #Elements of DGK
            NoReady;    # Boolean 
        DGK:=Difference(G,K1);
        LDGK:=Length(DGK);
        ElemK:=Elements(K1);
        OrderK:=Size(K1);
        KH:=Difference(K1,H1);
        i:=1;
        while i<=LDGK do        
            g:=DGK[i];
            NoReady:=true;
            j:=1;
            while j<=OrderK and NoReady do
            k:=ElemK[j];
                if Comm(k,g) in KH then
                    NoReady:=false;
                fi;
            j:=j+1;
            od;
            if j=(OrderK+1) and NoReady then
                return false;
            fi;
        i:=i+1;
        od;
    return true;
    end;

#end of functions

#PROGRAM
    
#We start checking if QG is a rational group algebra of a 
#finite group

    if not(IsFreeMagmaRing(QG)) or not(IsGroup(QG!.UnderlyingMagma)) 
        or not(QG!.LeftActingDomain!.Name="Rationals") then

        Print("The first input must be a rational group algebra \n");
        return fail;    
    fi;
#Checking (S1) and (S2)
    if not(IsSubgroup(QG!.UnderlyingMagma,K)) then
        Print("The group algebra does not correspond to the subgroups \n");
        return fail;
    elif not(IsSubgroup(K,H)) then
        Print("The second subgroup must be a subgroup in the first one \n");
        return fail;
    elif not(IsNormal(K,H)) then
        return false;
    elif not(IsCyclic(FactorGroup(K,H))) then
        return false;
    fi;

#Initialization
#Checking (S3) and computing the idempotent


    G:=QG!.UnderlyingMagma;
    Emb:= Embedding( G, QG );    
    
    if VerifyShoda(K,H)=true then
        e1:=eG(K,H);
        e2:=e1^2;
        if e2=e1 then
            return  e1;
        else
            B:=Basis(QG);
            ce1:=Coefficients(B,e1);
            ce2:=Coefficients(B,e2);
            alpha:=ce1[1]/ce2[1];
        fi;
    else 
        return false;
    fi;
#Here finish the main loop

return alpha*e1;
end);

#############################################################################
##
##  The function IsCompleteSetOfPCIs checks if the sum of the elements of QG is 1
##  It is suppose to be used to check if a list of PCIs of QG is complete.
##
#M IsCompleteSetOfPCIs( QG, ListPCIs )
##

InstallMethod( IsCompleteSetOfPCIs,"for list of primitive central idempotents", true, 
[IsFreeMagmaRing,IsList ], 0,
function( QG, ListPCIs )
    local i;
    if not(IsFreeMagmaRing(QG)) or not(QG!.LeftActingDomain!.Name="Rationals") then

        Print("The first input must be a rational group algebra \n");
        return fail;
    elif not(IsList(ListPCIs)) then
        Print("The second input must be a list of PCIs \n");
        return fail;
    fi;
    for i in [1..Length(ListPCIs)] do
        if not(ListPCIs[i] in QG) then
    
            Print("The ",i,"-th element of the second input should belong to the first input\n");
            return fail;    
        fi;
    od;    

    if Sum(ListPCIs)=One(QG) then
        return true;
    else
        return false;
    fi;
end);
