#############################################################################
##
#W  others.gi             The Wedderga package                Aurora Olivieri
#W                                                              Angel del Rio
#W                                                        Alexander Konovalov
##
#H  $Id$
##
#############################################################################

#############################################################################
##
#P  IsRationalGroupAlgebra( <R> )
##  
##  The function checks whether a group ring is a rational group algebra
InstallMethod( IsRationalGroupAlgebra,
    "for a group ring", 
    true,
    [IsGroupRing], 
    0,
    R -> IsRationals(LeftActingDomain(R)) 
); 

#############################################################################
##
##  The function PCIsFromShodaPairs computes the primitive central
##  idempotents of the form a e (G,K,H) where a is a rational number
##  and (H,K) is a SP. The sum of the PCIs obtained is 1 if and only
##  if G is monomial
##
#F PCIsFromShodaPairs( QG ) 
## 
InstallGlobalFunction(PCIsFromShodaPairs, 
function(QG) 
local   G,          #The group
        CCS,        #The conjugacy classes of subgroups
        LCCS,       #The lenght of CCS
        zero,       #0 of QG
        one,        #1 of QG
        Es,         #The list of primitive central idempotents  
        SEs,        #The sum of the elements of Es
        e,          #Idempotent
        H,          #Subgroup of G
        i,          #Counter
        SearchingKForSP;#Function to search a K for a given H  

#begin of functions 

#The following function search an element K such that 
#(K,H) is a SP
   
    SearchingKForSP:=function(H)
    local   
        NH,         #Normalizer of H in G
        Epi,        #NH --> NH/H        
        NHH,        #NH/H
        L,          #Centre of NHH
        K,          #The subgroup searched
        e,          #a.e(G,K,H) for some of the searched K
                    # and a rational
        KH,         #K/H
        X;          #a subset of NHH

        NH:=Normalizer(G,H);
        Epi:=NaturalHomomorphismByNormalSubgroup( NH, H ) ;
        NHH:=Image(Epi,NH);
        if IsAbelian(NHH) then 
            e:=PCIFromSP( QG, NH, H );
            if e<>false and SEs*e=zero then
                return e;
            fi;
        else 
            L:=Centre(NHH);
            if IsCyclic(L) then #This guaranties (S1) and (S2)
                X:=Difference(Elements(NHH),Elements(L));
                while X<>[] do
                    KH:=Subgroup(NHH,Union(L,[X[1]]));
                    K:=PreImages(Epi,KH);
                    e:=PCIFromSP( QG, K, H ); 
                    if e<>false and SEs*e=zero then
                        return e;
                    fi;
                    X:=Difference(X,KH);
                od;
            fi;
        fi;
    return false;                  
    end;     

#end of functions

#PROGRAM
    
#We start checking if QG is a rational group algebra of a 
#finite group
    if not IsRationalGroupAlgebra(QG) then
        Print("The input must be a rational group algebra \n");
        return fail;
    fi;

#Initialization

    G:=UnderlyingMagma(QG);
    CCS:=ConjugacyClassesSubgroups(G);
    LCCS:=Length(CCS);
    zero:=Zero(QG);
    one:=One(QG);
    Es:=[ ];
    SEs:=zero;
    i:=LCCS;   

#Main loop
    if Size(G)=1 then
        return [One(QG)];
    fi;
    while SEs<>one and i>=1 do 
        H:=Representative(CCS[i]); 
        e:=SearchingKForSP(H); 
        if e<>false then
            SEs:= SEs + e;
            Add(Es,e);
        fi;
        i:=i-1;
    od;
    
#Output
    if SEs<>one then 
        Print("Caution! It is not a complete set of primitive central idempotents \n");
    fi;
    return Es;
end);

#############################################################################
##
##  The function PCIsUsingCharacterTable uses the character table of G to compute the 
##  primitive central idempotents of QG with the classical method.
##
#F PCIsUsingCharacterTable( QG ) 
## 

InstallGlobalFunction( PCIsUsingCharacterTable, 
function(QG) 
local   G,      #The group
        ElemG,  #The elements of G
        OrderG, #Order of G
        zero,   #Zero of QG
        oneG,   #One of QG
        Emb,    #Embedding of G in QG
        I,      #The irreducible characters of G
        L,      #Length of I
        orbt,   #G-orbits of I 
        Lorbt,  #Length of orbt
        Nueva,  #Temporary value of (I minus the union of orbt)
        idem,   #Primitive Central Idempotent of QG
        eC,     #Primitive Central Idempotent of CG
        g,      #An element of G
        g1,     #g considered in QG
        Id,     #The list of primitive central idempotents of QG computed
        Qchi,   #Character Field
        GG,     #Gal(Qchi,Q)
        O,      #G-Orbit of a character of G
        char,   #The value of a character in a group element
        i, j, l;#Counters  

#First one check if QG is a rational group algebra of a finite group

    if not IsRationalGroupAlgebra(QG) then
        Print("The input must be a rational group algebra \n");
        return fail;
    fi;

#Initialization

    G:=UnderlyingMagma(QG);
    ElemG:=Elements(G);
    OrderG:=Size(G);
    
    zero:=Zero(QG);
    oneG:=One(G);
    Emb:=Embedding(G,QG);

#Computing the irreducible characters
    
    I:=Irr(G);
    L:=Length(I);    
    
#Computing the orbits
    
orbt:=[];
    
    Nueva:=I;
    while Nueva<>[] do
        Qchi:=CyclotomicField( Nueva[1]!.ValuesOfClassFunction);
        GG:=GaloisGroup( Qchi );
        O:=Orbit(GG , Nueva[1] ) ;
        AddSet( orbt, O);
        Nueva:=Difference(Nueva,O);
    od;
    Lorbt:=Length(orbt);
    
#Computing the PCIs of QG
    
    Id:=[];
    for i in [1..Lorbt] do
        O:=orbt[i];
        idem:=zero;
        for j in [1..Length(O)] do
            char:=O[j];
            
            #Computing PCI of CG
            
            eC:=zero;
            for l in [1..OrderG] do
                g:=ElemG[l];
                g1:=g^Emb;
                eC:=eC+ComplexConjugate(g^char)*g1;
            od;
            eC:=(oneG^char)*(OrderG^-1)*eC;
            
            #End of computation of PCI of CG
            
            idem:=idem+eC;
        od;
        AddSet(Id,idem);
    od;

return Id;
end);
 

#############################################################################
##
##  The function PCIsUsingConlon uses the function IrrConlon to compute the 
##  primitive central idempotents of QG associated to monomial representations
##  The result is the same as the function PCIsFromShodaPairs but it is slower
##
#F PCIsUsingConlon( QG )
##
InstallGlobalFunction(PCIsUsingConlon,
function(QG)
local   G, zero, one, IrrG, LIrrG, eGKHs, SeGKHs, i, eGKH, K, H; 

    if not IsRationalGroupAlgebra(QG) then
        Print("The input must be a rational group algebra \n");
        return fail;
    fi;

    
    G:=UnderlyingMagma(QG);
    zero:=Zero(QG);
    one:=One(QG);
    
    IrrG:=IrrConlon( G );
    LIrrG:=Length(IrrG);
    
    eGKHs:=[Epsilon(QG,G,G)];
    SeGKHs:=eGKHs[1];
    i:=2;
    while i<=LIrrG do
        K:=TestMonomial(IrrG[i])!.subgroup;
        H:=TestMonomial(IrrG[i])!.kernel;
        eGKH:=PCIFromSP(QG,K,H);    
        if eGKH*SeGKHs=zero then 
            SeGKHs:= SeGKHs + eGKH;
            Add(eGKHs,eGKH);
        fi;
        if SeGKHs=one then
            return eGKHs;
        fi;
    i:=i+1;
    od;
    if SeGKHs=one then
        return eGKHs;
    else 
        Print("Caution! It is not a complete set of primitive central idempotents \n");
        return eGKHs;
    fi;
end);
  
