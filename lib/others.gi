#############################################################################
##
#W  others.gi              The Wedderga package           Osnel Broche Cristo
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                              Ángel del Río
##
#H  $Id$
##
#############################################################################


#############################################################################
##
#A ShodaPairsAndIdempotents( QG )
##
## The attribute ShodaPairsAndIdempotents of the rational group algebra QG 
## returns a record with components ShodaPairs and PCIsBySP
## ShodaPairs = list of SP that covers the complete set of primitive 
##       central idempotents of QG realizable by SPs, 
## PCIsBySP = list of PCIs of QG realizable by SPs.
## 
InstallMethod( ShodaPairsAndIdempotents, 
    "for rational group algebra", 
    true, 
    [ IsSemisimpleRationalGroupAlgebra ], 
    0,
function(QG) 

local   G,          #The group
        CCS,        #The conjugacy classes of subgroups
        LCCS,       #The length of CCS
        one,        #1 of QG
        Es,         #The list of primitive central idempotents  
        SEs,        #The sum of the elements of Es
        e,          #Idempotent
        H,          #Subgroup of G
        i,          #Counter
        SearchingKForSP,#Function to search a K for a given H  
        SPs;

#begin of functions 

#The following function search an element K such that 
#(K,H) is a SP
   
    SearchingKForSP:=function( H )
    local   
        NH,         #Normalizer of H in G
        Epi,        #NH --> NH/H        
        NHH,        #NH/H
        L,          #Centre of NHH
        K,          #The subgroup searched
        e,          #a*e(G,K,H) for some of the searched K
                    # and a rational
        KH,         #K/H
        X;          #a subset of NHH

    NH:=Normalizer(G,H);
    Epi:=NaturalHomomorphismByNormalSubgroup( NH, H ) ;
    NHH:=Image(Epi,NH);
    if IsAbelian(NHH) then
        e := PrimitiveCentralIdempotentBySP( QG, NH, H );
        if e=fail then
            return fail; 
        else 
            return [[NH,H],e];
        fi;
    else 
        L:=Centre(NHH);
        if IsCyclic(L) then #This guaranties (S1) and (S2)
            X:=Difference(Elements(NHH),Elements(L));
            while X<>[] do
                KH:=Subgroup(NHH,Union(L,[X[1]]));
                K:=PreImages(Epi,KH);
                e:=PrimitiveCentralIdempotentBySP( QG, K, H ); 
                if e<>fail then
                    return [[K,H],e];
                fi;
                X:=Difference(X,KH);
            od;
        fi;
        return fail;                  
    fi;
    end;     

#end of functions

#PROGRAM
    
#We start checking if QG is a rational group algebra of a finite group

    if not IsSemisimpleRationalGroupAlgebra(QG) then
        Error("Wedderga: The input must be a rational group algebra!!!\n");
    fi;

#Initialization

    G:=UnderlyingMagma(QG);
    CCS:=ConjugacyClassesSubgroups(G);
    LCCS:=Length(CCS);
    one:=One(QG);
    Es:=[ ];
    SPs:=[];
    SEs:=Zero(QG);
    i:=LCCS;   

#Main loop
    if Size(G)=1 then
        return [[ [G,G] , One(QG)] ];
    fi;
    while SEs<>one and i>=1 do 
        H:=Representative(CCS[i]); 
        e:=SearchingKForSP( H ); 
        if e<>fail then
                if IsZero( SEs*e[2] ) then # if SEs*e <> zero then
                SEs:= SEs + e[2];
                Add(Es,e[2]);
                Add(SPs,e[1] );
            fi;    
        fi;
        i:=i-1;
    od;
    
#Output

  return rec( 
    ShodaPairs := SPs, 
    PCIsBySP := Es);

end);




#############################################################################
##
#F PrimitiveCentralIdempotentsBySP( QG )
##
## The function computes the primitive central idempotents of the form 
## a*e(G,K,H) where a is a rational number and (H,K) is a SP. 
## The sum of the PCIs obtained is 1 if and only if G is monomial
## This function is for rational group algebras (otherwise you will not
## use ShodaPairsAndIdempotents)
##
InstallGlobalFunction(PrimitiveCentralIdempotentsBySP, 
function(QG) 
local G;

G := UnderlyingMagma( QG );
if not(IsMonomial(G)) then 
   Print("Wedderga: Warning!!\nThe output is a NON-COMPLETE list of prim. central idemp.s of the input! \n");
fi;

return ShodaPairsAndIdempotents( QG ).PCIsBySP; 
end);


#############################################################################
##
#M PrimitiveCentralIdempotentBySP( QG, K, H )
##
## The function PrimitiveCentralIdempotentBySP checks if (H,K) is a SP of G 
## and in that case compute the primitive central idempotent of the form 
## alpha*e(G,K,H), where alpha is a rational number.
##
InstallMethod( PrimitiveCentralIdempotentBySP,
   "for pairs of subgroups", 
   true, 
   [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ],
   0,
function(QG,K,H)
local   eGKH,       # Function for eGKH
        K1, H1,     # Subgroups of G
        G,          # The group
        e1,         # eGKH
        e2,         # eGKH^2
        alpha;      # coeff e1 / coeff ce2

#begin of functions    

#The following function computes e(G,K,H)    

    eGKH:=function(K1,H1)
    local   alpha,  # Element of QG   
            Eps,    # \varepsilon(K1,H1), 
            Cen,    # Centralizer of Eps in G
            RTCen,  # Right transversal of Cen in G
            g;      # element of G
        
    Eps := IdempotentBySubgroups( QG, K1, H1 );
    Cen := Centralizer( G, Eps );
    RTCen := RightTransversal( G, Cen );

    return Sum( List( RTCen, g -> Eps^g ) ) ;
    end;
   
#end of functions       

#Program    

    G:=UnderlyingMagma(QG);
  
    if IsShodaPair( G, K, H ) then
        e1:=eGKH(K,H);
        e2:=e1^2;
        alpha := CoefficientsAndMagmaElements(e1)[2] / 
                 CoefficientsAndMagmaElements(e2)[2];        
        if IsOne( alpha ) then
            return e1;
        else    
            return alpha*e1;      
        fi;    
    else
        return fail;
    fi;            
end);        





#############################################################################
##
#M IsShodaPair( G, K, H )
##
## The function IsShodaPair verifies if (H,K) is a SP 
##
InstallMethod( IsShodaPair,
    "for pairs of subgroups", 
    true, 
    [ IsGroup, IsGroup, IsGroup ],
    0,
function( G, K, H )
local       DGK,        # G\K
            ElemK,      # Elements of K
            DKH,        # K\H
            k,          # Element of K
            g,          # Elements of DGK
            NoReady;    # Boolean

# Checking (S1) and (S2)
# First verifies if H, K are subgroups of G and K is a normal subgroup of K
if not ( IsSubgroup( G, K ) and IsSubgroup( K, H ) and IsNormal( K, H ) ) then
    Error("Wedderga: Each input should contain the next one and the last one \n",
          "should be normal in the second!!!\n");
elif not( IsCyclic( FactorGroup( K, H ) ) ) then
    return false;
fi;

#Now, checking (S3)
DGK := Difference( G, K );
ElemK := Elements( K );
DKH := Difference( K, H );
for g in DGK  do
    NoReady := true;
    for k in ElemK do
        if Comm( k, g ) in DKH then
            NoReady := false;
            break;
        fi;
    od;
    if NoReady then
        return false;
    fi;
od;

return true;
end);




#############################################################################
##
#F PrimitiveCentralIdempotentsByCharacterTable( QG )
##
## The function PrimitiveCentralIdempotentsByCharacterTable 
## uses the character table of G to compute the primitive 
## central idempotents of QG with the classical method.
##
InstallGlobalFunction( PrimitiveCentralIdempotentsByCharacterTable, 
function(QG) 
local G,      # The group
      OrderG, # Order of G
      zero,   # Zero of QG
      I,      # The irreducible characters of G
      rat,    # rational irreducible characters (Galois orbit sums)
      norms,  # norms of the orbit sums
      Id,     # The list of primitive central idempotents of QG computed
      nccl,
      i,
      chi,    # one entry in `rat'
      eC,     # one primitive central idempotent of QG
      j,      # loop over class positions
      F, c, elms, val, k, tbl, fusions, deg;

    # First check if `QG' is a rational group algebra of a finite group.
    if not ( IsFreeMagmaRing( QG ) and IsGroup( UnderlyingMagma( QG ) )
             and IsRationals( LeftActingDomain( QG ) ) ) then
      Error( "The input must be a rational group algebra" );
    fi;

    # Initialization
    G:= UnderlyingMagma( QG );
    elms:= Elements( G );
    OrderG:= Size( G );
    zero:= Zero( QG );
		       
    # Compute the irreducible characters.
    IsSupersolvable( G );
    tbl:= CharacterTable( G );
    I:= List( Irr( G ), ValuesOfClassFunction );

    # Compute the Galois orbit sums.
    rat:= RationalizedMat( I );
    norms:= List( rat, x -> ScalarProduct( tbl, x, x ) );

    # Compute the PCIs of QG.
    F:= FamilyObj( zero );
    fusions:= List( ConjugacyClasses( tbl ),
                    c -> List( Elements( c ),
                               x -> PositionSorted( elms, x ) ) );
    Id:= [];
    nccl:= Length( I );    
    for i in [ 1 .. Length( rat ) ] do
      chi:= rat[i];
      deg:= chi[1] / norms[i];
      eC:= 0 * [ 1 .. OrderG ] + 1;
      for j in [ 1 .. nccl ] do
        val:= chi[j] * deg / OrderG;
        for k in fusions[j] do
          eC[k]:= val;
        od;
      od;
      Add( Id, ElementOfMagmaRing( F, 0, eC, elms ) );
    od;

    return Id;
end);
 

#############################################################################
##
#F PrimitiveCentralIdempotentsUsingConlon( QG )
##
##  The function PrimitiveCentralIdempotentsUsingConlon uses the function IrrConlon to compute the 
##  primitive central idempotents of QG associated to monomial representations
##  The result is the same as the function PrimitiveCentralIdempotentsBySP but it is slower
##
InstallGlobalFunction(PrimitiveCentralIdempotentsUsingConlon,
function(QG)
local G, zero, one, IrrG, LIrrG, eGKHs, SeGKHs, i, eGKH, K, H; 

    if not IsSemisimpleRationalGroupAlgebra(QG) then
        Error("Wedderga: The input must be a rational group algebra \n");
    fi;
    
    G := UnderlyingMagma( QG );
    zero := Zero( QG );
    one := One( QG );
    
    IrrG := IrrConlon( G );
    LIrrG := Length(IrrG);
    
    eGKHs:=[ IdempotentBySubgroups(QG,G,G) ];
    SeGKHs:=eGKHs[1];
    i:=2;
    while i<=LIrrG do
        #
        # ??? Why sometimes here the component inducedFrom appears ???
        #
        K:=TestMonomial(IrrG[i]).subgroup;
        H:=TestMonomial(IrrG[i]).kernel;
        eGKH:=PrimitiveCentralIdempotentBySP( QG, K, H );    
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
        Print("Warning!!! This is not a complete set of primitive central idempotents !!!\n");
        return eGKHs;
    fi;
end);
  
  
#############################################################################
##
#E
##
