#############################################################################
##
#W  chairr.gi             The Wedderga package                Aurora Olivieri
#W                                                              Angel del Rio
#W                                                        Alexander Konovalov
##
#H  $Id$
##
#############################################################################


#############################################################################
##
#F IrrSSP( G )
##
InstallGlobalFunction(IrrSSP, 
function(G)
local   caracter, K1, H1, QG, KHs, KH, cara;

    caracter:=function(K1,H1)
    local   Epi, KH, y, CCKH, ind,
            Lista, ChiKH, ChiK, GG;
        if Size(G)=1 then
            return ClassFunction(G,[1]);
        fi;
        
        ind:=Index(K1,H1); 
        if ind=1 then
            if G=K1 then
                return InducedClassFunction(ClassFunction(K1,List(ConjugacyClasses(K1),i->1)),G);
            fi;
        return InducedClassFunction(ClassFunction(K1,List(ConjugacyClasses(K1),i->1)),G);
        fi;
        
        Epi:=NaturalHomomorphismByNormalSubgroup( K1, H1 ) ;
        KH:=Image(Epi,K1);
        y:=Product(IndependentGeneratorsOfAbelianGroup(KH)); 
        CCKH:=ConjugacyClasses(KH);
            
        Lista:=List([1..ind],i->y^i);
        ChiKH:=ClassFunction(KH,List(CCKH,i->E(ind)^Position(Lista,Representative(i))));
        ChiK:=ClassFunction(K1,List(ConjugacyClasses(K1),i->(Representative(i)^Epi)^ChiKH));
        
        #Qchi:=CyclotomicField(ind);
        GG:=GaloisGroup( CF(ind) );
        
    return Orbit(GG , InducedClassFunction(ChiK,G));
    end;

QG:=GroupRing(Rationals, G);
KHs:=StronglyShodaPairs( QG );
cara:=[caracter(G,G)];
KHs:=Difference(KHs,[[G,G]]);
for KH in KHs do
   cara:=Concatenation(cara,caracter(KH[1],KH[2]));
od;
return cara;
end);


#############################################################################
##
##  
##  
##  
##
#F IrrSSP2( G )
##
InstallGlobalFunction(IrrSSP2, 
function(G)
local   caracter, K1, H1, QG, KHs, KH, cara;

    caracter:=function(K1,H1)
    local   Epi, KH, y, CCKH, ind,
            Lista, ChiKH, CCK, ChiK, l, i, x, Orb;
            #Qchi,GG;
        if Size(G)=1 then
            return ClassFunction(G,[1]);
        fi;
        
        ind:=Index(K1,H1); 
        if ind=1 then
            if G=K1 then
                return InducedClassFunction(ClassFunction(K1,List(ConjugacyClasses(K1),i->1)),G);
            fi;
        return InducedClassFunction(ClassFunction(K1,List(ConjugacyClasses(K1),i->1)),G);
        fi;
        
        Epi:=NaturalHomomorphismByNormalSubgroup( K1, H1 ) ;
        KH:=Image(Epi,K1);
        y:=Product(IndependentGeneratorsOfAbelianGroup(KH)); 
        CCKH:=List(ConjugacyClasses(KH),i->Representative(i));
            
        Lista:=List([1..ind],i->y^i);
        l:=PrimeResidues( ind );
        Orb:=[];
        CCK:=List(ConjugacyClasses(K1),i->Representative(i));
        for i in l do
            ChiKH:=ClassFunction(KH,List(CCKH,x->(E(ind)^i)^Position(Lista,x)));
            ChiK:=ClassFunction(K1,List(CCK,x->(x^Epi)^ChiKH));
            if not(InducedClassFunction(ChiK,G) in Orb) then
                Add(Orb,InducedClassFunction(ChiK,G));
            fi;
        od;
        
    return Orb;
    end;

QG:=GroupRing(Rationals, G);
KHs:=StronglyShodaPairs( QG );
cara:=[caracter(G,G)];
KHs:=Difference(KHs,[[G,G]]);
for KH in KHs do
   cara:=Concatenation(cara,caracter(KH[1],KH[2]));
od;
return cara;
end);
