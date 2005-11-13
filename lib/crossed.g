#############################################################################
##
#W  crossed.g              The Wedderga package           Osnel Broche Cristo
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                              Ángel del Río
##
#H  $Id$
##
#############################################################################


CrossedProductBySSP:= function( G, K, H )
local   N,          # Normalizer of H in G
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
        act,        # Action for the crossed product
        coc;        # Twisting for the crossed product
  
N   := Normalizer(G,H);
ind := Index(G,N);
if N=K then
    ok := Index( K, H );
    if ind=1 then # G=N
        Info( InfoPCI, 2, "N_G(H) = K = G, returning CF(", ok, ")");
        return CF(ok);
    else
        Info( InfoPCI, 2, "N_G(H) = K <> G, returning M_", 
              ind, "( CF(", ok, ") )");
        return FullMatrixAlgebra( CF(ok), ind );
    fi;                          
else # if N_G(H) <> K
    Epi := NaturalHomomorphismByNormalSubgroup( N, H ) ;
    NH  := Image(Epi,N);
    KH  := Image(Epi,K);
    repeat
        k  := Random(KH);
        ok := Order(k);
    until ok = Size(KH);
    Potk:= [ k ];
    for i in [ 2 .. ok ] do
        Potk[i] := Potk[i-1]*k; 
    od;
    Epi2:=NaturalHomomorphismByNormalSubgroup( NH, KH ) ;
    NdK:=Image(Epi2,NH);
        
      act := function(a) 
             return MappingByFunction( CF(ok), CF(ok), x -> 
               GaloisCyc(x, Position(Potk,k^PreImagesRepresentative(Epi2,a))));
             end;
               
      coc := function(a,b)
             return E(ok)^Position( Potk,
                                    PreImagesRepresentative(Epi2,a*b)^-1 *
                                    PreImagesRepresentative(Epi2,a) *
                                    PreImagesRepresentative(Epi2,b) );
             end;        
    if ind=1 then
      Info( InfoPCI, 2, "N_G(H) <> K, returning crossed product");
      return CrossedProduct(CF(ok), NdK, act, coc);
    else
      Info( InfoPCI, 2, 
        "N_G(H) <> K, returning matrix algebra over crossed product");
      return FullMatrixAlgebra( CrossedProduct(CF(ok), NdK, act, coc), ind );
    fi;  
fi;      
end;
