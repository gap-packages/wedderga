#############################################################################
##
#W  wedderga.gd           The Wedderga package            Osnel Broche Cristo
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                              Ángel del Río
##
#H  $Id$
##
#############################################################################


#############################################################################
##
##  InfoPCI
##  
##  We declare new Info class for Wedderga algorithms. 
##  It has 3 levels - 0, 1 (default) and 2
##  To change Info level to k, use command SetInfoLevel(InfoPCI, k)
DeclareInfoClass("InfoPCI");

DeclareProperty( "IsSemisimpleRationalGroupAlgebra", IsGroupRing );
DeclareProperty( "IsSemisimpleFiniteGroupAlgebra", IsGroupRing );
DeclareProperty( "IsZeroCharacteristicGroupAlgebra", IsGroupRing );

#################### main.gi #####################

DeclareAttribute( "WeddDecomp", IsGroupRing );
DeclareOperation( "WedderburnDecomposition", [IsGroupRing] );

DeclareAttribute( "WeddDecompInfo", IsGroupRing );
DeclareOperation( "WedderburnDecompositionInfo", [ IsGroupRing ]  );

DeclareOperation( "SimpleAlgebraByStronglySP", [ IsGroupRing, 
                                       IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStronglySPNC", [ IsGroupRing,
                                       IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStronglySPInfo", [ IsGroupRing, 
                                            IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStronglySPInfoNC", [ IsGroupRing, 
                                            IsGroup, IsGroup, IsList ] );

DeclareAttribute( "StronglyShodaPairs", IsGroup and IsFinite );

DeclareAttribute( "StronglyShodaPairsAndIdempotents", IsGroupRing );
DeclareGlobalFunction( "SearchingKForSSP" );
DeclareOperation( "eG",   [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] ); 

DeclareGlobalFunction( "PrimitiveCentralIdempotentsByStronglySP" );

DeclareOperation( "AddCrossedProductBySSP", [ IsGroup, IsGroup, IsGroup]);
DeclareOperation( "AddCrossedProductBySST", 
                  [ IsInt, IsInt, IsCyclotomicField, IsGroup, IsList ] );
DeclareAttribute( "WeddDecompData", IsGroup );
DeclareOperation( "GenWeddDecomp", [ IsGroupAlgebra ] );
DeclareOperation( "SimpleAlgebraByData", [ IsList ] );
DeclareOperation( "WDecomp", [ IsGroupAlgebra ] );
DeclareOperation( "SimpleAlgebraInfoByData", [ IsList ] );
DeclareOperation( "WDecompInfo", [ IsGroupAlgebra ] );

#################### idempot.gi #####################

DeclareOperation( "CentralElementBySubgroups",   
        [ IsGroupRing , IsGroup, IsGroup, IsList, IsList ] ); 
        
DeclareOperation( "IdempotentBySubgroups", 
        [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList ] );
DeclareOperation( "IdempotentBySubgroups", 
            [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ] );
DeclareOperation( "IdempotentBySubgroups", 
            [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsPosInt ] );


DeclareOperation( "AverageSum", [ IsGroupRing, IsObject ] );

#################### auxiliar.gi #####################

DeclareOperation("IsCompleteSetOfPCIs",[IsRing,IsList]);

DeclareOperation( "IsStronglyShodaPair", [ IsGroup, IsGroup, IsGroup ] );

DeclareOperation( "CyclotomicClasses", [ IsPosInt, IsPosInt ] );
DeclareOperation( "BigPrimitiveRoot", [ IsPosInt ] );
DeclareOperation( "BigTrace", [ IsPosInt, IsField, IsObject ] ); 
DeclareProperty( "IsStronglyMonomial", IsGroup );

DeclareOperation( "IsCyclotomicClass", [ IsPosInt, IsPosInt, IsList ] );

#################### others.gi #####################

DeclareAttribute( "ShodaPairsAndIdempotents", IsSemisimpleRationalGroupAlgebra );
DeclareGlobalFunction( "PrimitiveCentralIdempotentsBySP" );
DeclareOperation( "PrimitiveCentralIdempotentBySP", 
                        [IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] );
DeclareOperation( "PCIBySP", 
                        [IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] );

DeclareOperation( "IsShodaPair", [ IsGroup, IsGroup, IsGroup ]);
DeclareOperation( "IsSP", [ IsGroup, IsGroup, IsGroup ]);


DeclareGlobalFunction( "PrimitiveCentralIdempotentsUsingConlon" );
DeclareGlobalFunction( "PrimitiveCentralIdempotentsByCharacterTable" );




#################### bw.gi #####################


DeclareOperation( "LinCharByStronglySP", [ IsGroup, IsGroup ] );
DeclareOperation( "BW", [IsGroup and IsFinite] );
DeclareOperation( "BWNoStMon", [IsGroup and IsFinite] );
DeclareOperation( "red", [IsPosInt, IsPosInt] );
DeclareOperation( "GalToInt", [IsGroup] );
DeclareGlobalFunction( "cocycle" );
DeclareOperation( "SimpleAlgebraByStronglySTInfo", 
                            [IsPosInt , IsPosInt , IsField , IsGroup , IsList] );


#############################################################################
##
#E
##
