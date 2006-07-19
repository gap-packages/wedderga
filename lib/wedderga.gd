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
DeclareProperty( "IsSemisimpleZeroCharacteristicGroupAlgebra", IsGroupRing );
DeclareProperty( "IsCFGroupAlgebra", IsGroupRing );
DeclareProperty( "IsSemisimpleANFGroupAlgebra", IsGroupRing );

#################### main.gi #####################

DeclareOperation( "WedderburnDecomposition", [IsSemisimpleANFGroupAlgebra] );
DeclareOperation( "WedderburnDecomposition", [IsSemisimpleFiniteGroupAlgebra] );
DeclareAttribute( "WeddDecomp", IsSemisimpleANFGroupAlgebra );
DeclareAttribute( "WedderburnDecompositionInfo", IsSemisimpleANFGroupAlgebra );
DeclareAttribute( "WedderburnDecompositionInfo", IsSemisimpleFiniteGroupAlgebra );

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
DeclareOperation( "eGsum",   [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] ); 

DeclareGlobalFunction( "PrimitiveCentralIdempotentsByStronglySP" );

DeclareOperation( "AddCrossedProductBySSP", [ IsGroup, IsGroup, IsGroup]);
DeclareOperation( "AddCrossedProductBySST", 
                  [ IsInt, IsInt, IsField, IsGroup, IsList ] );
DeclareAttribute( "WeddDecompData", IsGroup );
DeclareOperation( "GenWeddDecomp", [ IsGroupAlgebra ] );
DeclareOperation( "SimpleAlgebraByData", [ IsList ] );
DeclareOperation( "SimpleAlgebraInfoByData", [ IsList ] );
DeclareOperation( "SimpleAlgebraByCharacter", [ IsSemisimpleANFGroupAlgebra, IsCharacter ] );
DeclareOperation( "SimpleAlgebraByCharacter", [ IsSemisimpleFiniteGroupAlgebra, IsCharacter]);
DeclareOperation( "SimpleAlgebraByCharacterInfo", [ IsSemisimpleANFGroupAlgebra, IsCharacter ] );
DeclareOperation( "SimpleAlgebraByCharacterInfo", [ IsSemisimpleFiniteGroupAlgebra, IsCharacter ] );

#################### idempot.gi #####################

DeclareOperation( "CentralElementBySubgroups",   
        [ IsGroupRing, IsGroup, IsGroup, IsList, IsList, IsList ] );
        
DeclareOperation( "IdempotentBySubgroups", 
        [ IsGroupRing, IsGroup, IsGroup, IsList, IsList, IsList ] );

DeclareOperation( "AverageSum", [ IsGroupRing, IsObject ] );

#################### auxiliar.gi #####################

DeclareOperation("IsCompleteSetOfPCIs",[IsRing,IsList]);
DeclareOperation("IsCompleteSetOfOrthIdemps",[IsRing,IsList]);

DeclareOperation( "IsStronglyShodaPair", [ IsGroup, IsGroup, IsGroup ] );

DeclareOperation( "CyclotomicClasses", [ IsPosInt, IsPosInt ] );
DeclareOperation( "BigPrimitiveRoot", [ IsPosInt ] );
DeclareOperation( "BigTrace", [ IsPosInt, IsField, IsObject ] ); 
DeclareProperty( "IsStronglyMonomial", IsGroup );

DeclareOperation( "IsCyclotomicClass", [ IsPosInt, IsPosInt, IsList ] );
DeclareAttribute( "IsCyclGroupAlgebra", IsGroupRing );
DeclareOperation( "SizeOfSplittingField", [IsCharacter, IsPosInt] ); 

#################### others.gi #####################

DeclareAttribute( "ShodaPairsAndIdempotents", IsGroupRing );
DeclareGlobalFunction( "PrimitiveCentralIdempotentsBySP" );
DeclareOperation( "PrimitiveCentralIdempotentBySP", 
                        [IsGroupRing, IsGroup, IsGroup ] );

DeclareOperation( "IsShodaPair", [ IsGroup, IsGroup, IsGroup ]);

DeclareGlobalFunction( "PrimitiveCentralIdempotentsUsingConlon" );
#DeclareGlobalFunction( "PrimitiveCentralIdempotentsByCharacterTable" );
DeclareOperation( "PrimitiveCentralIdempotentsByCharacterTable", [ IsSemisimpleANFGroupAlgebra ] );
DeclareOperation( "PrimitiveCentralIdempotentsByCharacterTable", [ IsSemisimpleFiniteGroupAlgebra ] );

#################### bw.gi #####################

DeclareOperation( "LinCharByKernel", [ IsGroup, IsGroup ] );
DeclareOperation( "BWNoStMon", [IsGroup and IsFinite] );
DeclareOperation( "ReductionModnZ", [IsPosInt, IsPosInt] );
DeclareOperation( "GalToInt", [IsGroup] );
DeclareGlobalFunction( "CocycleByData" );
DeclareOperation( "SimpleAlgebraByStronglySTInfo", 
                            [IsPosInt , IsPosInt , IsField , IsGroup , IsList] );

#############################################################################
##
#E
##
