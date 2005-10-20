#############################################################################
##
#W  wedderga.gd           The Wedderga package                Aurora Olivieri
#W                                                              Angel del Rio
#W                                                        Alexander Konovalov
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

#################### main.gi #####################

DeclareAttribute( "WedderburnDecomposition", IsSemisimpleFiniteGroupAlgebra );
DeclareAttribute( "WedderburnDecompositionInfo", IsGroupRing  );

DeclareOperation( "SimpleAlgebra", [ IsSemisimpleFiniteGroupAlgebra, 
                                       IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraNC", [ IsSemisimpleFiniteGroupAlgebra,
                                       IsGroup, IsGroup, IsList ] ); 

DeclareOperation( "SimpleAlgebraInfo", [ IsSemisimpleRationalGroupAlgebra, 
                                               IsGroup, IsGroup ] );
DeclareOperation( "SimpleAlgebraInfoNC", [ IsSemisimpleRationalGroupAlgebra, 
                                            IsGroup, IsGroup ] );


DeclareOperation( "SimpleAlgebraInfo", [ IsSemisimpleFiniteGroupAlgebra, 
                                            IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraInfoNC", [ IsSemisimpleFiniteGroupAlgebra, 
                                            IsGroup, IsGroup, IsList ] );


DeclareAttribute( "StronglyShodaPairs", IsGroup and IsFinite );

DeclareAttribute( "StronglyShodaPairsAndIdempotents", IsSemisimpleRationalGroupAlgebra );
DeclareAttribute( "StronglyShodaPairsAndIdempotents", IsSemisimpleFiniteGroupAlgebra );

DeclareGlobalFunction( "SearchingKForSSP" );
DeclareGlobalFunction( "eG" );

DeclareGlobalFunction( "PrimitiveCentralIdempotentsFromSSP" );

#################### idempot.gi #####################

DeclareOperation( "eGKH",  [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] );
DeclareOperation( "eGKHc", [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList ] ); 
DeclareOperation( "eGKHc", [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList ] );

DeclareOperation( "Epsilon", [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList ] );

DeclareOperation( "Hat", [ IsGroupRing, IsObject ] );

#################### auxiliar.gi #####################

DeclareOperation("IsCompleteSetOfPCIs",[IsFreeMagmaRing,IsList]);

DeclareOperation( "IsStronglyShodaPair", [ IsGroup, IsGroup, IsGroup ] );

DeclareOperation( "CentralizerG", [ IsFreeMagmaRing, IsElementOfFreeMagmaRing ] );
DeclareOperation( "Conjugate", [ IsFreeMagmaRing, IsElementOfFreeMagmaRing, IsObject ] ); 

DeclareOperation( "CyclotomicClasses", [ IsPosInt, IsPosInt ] );
DeclareOperation( "BigPrimitiveRoot", [ IsPosInt ] );
DeclareOperation( "BigTrace", [ IsPosInt, IsField, IsObject ] ); 

#################### others.gi #####################

DeclareGlobalFunction( "PCIsFromShodaPairs" );
DeclareOperation( "PCIFromSP", [IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] );

DeclareOperation( "IsShodaPair", [ IsGroup, IsGroup, IsGroup ]);

DeclareGlobalFunction( "PCIsUsingConlon" );
DeclareGlobalFunction( "PCIsUsingCharacterTable" );