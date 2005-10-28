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

DeclareOperation( "eG",   [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] ); 

DeclareAttribute( "StronglyShodaPairsAndIdempotents", IsGroupRing );

DeclareGlobalFunction( "SearchingKForSSP" );

DeclareGlobalFunction( "PrimitiveCentralIdempotentsByStronglySP" );

#################### idempot.gi #####################

DeclareOperation( "CentralElementBySubgroups",   [ IsGroupRing, IsGroup, IsGroup, IsList, IsList ] ); 

DeclareOperation( "IdempotentBySubgroups", [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList ] );

DeclareOperation( "AverageSum", [ IsGroupRing, IsObject ] );

#################### auxiliar.gi #####################

DeclareOperation("IsCompleteSetOfPCIs",[IsFreeMagmaRing,IsList]);

DeclareOperation( "IsStronglyShodaPair", [ IsGroup, IsGroup, IsGroup ] );

DeclareOperation( "CyclotomicClasses", [ IsPosInt, IsPosInt ] );
DeclareOperation( "BigPrimitiveRoot", [ IsPosInt ] );
DeclareOperation( "BigTrace", [ IsPosInt, IsField, IsObject ] ); 

#################### others.gi #####################

DeclareGlobalFunction( "ListOfPrimitiveCentralIdempotentsBySP" );
DeclareOperation( "PrimitiveCentralIdempotentBySP", [IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] );

DeclareOperation( "IsShodaPair", [ IsGroup, IsGroup, IsGroup ]);

DeclareGlobalFunction( "PrimitiveCentralIdempotentsUsingConlon" );
DeclareGlobalFunction( "PrimitiveCentralIdempotentsByCharacterTable" );


#############################################################################
##
#E
##
