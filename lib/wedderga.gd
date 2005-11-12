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

DeclareAttribute( "WeddDecomp", IsGroupRing );
DeclareOperation( "WedderburnDecomposition", [IsGroupRing] );

DeclareAttribute( "WeddDecompInfo", IsGroupRing );
DeclareOperation( "WedderburnDecompositionInfo", [ IsGroupRing ]  );

DeclareOperation( "SimpleAlgebraByStronglySP", [ IsSemisimpleFiniteGroupAlgebra, 
                                       IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStronglySPNC", [ IsSemisimpleFiniteGroupAlgebra,
                                       IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStronglySP", [ IsSemisimpleFiniteGroupAlgebra, 
                                       IsGroup, IsGroup, IsPosInt ] );
DeclareOperation( "SimpleAlgebraByStronglySPNC", [ IsSemisimpleFiniteGroupAlgebra,
                                       IsGroup, IsGroup, IsPosInt ] );
DeclareOperation( "SimpleAlgebraByStronglySP", [ IsSemisimpleRationalGroupAlgebra, 
                                       IsGroup, IsGroup] );
DeclareOperation( "SimpleAlgebraByStronglySPNC", [ IsSemisimpleRationalGroupAlgebra,
                                       IsGroup, IsGroup ] ); 

DeclareOperation( "SimpleAlgebraByStronglySPInfo", [ IsSemisimpleRationalGroupAlgebra, 
                                               IsGroup, IsGroup ] );
DeclareOperation( "SimpleAlgebraByStronglySPInfoNC", [ IsSemisimpleRationalGroupAlgebra, 
                                            IsGroup, IsGroup ] );

DeclareOperation( "SimpleAlgebraByStronglySPInfo", [ IsSemisimpleFiniteGroupAlgebra, 
                                            IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStronglySPInfo", [ IsSemisimpleFiniteGroupAlgebra, 
                                            IsGroup, IsGroup, IsPosInt ] );
DeclareOperation( "SimpleAlgebraByStronglySPInfoNC", [ IsSemisimpleFiniteGroupAlgebra, 
                                            IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStronglySPInfoNC", [ IsSemisimpleFiniteGroupAlgebra, 
                                            IsGroup, IsGroup, IsPosInt ] );

DeclareAttribute( "StronglyShodaPairs", IsGroup and IsFinite );

DeclareAttribute( "StronglyShodaPairsAndIdempotents", IsGroupRing );
DeclareGlobalFunction( "SearchingKForSSP" );
DeclareOperation( "eG",   [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] ); 

DeclareGlobalFunction( "PrimitiveCentralIdempotentsByStronglySP" );

#################### idempot.gi #####################

DeclareOperation( "CentralElementBySubgroups",   [ IsGroupRing, IsGroup, IsGroup, IsList, IsList ] ); 

DeclareOperation( "IdempotentBySubgroups", [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList ] );

DeclareOperation( "AverageSum", [ IsGroupRing, IsObject ] );

#################### auxiliar.gi #####################

DeclareOperation("IsCompleteSetOfPCIs",[IsRing,IsList]);

DeclareOperation( "IsStronglyShodaPair", [ IsGroup, IsGroup, IsGroup ] );

DeclareOperation( "CyclotomicClasses", [ IsPosInt, IsPosInt ] );
DeclareOperation( "BigPrimitiveRoot", [ IsPosInt ] );
DeclareOperation( "BigTrace", [ IsPosInt, IsField, IsObject ] ); 

DeclareProperty( "IsStronglyMonomial", IsGroup );

#################### others.gi #####################

DeclareGlobalFunction( "PrimitiveCentralIdempotentsBySP" );
DeclareOperation( "PrimitiveCentralIdempotentBySP", [IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] );

DeclareOperation( "IsShodaPair", [ IsGroup, IsGroup, IsGroup ]);

DeclareGlobalFunction( "PrimitiveCentralIdempotentsUsingConlon" );
DeclareGlobalFunction( "PrimitiveCentralIdempotentsByCharacterTable" );

#############################################################################
##
#E
##
