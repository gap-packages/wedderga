#############################################################################
##
#W  wedderga.gd           The Wedderga package            Osnel Broche Cristo
#W                                                         Olexandr Konovalov
#W                                                            Aurora Olivieri
#W                                                           Gabriela Olteanu
#W                                                              Ángel del Río
#W                                                          Inneke Van Gelder
##
#############################################################################


#############################################################################
##
##  InfoWedderga
##  
##  We declare new Info class for Wedderga algorithms. 
##  It has 3 levels - 0, 1 (default) and 2
##  To change Info level to k, use command SetInfoLevel(InfoWedderga, k)
DeclareInfoClass("InfoWedderga");

DeclareProperty( "IsSemisimpleRationalGroupAlgebra", IsGroupRing );
DeclareProperty( "IsSemisimpleFiniteGroupAlgebra", IsGroupRing );
DeclareProperty( "IsSemisimpleZeroCharacteristicGroupAlgebra", IsGroupRing );
DeclareProperty( "IsCFGroupAlgebra", IsGroupRing );
DeclareProperty( "IsSemisimpleANFGroupAlgebra", IsGroupRing );

#################### ExtremeSSPs.gi #####################

DeclareGlobalFunction( "IsMaximalAbelianFactorGroup" );
DeclareGlobalFunction( "SearchingNNKForSSP" );
DeclareAttribute( "ExtSSPAndDim", IsGroup and IsFinite );
DeclareGlobalFunction( "ExtremelyStrongShodaPairs" );
DeclareProperty( "IsNormallyMonomial", IsGroup );
DeclareAttribute( "PrimitiveCentralIdempotentsByExtSSP", IsGroupRing );
DeclareGlobalFunction( "PrimitiveCentralIdempotentsByESSP");
DeclareOperation( "IsExtremelyStrongShodaPair", [ IsGroup, IsGroup, IsGroup ] );

#################### main.gi #####################

DeclareOperation( "WedderburnDecomposition", [IsSemisimpleANFGroupAlgebra] );
DeclareOperation( "WedderburnDecomposition", [IsSemisimpleFiniteGroupAlgebra] );
DeclareAttribute( "WeddDecomp", IsSemisimpleANFGroupAlgebra );
DeclareAttribute( "WedderburnDecompositionInfo", IsSemisimpleANFGroupAlgebra );
DeclareAttribute( "WedderburnDecompositionInfo", IsSemisimpleFiniteGroupAlgebra );

DeclareOperation( "SimpleAlgebraByStrongSP", [ IsGroupRing, 
                                       IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStrongSPNC", [ IsGroupRing,
                                       IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStrongSPInfo", [ IsGroupRing, 
                                            IsGroup, IsGroup, IsList ] );
DeclareOperation( "SimpleAlgebraByStrongSPInfoNC", [ IsGroupRing, 
                                            IsGroup, IsGroup, IsList ] );

DeclareAttribute( "StrongShodaPairs", IsGroup and IsFinite );

DeclareAttribute( "StrongShodaPairsAndIdempotents", IsGroupRing );
DeclareAttribute( "SSPNonESSPAndTheirIdempotents", IsGroupRing );
DeclareOperation( "Idempotent_eGsum",   [ IsSemisimpleRationalGroupAlgebra, IsGroup, IsGroup ] ); 

DeclareAttribute( "PrimitiveCentralIdempotentsByStrongSP", IsGroupRing );

DeclareOperation( "AddCrossedProductBySSP", [ IsGroup, IsGroup, IsGroup]);
DeclareOperation( "AddCrossedProductBySST", 
                  [ IsInt, IsInt, IsField, IsGroup, IsList ] );
DeclareAttribute( "WeddDecompData", IsGroup );
DeclareOperation( "GenWeddDecomp", [ IsGroupRing ] );
DeclareOperation( "SimpleAlgebraByData", [ IsList ] );
DeclareOperation( "SimpleAlgebraInfoByData", [ IsList ] );
DeclareOperation( "SimpleAlgebraByCharacter", [ IsSemisimpleANFGroupAlgebra, IsCharacter ] );
DeclareOperation( "SimpleAlgebraByCharacter", [ IsSemisimpleFiniteGroupAlgebra, IsCharacter]);
DeclareOperation( "SimpleAlgebraByCharacterInfo", [ IsSemisimpleANFGroupAlgebra, IsCharacter ] );
DeclareOperation( "SimpleAlgebraByCharacterInfo", [ IsSemisimpleFiniteGroupAlgebra, IsCharacter ] );

DeclareOperation( "PrimitiveIdempotentsNilpotent", [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList] );
DeclareOperation( "PrimitiveIdempotentsTrivialTwisting", [ IsSemisimpleFiniteGroupAlgebra, IsGroup, IsGroup, IsList, IsList] );

#################### idempot.gi #####################

DeclareOperation( "CentralElementBySubgroups",   
        [ IsGroupRing, IsGroup, IsGroup, IsList, IsList, IsList ] );
        
DeclareOperation( "IdempotentBySubgroups", 
        [ IsGroupRing, IsGroup, IsGroup, IsList, IsList, IsList ] );

DeclareOperation( "AverageSum", [ IsGroupRing, IsObject ] );

#################### auxiliar.gi #####################

DeclareOperation("IsCompleteSetOfPCIs",[IsRing,IsList]);
DeclareOperation("IsCompleteSetOfOrthogonalIdempotents",[IsRing,IsList]);

DeclareOperation( "IsStrongShodaPair", [ IsGroup, IsGroup, IsGroup ] );

DeclareOperation( "CyclotomicClasses", [ IsPosInt, IsPosInt ] );
DeclareOperation( "BigPrimitiveRoot", [ IsPosInt ] );
DeclareOperation( "BigTrace", [ IsPosInt, IsField, IsObject ] ); 
DeclareProperty( "IsStronglyMonomial", IsGroup );

DeclareOperation( "IsCyclotomicClass", [ IsPosInt, IsPosInt, IsList ] );
DeclareAttribute( "IsCyclGroupAlgebra", IsGroupRing );
DeclareOperation( "SizeOfSplittingField", [IsCharacter, IsPosInt] ); 

DeclareOperation( "SquareRootMod", [ IsPosInt, IsPosInt ]); 
DeclareOperation( "SquaresMod", [IsPosInt] );
DeclareOperation( "SolveEquation2@", [IsPosInt] );
DeclareOperation( "SolveEquation3@", [IsPosInt] );
DeclareOperation( "SolveEquation@", [IsField] );
DeclareOperation( "PrimRootOfUnity", [IsField, IsPosInt] );
DeclareOperation( "MakeMatrixByBasis", [ IsMapping, IsBasis] );
DeclareOperation( "ReturnGalElement", [ IsObject , IsGroup, IsGroup, IsGroup, IsField, IsObject] );
DeclareOperation( "LeftMultiplicationBy", [ IsObject , IsField] );
DeclareOperation( "MakeLinearCombination", [IsAlgebra, IsList, IsList] );
DeclareOperation( "Product3Lists", [IsList] );

DeclareOperation( "IsTwistingTrivial", [IsGroup, IsGroup, IsGroup] );

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

DeclareOperation( "CodeWordByGroupRingElement", [ IsField, IsSet, IsObject ] );
DeclareOperation( "CodeByLeftIdeal", [ IsField, IsGroup, IsSet, IsRing ] );

#################### bw.gi #####################

DeclareOperation( "LinCharByKernel", [ IsGroup, IsGroup ] );
DeclareOperation( "BWNoStMon", [IsGroup and IsFinite] );
DeclareOperation( "ReductionModnZ", [IsPosInt, IsPosInt] );
DeclareOperation( "GalToInt", [IsGroup] );
DeclareGlobalFunction( "CocycleByData" );
DeclareOperation( "SimpleAlgebraByStrongSTInfo", 
                            [IsPosInt , IsPosInt , IsField , IsGroup , IsList] );

#############################################################################
##
#E
##
