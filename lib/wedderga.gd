#############################################################################
##
#W  wedderga.gd           The Wedderga package                Aurora Olivieri
#W                                                              Angel del Rio
#W                                                        Alexander Konovalov
##
#H  $Id$
##
#############################################################################


DeclareGlobalFunction( "PCIsFromSSP" );
DeclareOperation("IsCompleteSetOfPCIs",[IsFreeMagmaRing,IsList]);

DeclareGlobalFunction( "StronglyShodaPairs" );
DeclareGlobalFunction( "eGKHsFromKHs" );

DeclareGlobalFunction( "SimpleAlgebraFromSSP" );
DeclareGlobalFunction( "SimpleFactorsFromListOfSSP" );

DeclareGlobalFunction( "PCIsUsingConlon" );
DeclareGlobalFunction( "PCIsFromShodaPairs" );
DeclareGlobalFunction( "PCIsUsingCharacterTable" );


DeclareOperation( "Epsilon", [IsFreeMagmaRing,IsGroup, IsGroup ] );
#DeclareOperation( "CentralizerG", [IsFreeMagmaRing,IsElementOfFreeMagmaRing ] );
DeclareOperation( "eGKH", [IsFreeMagmaRing,IsGroup, IsGroup ] );
DeclareOperation( "EpsilonCyclic", [IsFreeMagmaRing,IsGroup, IsGroup ] );
DeclareGlobalFunction( "eGKHFromSSP" );
#DeclareOperation( "VerifyShoda", [IsGroup, IsGroup, IsGroup ] );
DeclareOperation( "PCIFromSP", [IsFreeMagmaRing,IsGroup, IsGroup ] );

#viejas
#DeclareGlobalFunction( "PCIsFromSSP1" );
#DeclareGlobalFunction( "PCIsFromSSP2" );
#DeclareGlobalFunction( "PairsOfSubgroupsforPCIs" );