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


DeclareProperty("IsRationalGroupAlgebra", IsGroupRing);

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
DeclareGlobalFunction( "eGKHFromSSP" );
#DeclareOperation( "VerifyShoda", [IsGroup, IsGroup, IsGroup ] );
DeclareOperation( "PCIFromSP", [IsFreeMagmaRing,IsGroup, IsGroup ] );

#viejas
#DeclareGlobalFunction( "PCIsFromSSP1" );
#DeclareGlobalFunction( "PCIsFromSSP2" );
#DeclareGlobalFunction( "PairsOfSubgroupsforPCIs" );