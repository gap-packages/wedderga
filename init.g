#############################################################################
##
#W  init.g                The Wedderga package                Aurora Olivieri
#W                                                              Angel del Rio
#W                                                        Alexander Konovalov
##
#H  $Id$
##
#############################################################################

# announce the package version 
DeclarePackage( "wedderga","4.0", ReturnTrue );

# install the documentation
#DeclarePackageAutoDocumentation( "aclib", "doc" );

# require other packages???????

# read .gd files
ReadPkg( "wedderga", "gap/wedderga.gd" );
