#############################################################################
##
#W  init.g                The Wedderga package            Osnel Broche Cristo
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                              Ángel del Río
##
#H  $Id$
##
#############################################################################

# read Wedderga declarations
ReadPackage( "wedderga/lib/wedderga.gd" );
ReadPackage( "wedderga/lib/crossed.gd" );

# read the other part of code
ReadPackage("wedderga/lib/wedderga.g");

# set the default InfoLevel
SetInfoLevel( InfoPCI, 1 );
