#############################################################################
##
#W  init.g                The Wedderga package            Osnel Broche Cristo
#W                                                         Olexandr Konovalov
#W                                                            Aurora Olivieri
#W                                                           Gabriela Olteanu
#W                                                              Ángel del Río
#W                                                          Inneke Van Gelder
##
#############################################################################

#I introducing globally the NC versions of PreImages...  
if not IsBound( PreImagesNC ) then 
    BindGlobal( "PreImagesNC", PreImages ); 
fi; 
if not IsBound( PreImagesElmNC ) then 
    BindGlobal( "PreImagesElmNC", PreImagesElm ); 
fi; 
if not IsBound( PreImagesRepresentativeNC ) then 
    BindGlobal( "PreImagesRepresentativeNC", PreImagesRepresentative ); 
fi; 

# read Wedderga declarations
ReadPackage( "wedderga", "lib/wedderga.gd" );
ReadPackage( "wedderga", "lib/crossed.gd" );
ReadPackage( "wedderga", "lib/div-alg.gd" );

# set the default InfoLevel
SetInfoLevel( InfoWedderga, 1 );
