LoadPackage( "guava" : OnlyNeeded ); # to avoid loading SONATA
LoadPackage( "wedderga" );

TestDirectory(DirectoriesPackageLibrary( "wedderga", "tst" ),
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error

