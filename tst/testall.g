LoadPackage( "guava" : OnlyNeeded ); # to avoid loading SONATA
LoadPackage( "irredsol" ); # to avoid variations in tests with no packages
LoadPackage( "wedderga" );

TestDirectory(DirectoriesPackageLibrary( "wedderga", "tst" ),
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error

