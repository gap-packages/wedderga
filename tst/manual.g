#############################################################################
##  
#W  manual.g                 The LAGUNA package                  Viktor Bovdi
#W                                                        Alexander Konovalov
#W                                                         Richard Rossmanith
#W                                                            Csaba Schneider
##
#H  $Id$
##
#############################################################################

ReadPackage("GAPDoc","lib/Examples.g");

#############################################################################
##
##  WEDDERGATestManual()
##
WEDDERGATestManual:=function()
local path, TstPkgExamples; 

TstPkgExamples := function ( path, main, files )
  local  str, r, examples, temp_dir, file, otf;
 
  str := ComposedXMLString( path, 
                            Concatenation( main, ".xml" ), 
			    files );
  r := ParseTreeXMLString( str );
  
  examples := Concatenation( 
    "gap> START_TEST( \"Test by GapDoc\" );\n",
    TstExamples( r ),
    "\ngap> STOP_TEST( \"test\", 10000 );\n",
    "Test by GapDoc\nGAP4stones: fail\n" );
  
  temp_dir := DirectoryTemporary( "gapdoc" );
  file := Filename( temp_dir, "testfile" );
  otf := OutputTextFile( file, true );
  SetPrintFormattingStatus( otf, false );
  AppendTo( otf, examples );
  CloseStream( otf );
  
  ReadTest( file );
  
  RemoveFile( file );
  RemoveFile( temp_dir![1] );
  end;

path:=DirectoriesPackageLibrary("wedderga","doc");
Info(InfoPCI, 1, "********** Testing intro.xml **********\n" );   
TstPkgExamples(path,"manual", [ "intro.xml" ] ); 
Info(InfoPCI, 1, "********** Testing decomp.xml **********\n" );   
TstPkgExamples(path,"manual", [ "decomp.xml" ] ); 
Info(InfoPCI, 1, "********** Testing SSP.xml **********\n" );   
TstPkgExamples(path,"manual", [ "SSP.xml" ] );                                                                    
Info(InfoPCI, 1, "********** Testing idempot.xml **********\n" );   
TstPkgExamples(path,"manual", [ "idempot.xml" ] );                                                                     
Info(InfoPCI, 1, "********** Testing crossed.xml **********\n" );   
TstPkgExamples(path,"manual", [ "crossed.xml" ] );                                                                     
Info(InfoPCI, 1, "********** Testing auxiliar.xml **********\n" );   
TstPkgExamples(path,"manual", [ "auxiliar.xml" ] );                                                                    
Info(InfoPCI, 1, "********** Testing theory.xml **********\n" );   
TstPkgExamples(path,"manual", [ "theory.xml" ] );       
end;


#############################################################################
##
#E
##