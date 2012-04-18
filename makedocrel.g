###########################################################################
##
#W makedocrel.g            The Wedderga package         Alexander Konovalov
##
###########################################################################

ExtractManualExamples:=function( pkgname, main, files )
local path, tst, i, s, name, output;
path:=Directory( 
        Concatenation(PackageInfo(pkgname)[1].InstallationPath, "/doc") );
Print("===============================================================\n");
Print("Extracting manual examples for ", pkgname, " package\n" );
Print("===============================================================\n");
tst:=ManualExamples( path, main, files, "Chapter" );
for i in [ 1 .. Length(tst) ] do 
  Print( "Processing '", pkgname, 
         "' chapter number ", i, " of ", Length(tst), "\c" );
  if Length( tst[i] ) > 0 then
    s := String(i);
    if Length(s)=1 then 
      # works for <100 chapters
      s:=Concatenation("0",s); 
    fi;
    name := Filename( 
              Directory( 
                Concatenation( PackageInfo(pkgname)[1].InstallationPath, 
                               "/tst" ) ), 
                Concatenation( pkgname, s, ".tst" ) );
    output := OutputTextFile( name, false ); # to empty the file first
    SetPrintFormattingStatus( output, false ); # to avoid line breaks
    PrintTo( output, tst[i] );
    CloseStream(output);
    # one superfluous check
    if tst[i] <> StringFile( name ) then
      Error("Saved file does not match original examples string!!!\n");  
    else
      Print(" - OK! \n" );
    fi;
  else
    Print(" - no examples to save! \n" );    
  fi;  
od;
Print("===============================================================\n");
end;

###########################################################################

WEDDERGAMANUALFILES:=[ 
"../PackageInfo.g"
];

MakeGAPDocDoc( "doc",     # path to the directory containing the main file
               "manual",  # the name of the main file (without extension)
                          # list of (probably source code) files relative 
                          # to path which contain pieces of documentation 
                          # which must be included in the document
               WEDDERGAMANUALFILES, 
               "Wedderga",# the name of the book used by GAP's online help
               "../../..",# optional: relative path to the main GAP root 
                          # directory to produce HTML files with relative 
                          # paths to external books.
               "MathJax"  # optional: use "MathJax", "Tth" and/or "MathML"
                          # to produce additional variants of HTML files
               );; 

# Copy the *.css and *.js files from the styles directory of the GAPDoc 
# package into the directory containing the package manual.
CopyHTMLStyleFiles( "doc" );

# Create the manual.lab file which is needed if the main manuals or another 
# package is referring to your package
GAPDocManualLab( "Wedderga" );; 
 
ExtractManualExamples( "wedderga", "manual.xml", WEDDERGAMANUALFILES );

QUIT;
###########################################################################
##
#E
##