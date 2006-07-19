#############################################################################
##
##  WEDDERGABuildManual()
##
WEDDERGABuildManual:=function()
local wedderga_path, wedderga_main, wedderga_files, wedderga_bookname;
wedderga_path:=Concatenation(
               GAPInfo.PackagesInfo.("wedderga")[1].InstallationPath,"/doc/");
wedderga_main:="manual.xml";
wedderga_files:=["intro.xml", "decomp.xml", "SSP.xml", "idempot.xml",
                 "crossed.xml", "auxiliar.xml", "theory.xml" ];
wedderga_bookname:="wedderga";
MakeGAPDocDoc(wedderga_path, wedderga_main, wedderga_files, wedderga_bookname);  
end;


#############################################################################
##
##  WEDDERGABuildManualHTML()
##
WEDDERGABuildManualHTML:=function()
local wedderga_path, wedderga_main, wedderga_files, str, r, h;
wedderga_path:=Concatenation(
               GAPInfo.PackagesInfo.("wedderga")[1].InstallationPath,"/doc/");
wedderga_main:="manual.xml";
wedderga_files:=["intro.xml", "decomp.xml", "SSP.xml", "idempot.xml",
                 "crossed.xml", "auxiliar.xml", "theory.xml" ];
str:=ComposedXMLString(wedderga_path, wedderga_main, wedderga_files);
r:=ParseTreeXMLString(str);
CheckAndCleanGapDocTree(r);
h:=GAPDoc2HTML(r, wedderga_path);
GAPDoc2HTMLPrintHTMLFiles(h, wedderga_path);
end;
