#############################################################################
##
#W  PackageInfo.g         The Wedderga package                Aurora Olivieri
#W                                                              Angel del Rio
#W                                                        Alexander Konovalov
##
#H  $Id$
##
#############################################################################

SetPackageInfo( rec(

PackageName := "Wedderga",
Subtitle := "Wedderga",
Version := "4.0",
Date := "04/04/2004",
ArchiveURL := "http://.../wedderga-4.0",
ArchiveFormats := ".zoo .tar.gz .tar.bz2 -win.zip",

#TextFiles := ["init.g", ......],
#BinaryFiles := ["doc/manual.dvi", ......],

Persons := [
     rec(
     LastName := "Olivieri",
     FirstNames := "Aurora",
     IsAuthor := true,
     IsMaintainer := true,
     Email := "",
     WWWHome := "",
     PostalAddress := "",
     Institution := ""
     ),     
     rec(
     LastName := "del Rio",
     FirstNames := "Angel",
     IsAuthor := true,
     IsMaintainer := true,
     Email := "",
     PostalAddress := "",
     ),
     rec(
     LastName := "Konovalov",
     FirstNames := "Alexander",
     IsAuthor := true,
     IsMaintainer := true,
     Email := "konovalov@member.ams.org",
     WWWHome := "http://ukrgap.exponenta.ru/konoval.htm",
     PostalAddress := "P.O.Box 1317, Central Post Office, Zaporozhye, 69000 Ukraine",
     Institution := "Department of Mathematics, Zaporozhye State University"
     )
],

Status := "",
CommunicatedBy := "",
AcceptDate := "",

README_URL := "http://.../README.wedderga",
PackageInfoURL := "http://.../PackageInfo.g",
AbstractHTML := "The <span class=\"pkgname\">Wedderga</span> package ...",
PackageWWWHome := "http://.../wedderga.htm",
                  
PackageDoc := rec(
  BookName := "Wedderga",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile := "doc/manual.pdf",
  SixFile := "doc/manual.six",
  LongTitle := "Wedderga",
  Autoload := false
),


Dependencies := rec(
  GAP := ">=4.4",
  NeededOtherPackages := [["GAPDoc", ">= 0.999"]],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

AvailabilityTest := ReturnTrue,
Autoload := false,
#TestFile := "tst/testall.g",

Keywords := ["Wedderburn decomposition"]

));