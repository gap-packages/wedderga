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
Date := "21/09/2004",
ArchiveURL := "http://.../wedderga-4.0",
ArchiveFormats := ".zoo .tar.gz .tar.bz2 -win.zip",

#TextFiles := ["init.g", ......],
#BinaryFiles := ["doc/manual.dvi", ......],

Persons := [
     rec(
       LastName      := "Olivieri",
       FirstNames    := "Aurora",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "olivieri@usb.ve",
       WWWHome       := "",
       PostalAddress := Concatenation( [
                        "Departamento de Matematicas\n",
                        "Universidad Simon Bolivar\n",
                        "Apartado Postal 89000\n", 
                        "Caracas 1080-A, Venezuela" ] ),
       Place         := "Caracas",
       Institution   := "Universidad Simon Bolivar"
     ),     
     rec(
       LastName      := "del Rio",
       FirstNames    := "Angel",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "adelrio@um.es",
       WWWHome       := "",
       PostalAddress := Concatenation( [
                        "Departamento de Matematicas\n",
                        "Universidad de Murcia\n", 
                        "30100 Murcia, Spain" ] ),
       Place         := "Murcia",
       Institution   := "Universidad de Murcia"
     ),
     rec(
       LastName      := "Konovalov",
       FirstNames    := "Alexander",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "konovalov@member.ams.org",
       WWWHome       := "http://ukrgap.exponenta.ru/konoval.htm",
       PostalAddress := Concatenation( [
                        "P.O.Box 1317\n",
                        "Central Post Office\n", 
                        "Zaporozhye, 69000 Ukraine" ] ),
       Place         := "Zaporozhye",
       Institution   := "Zaporozhye State University"
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
  NeededOtherPackages := [],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

AvailabilityTest := ReturnTrue,
Autoload := false,
#TestFile := "tst/testall.g",

Keywords := ["Wedderburn decomposition of rational group algebrasH"]

));