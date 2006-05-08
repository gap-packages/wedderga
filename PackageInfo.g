#############################################################################
##
#W  PackageInfo.g         The Wedderga package            Osnel Broche Cristo
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                              Ángel del Río
##
#H  $Id$
##
#############################################################################

SetPackageInfo( rec(

PackageName    := "Wedderga",
Subtitle       := Concatenation( [
                  "Central idempotents and simple components",
                  "of group algebras" ] ),
Version        := "4.0",
Date           := "08/05/2006",
ArchiveURL     := "http://.../wedderga-4.0",
ArchiveFormats := ".zoo .tar.gz .tar.bz2 -win.zip",

#TextFiles     := ["init.g", ......],
#BinaryFiles   := ["doc/manual.dvi", ......],

Persons :=
 [
     rec(
       LastName      := "Broche Cristo",
       FirstNames    := "Osnel",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "osnelier@ime.usp.br",
       PostalAddress := Concatenation( [
                        "Instituto de Matemática\n",
                        "Universidade de Sao Paulo\n",
                        "Caixa Postal 66281, Sao Paulo\n", 
                        "CEP 05315-970, Brazil" ] ),
       Place         := "Sao Paulo",
       Institution   := "Universidade de Sao Paulo"
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
       Institution   := "Zaporozhye National University"
     ),
    rec(
       LastName      := "Olteanu",
       FirstNames    := "Gabriela",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "golteanu@um.es",
       PostalAddress := Concatenation( [
                        "Departamento de Matemáticas\n",
                        "Universidad de Murcia\n", 
                        "30100 Murcia, Spain" ] ),
       Place         := "Murcia",
       Institution   := "Universidad de Murcia"
     ),
     rec(
       LastName      := "Olivieri",
       FirstNames    := "Aurora",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "olivieri@usb.ve",
       PostalAddress := Concatenation( [
                        "Departamento de Matemáticas\n",
                        "Universidad Simón Bolívar\n",
                        "Apartado Postal 89000\n", 
                        "Caracas 1080-A, Venezuela" ] ),
       Place         := "Caracas",
       Institution   := "Universidad Simón Bolívar"
     ),     
     rec(
       LastName      := "del Río",
       FirstNames    := "Ángel",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "adelrio@um.es",
       WWWHome       := "http://www.um.es/adelrio",
       PostalAddress := Concatenation( [
                        "Departamento de Matemáticas\n",
                        "Universidad de Murcia\n", 
                        "30100 Murcia, Spain" ] ),
       Place         := "Murcia",
       Institution   := "Universidad de Murcia"
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
  GAP                    := ">=4.4",
  NeededOtherPackages    := [],
  SuggestedOtherPackages := [["laguna", "3.3"]],
  ExternalConditions     := []
),

AvailabilityTest := ReturnTrue,
Autoload         := false,
#TestFile        := "tst/testall.g",

Keywords := ["Wedderburn decomposition", "simple components", 
             "central idempotents", "group algebras"]

));
