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
Version        := "4.3",
Date           := "27/12/2007",
ArchiveURL     := "http://www.um.es/adelrio/wedderga/wedderga-4.3",
ArchiveFormats := ".tar.gz .tar.bz2 -win.zip",

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
                        "Departamento de Matemática\n",
                        "Instituto de Ciências Exatas\n",
                        "Universidade Federal de Juiz de Fora\n",
                        "Campus-Cidade Universitária\n",
                        "36036-900, Juiz de Fora, Brazil" ] ),
       Place         := "Sao Paulo",
       Institution   := "Universidade de Sao Paulo"
     ),
     rec(
       LastName      := "Konovalov",
       FirstNames    := "Alexander",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "konovalov@member.ams.org",
       WWWHome       := "http://www.cs.st-andrews.ac.uk/~alexk/",
       PostalAddress := Concatenation( [
                        "School of Computer Science\n",
                        "University of St Andrews\n",
                        "Jack Cole Building, North Haugh,\n",
                        "St Andrews, Fife, KY16 9SX, Scotland" ] ),
       Place         := "St Andrews",
       Institution   := "University of St Andrews"
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
       LastName      := "del Rio",
       FirstNames    := "Angel",
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

Status := "deposited",
# CommunicatedBy := "",
# AcceptDate := "",

README_URL := "http://www.um.es/adelrio/wedderga/README.wedderga",
PackageInfoURL := "http://www.um.es/adelrio/wedderga/PackageInfo.g",
AbstractHTML := "<span class=\"pkgname\">Wedderga</span> is the package to compute the simple components of the Wedderburn decomposition of semisimple group algebras of finite groups over finite fields and over subfields of finite cyclotomic extensions of the rational. It also contains functions that produce the primitive central idempotents of semisimple group algebras. Other functions of <span class=\"pkgname\">Wedderga</span> allows to construct crossed products over a group with coefficients in an associative ring with identity and the multiplication determined by a given action and twisting.",
PackageWWWHome := "http://www.um.es/adelrio/wedderga.htm",
                  
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
  NeededOtherPackages    := [ ["GAPDoc", ">= 1.1"] ],
  SuggestedOtherPackages := [ ["laguna", "3.4"] ],
  ExternalConditions     := []
),

AvailabilityTest := ReturnTrue,
Autoload         := false,
#TestFile        := "tst/testall.g",

Keywords := ["Wedderburn decomposition", "simple components", 
             "central idempotents", "group algebras"]

));