#############################################################################
##
#W  PackageInfo.g         The Wedderga package            Osnel Broche Cristo
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                           Gabriela Olteanu
#W                                                              �ngel del R�o
##
#############################################################################

SetPackageInfo( rec(

PackageName    := "Wedderga",
Subtitle       := Concatenation( [
                  "Wedderburn Decomposition of Group Algebras" ] ),
Version        := "4.5.1",
Date           := "31/05/2012",
##  <#GAPDoc Label="PKGVERSIONDATA">
##  <!ENTITY VERSION "4.4.3">
##  <!ENTITY RELEASEDATE "31 May 2012">
##  <!ENTITY RELEASEYEAR "2012">
##  <#/GAPDoc>

PackageWWWHome := "http://www.cs.st-andrews.ac.uk/~alexk/wedderga/",

ArchiveURL := Concatenation( ~.PackageWWWHome, "wedderga-", ~.Version ),
ArchiveFormats := ".tar.gz",

Persons :=
 [
     rec(
       LastName      := "Broche Cristo",
       FirstNames    := "Osnel",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "osnel@ufla.br",
       PostalAddress := Concatenation( [
                        "Departamento de Ci�ncias Exatas\n",
                        "Universidade Federal de Lavras - UFLA\n",
                        "Campus Universit�rio - Caixa Postal 3037\n",
                        "37200-000, Lavras - MG, Brazil" ] ),
       Place         := "Lavras - MG",
       Institution   := "Universidade Federal de Lavras - UFLA"
     ),
     rec(
       LastName      := "Konovalov",
       FirstNames    := "Alexander",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "alexk@mcs.st-andrews.ac.uk",
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
       Email         := "golteanu@um.es, olteanu@math.ubbcluj.ro",
       PostalAddress := Concatenation( [
                        "Department of Mathematics and Computer Science\n",
                        "North University of Baia Mare\n",
                        "Victoriei 76, 430122 Baia Mare, Romania" ] ),
       Place         := "Baia Mare",
       Institution   := "North University of Baia Mare"
     ),
     rec(
       LastName      := "Olivieri",
       FirstNames    := "Aurora",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "olivieri@usb.ve",
       PostalAddress := Concatenation( [
                        "Departamento de Matem�ticas\n",
                        "Universidad Sim�n Bol�var\n",
                        "Apartado Postal 89000\n", 
                        "Caracas 1080-A, Venezuela" ] ),
       Place         := "Caracas",
       Institution   := "Universidad Sim�n Bol�var"
     ),     
     rec(
       LastName      := "del Rio",
       FirstNames    := "Angel",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "adelrio@um.es",
       WWWHome       := "http://www.um.es/adelrio",
       PostalAddress := Concatenation( [
                        "Departamento de Matem�ticas\n",
                        "Universidad de Murcia\n", 
                        "30100 Murcia, Spain" ] ),
       Place         := "Murcia",
       Institution   := "Universidad de Murcia"
     )
],

Status := "accepted",
CommunicatedBy := "Gerhard Hiss (Aachen)",
AcceptDate := "01/2008",

README_URL := 
  Concatenation( ~.PackageWWWHome, "README" ),
PackageInfoURL := 
  Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
  
AbstractHTML := "<span class=\"pkgname\">Wedderga</span> is the package to compute the simple components of the Wedderburn decomposition of semisimple group algebras of finite groups over finite fields and over subfields of finite cyclotomic extensions of the rational. It also contains functions that produce the primitive central idempotents of semisimple group algebras. Other functions of <span class=\"pkgname\">Wedderga</span> allows to construct crossed products over a group with coefficients in an associative ring with identity and the multiplication determined by a given action and twisting.",
                  
PackageDoc := rec(
  BookName := "Wedderga",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile := "doc/manual.pdf",
  SixFile := "doc/manual.six",
  LongTitle := "Wedderga with div-alg",
  Autoload := true
),

Dependencies := rec(
  GAP                    := ">=4.5",
  NeededOtherPackages    := [ ["GAPDoc", ">= 1.5.1"] ],
  SuggestedOtherPackages := [ ["laguna", "3.4"] ],
  ExternalConditions     := []
),

AvailabilityTest := ReturnTrue,
Autoload         := false,
TestFile        := "tst/testall.g",

Keywords := ["Wedderburn decomposition", "simple components", 
             "central idempotents", "group algebras"]

));
