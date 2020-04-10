#############################################################################
##
#W  PackageInfo.g         The Wedderga package            Gurmeet Kaur Bakshi
#W                                                        Osnel Broche Cristo
#W                                                               Allen Herman
#W                                                        Alexander Konovalov
#W                                                            Aurora Olivieri
#W                                                           Gabriela Olteanu
#W                                                              Ángel del Río
#W                                                        Sugandha Maheshwary
#W                                                          Inneke Van Gelder
##
#############################################################################

SetPackageInfo( rec(

PackageName    := "Wedderga",
Subtitle       := Concatenation( [
                  "Wedderburn Decomposition of Group Algebras" ] ),
Version        := "4.9.5-divalg-dev",
Date           := "30/11/2018", # dd/mm/yyyy format
License        := "GPL-2.0-or-later",
##  <#GAPDoc Label="PKGVERSIONDATA">
##  <!ENTITY VERSION "4.9.5">
##  <!ENTITY RELEASEDATE "30 November 2018">
##  <!ENTITY RELEASEYEAR "2018">
##  <#/GAPDoc>

SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/gap-packages/", LowercaseString(~.PackageName) ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", LowercaseString(~.PackageName) ),
README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", LowercaseString(~.PackageName), "-", ~.Version ),

ArchiveFormats := ".tar.gz",

Persons :=
 [
      
  rec(
       LastName      := "Bakshi",
       FirstNames    := "Gurmeet Kaur",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "gkbakshi@pu.ac.in",
       PostalAddress := Concatenation( [
                        "Center for Advanced Study in Mathematics,\n",
                        "Panjab University,\n",
                        "Chandigarh,  \n",
                        "160014, India " ] ),
       Place         := "Chandigarh",
       Institution   := "Panjab University"
     ),
     rec(
       LastName      := "Broche Cristo",
       FirstNames    := "Osnel",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "osnel@ufla.br",
       PostalAddress := Concatenation( [
                        "Departamento de Ciências Exatas\n",
                        "Universidade Federal de Lavras - UFLA\n",
                        "Campus Universitário - Caixa Postal 3037\n",
                        "37200-000, Lavras - MG, Brazil" ] ),
       Place         := "Lavras - MG",
       Institution   := "Universidade Federal de Lavras - UFLA"
     ),
rec(
       LastName      := "Herman",
       FirstNames    := "Allen",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "aherman@math.uregina.ca",
       PostalAddress := Concatenation( [
                        "Department of Mathematics and Statistics\n",
                        "University of Regina\n",
                        "3737 Wascana Parkway\n",
                        "Regina, SK, S0G 0E0 Canada"] ),
       Place         := "Regina",
       Institution   := "University of Regina"
     ),
     rec(
       LastName      := "Konovalov",
       FirstNames    := "Alexander",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "alexander.konovalov@st-andrews.ac.uk",
       WWWHome       := "https://alexk.host.cs.st-andrews.ac.uk",
       PostalAddress := Concatenation( [
                        "School of Computer Science\n",
                        "University of St Andrews\n",
                        "Jack Cole Building, North Haugh,\n",
                        "St Andrews, Fife, KY16 9SX, Scotland" ] ),
       Place         := "St Andrews",
       Institution   := "University of St Andrews"
    ),
    
    rec(
       LastName      := "Maheshwary",
       FirstNames    := "Sugandha",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "sugandha@iisermohali.ac.in",
       PostalAddress := Concatenation( [
                        "Department of Mathematical Sciences\n",
                        "IISER Mohali\n",
                        "Knowledge city, \n",
                        "Sector 81, \n",
                        " SAS Nagar, Manauli, \n",
                        " PO 140306, India" ] ),
       Place         := "Mohali",
       Institution   := "Indian Institute of Science Education and Research Mohali"
     ),
    
    rec(
       LastName      := "Olteanu",
       FirstNames    := "Gabriela",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "gabriela.olteanu@econ.ubbcluj.ro",
       WWWHome       := "http://math.ubbcluj.ro/~olteanu",
       PostalAddress := Concatenation( [
                        "Department of Statistics-Forecasts-Mathematics\n",
                        "Faculty of Economics and Business Administration\n",
                        "Babes-Bolyai University\n",
                        "Str. T. Mihali 58-60, 400591 Cluj-Napoca, Romania" ] ),
       Place         := "Cluj-Napoca",
       Institution   := "Babes-Bolyai University"
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
     ),     
     rec(
       LastName      := "Van Gelder",
       FirstNames    := "Inneke",
       IsAuthor      := true,
       IsMaintainer  := true,
       Email         := "ivgelder@vub.ac.be",
       WWWHome       := "http://homepages.vub.ac.be/~ivgelder",
       PostalAddress := Concatenation( [
                        "Vrije Universiteit Brussel\n",
                        "Departement Wiskunde\n",
                        "Pleinlaan 2\n", 
                        "1050 Brussels , Belgium" ] ),
       Place         := "Brussels",
       Institution   := "Vrije Universiteit Brussel"
     )
],

Status := "accepted",
CommunicatedBy := "Gerhard Hiss (Aachen)",
AcceptDate := "01/2008",

AbstractHTML := "<span class=\"pkgname\">Wedderga</span> is the package to compute the simple components of the Wedderburn decomposition of semisimple group algebras of finite groups over finite fields and over subfields of finite cyclotomic extensions of the rationals. It also contains functions that produce the primitive central idempotents of semisimple group algebras and functions for computing Schur indices. Other functions of <span class=\"pkgname\">Wedderga</span> allow one to construct crossed products over a group with coefficients in an associative ring with identity and the multiplication determined by a given action and twisting.",
                  
PackageDoc := rec(
  BookName := "Wedderga",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile := "doc/manual.pdf",
  SixFile := "doc/manual.six",
  LongTitle := "Wedderga",
  Autoload := true
),

Dependencies := rec(
  GAP                    := ">=4.8",
  NeededOtherPackages    := [ ["GAPDoc", ">= 1.5.1"] ],
  SuggestedOtherPackages := [ ["laguna", ">= 3.4"], ["GUAVA", ">= 3.12"] ],
  ExternalConditions     := []
),

AvailabilityTest := ReturnTrue,
TestFile        := "tst/testall.g",

Keywords := ["Wedderburn decomposition", "simple components", 
             "central idempotents", "group algebras"]

));
