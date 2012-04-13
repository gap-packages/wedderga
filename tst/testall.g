LoadPackage("wedderga");

dir := DirectoriesPackageLibrary("wedderga","tst");

testsfiles := [ 
"wedderga01.tst",
"wedderga02.tst",
"wedderga03.tst",
"wedderga04.tst",
"wedderga05.tst",
"wedderga06.tst"
];

Print("=================================================================\n");
for ff in testsfiles do
  fn := Filename(dir, ff );
  Print("*** TESTING ", fn, "\n");
  ReadTest( fn );
  Print("=================================================================\n");
od;  
Print("*** FINISHED WEDDERGA PACKAGE TESTS\n");
Print("=================================================================\n");