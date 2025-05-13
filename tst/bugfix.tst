# Andreas Bachle, Mauricio Caceido, and Inneke van Gelder reported a problem with 
# LocalIndexAtTwo for quaternion algebras over local fields in November 2014.  
# Allen Herman provided a bugfix November 21, 2014 that will be included in 
# Wedderga 4.7.2. 

gap> F:=Field([ER(-39)]);
NF(39,[ 1, 2, 4, 5, 8, 10, 11, 16, 20, 22, 25, 32 ])
gap> G:=SmallGroup(8,4);        
<pc group of size 8 with 3 generators>
gap> R:=GroupRing(F,G);        
<algebra-with-one over NF(39,[ 1, 2, 4, 5, 8, 10, 11, 16, 20, 22, 25, 32 
 ]), with 3 generators>
gap> W:=WedderburnDecompositionWithDivAlgParts(R);
[ [ 1, NF(39,[ 1, 2, 4, 5, 8, 10, 11, 16, 20, 22, 25, 32 ]) ], 
  [ 1, NF(39,[ 1, 2, 4, 5, 8, 10, 11, 16, 20, 22, 25, 32 ]) ], 
  [ 1, NF(39,[ 1, 2, 4, 5, 8, 10, 11, 16, 20, 22, 25, 32 ]) ], 
  [ 1, NF(39,[ 1, 2, 4, 5, 8, 10, 11, 16, 20, 22, 25, 32 ]) ], 
  [ 1, rec( Center := NF(39,[ 1, 2, 4, 5, 8, 10, 11, 16, 20, 22, 25, 32 ]), 
          DivAlg := true, LocalIndices := [ [ 2, 2 ] ], SchurIndex := 2 ) ] ]

# Andreas Bachle and Inneke van Gelder reported a problem when there is 
# a simple component that is a division algebra with index 2 whose local
# index at infinity is 1, caused by "DefiningCharacterOfCyclotomicAlgebra" 
# looking for a field of values that was larger than it should be for these 
# groups. This was fixed by Allen Herman in Wedderga 4.6.2 (September 2014).
gap> w := WedderburnDecompositionWithDivAlgParts(GroupRing(Rationals, SmallGroup(288,69)));;
gap> w1 := SortedList( Filtered( w, t -> not IsRecord(t[2]) ) );
[ [ 1, Rationals ], [ 1, Rationals ], [ 1, GaussianRationals ], 
  [ 2, Rationals ], [ 2, NF(8,[ 1, 3 ]) ], [ 2, NF(9,[ 1, 8 ]) ], 
  [ 3, Rationals ], [ 3, Rationals ], [ 3, GaussianRationals ], 
  [ 4, Rationals ], [ 4, NF(9,[ 1, 8 ]) ], [ 6, Rationals ] ]
gap> w2 := SortedList( Filtered( w, t -> IsRecord(t[2]) ) );
[ [ 1, 
      rec( Center := Rationals, DivAlg := true, 
          LocalIndices := [ [ 3, 2 ], [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 1, 
      rec( Center := NF(8,[ 1, 7 ]), DivAlg := true, 
          LocalIndices := [ [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 1, 
      rec( Center := NF(9,[ 1, 8 ]), DivAlg := true, 
          LocalIndices := [ [ 3, 2 ], [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 2, 
      rec( Center := Rationals, DivAlg := true, 
          LocalIndices := [ [ 3, 2 ], [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 2, 
      rec( Center := NF(9,[ 1, 8 ]), DivAlg := true, 
          LocalIndices := [ [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 3, 
      rec( Center := Rationals, DivAlg := true, 
          LocalIndices := [ [ 3, 2 ], [ infinity, 2 ] ], SchurIndex := 2 ) ] ]
gap> WedderburnDecompositionWithDivAlgParts(GroupRing(Rationals, SmallGroup(336,118)));
[ [ 1, Rationals ], [ 1, Rationals ], [ 2, Rationals ], [ 2, NF(7,[ 1, 6 ]) ],
  [ 2, NF(21,[ 1, 20 ]) ], [ 6, NF(7,[ 1, 6 ]) ], [ 3, Rationals ], 
  [ 3, Rationals ], 
  [ 1, 
      rec( Center := NF(8,[ 1, 7 ]), DivAlg := true, 
          LocalIndices := [ [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 2, 
      rec( Center := Rationals, DivAlg := true, 
          LocalIndices := [ [ 3, 2 ], [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 2, 
      rec( Center := NF(7,[ 1, 6 ]), DivAlg := true, 
          LocalIndices := [ [ infinity, 2 ], [ 2, 2 ] ], SchurIndex := 2 ) ], 
  [ 2, 
      rec( Center := NF(21,[ 1, 20 ]), DivAlg := true, 
          LocalIndices := [ [ infinity, 2 ] ], SchurIndex := 2 ) ] ]
gap> WedderburnDecompositionWithDivAlgParts(GroupRing(Rationals, SmallGroup(432,37))); 
[ [ 1, Rationals ], [ 1, Rationals ], [ 2, Rationals ], [ 2, NF(9,[ 1, 8 ]) ],
  [ 2, NF(27,[ 1, 26 ]) ], [ 6, Rationals ], [ 6, NF(9,[ 1, 8 ]) ], 
  [ 3, Rationals ], [ 3, Rationals ], 
  [ 1, 
      rec( Center := NF(8,[ 1, 7 ]), DivAlg := true, 
          LocalIndices := [ [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 2, 
      rec( Center := Rationals, DivAlg := true, 
          LocalIndices := [ [ 3, 2 ], [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 2, 
      rec( Center := NF(9,[ 1, 8 ]), DivAlg := true, 
          LocalIndices := [ [ infinity, 2 ] ], SchurIndex := 2 ) ], 
  [ 2, 
      rec( Center := NF(27,[ 1, 26 ]), DivAlg := true, 
          LocalIndices := [ [ infinity, 2 ] ], SchurIndex := 2 ) ] ]

# Fix for a bug reported by Ángel del Río on 14/11/2014
gap> G:=SmallGroup(16,7);;
gap> QG:=GroupRing(Rationals,G);;
gap> pci := PrimitiveCentralIdempotentsByCharacterTable(QG);;
gap> Length(pci);
6
gap> Length(pci)=Length(RationalClasses(G));
true

# Fix the issue with LocalIndexAtOddP and LocalIndexAtTwo (Feb 2020)
gap> G:=SmallGroup(81,3);
<pc group of size 81 with 4 generators>
gap> QG:=GroupRing(Rationals,G);
<algebra-with-one over Rationals, with 4 generators>
gap> WedderburnDecompositionWithDivAlgParts(QG);
[ [ 1, Rationals ], [ 1, CF(3) ], [ 1, CF(3) ], [ 1, CF(3) ], [ 1, CF(3) ],
  [ 1, CF(9) ], [ 1, CF(9) ], [ 1, CF(9) ], [ 3, CF(3) ], [ 3, CF(3) ],
  [ 3, CF(3) ] ]
gap> A:=[1,CF(7),4,[2,3,2]];
[ 1, CF(7), 4, [ 2, 3, 2 ] ]
gap> CyclotomicAlgebraWithDivAlgPart(A);
[ 1, rec( Center := CF(7), DivAlg := true, LocalIndices := [ [ 2, 2 ] ],
      SchurIndex := 2 ) ]

# Fix returning wrong result determining commutativity of a crossed product (issue #96)
gap> IsCommutative(WedderburnDecomposition(GroupRing(Rationals, QuaternionGroup(8)))[5]);
Error, no method found to check commutativity of a crossed product
