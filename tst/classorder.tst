# A character table may list the conjugacy classes in another order than
# ConjugacyClasses of its group, and Irr need not start with the trivial
# character (GAP reference manual, "The Interface between Character Tables
# and Groups").
gap> START_TEST( "classorder.tst" );

# IdLast( G ) is a copy of G whose class list ends with the identity class.
# Its character table starts with the identity class.
gap> IdLast := function( G )
>   local reps, H;
>   reps := List( ConjugacyClasses( G ), Representative );
>   reps := Concatenation( Filtered( reps, x -> not IsOne( x ) ), [ One( G ) ] );
>   H := GroupWithGenerators( GeneratorsOfGroup( G ) );
>   SetConjugacyClasses( H, List( reps, x -> ConjugacyClass( H, x ) ) );
>   return H;
> end;;

# LinCharByKernel
gap> K := IdLast( SymmetricGroup( 3 ) );;
gap> chi := LinCharByKernel( K, DerivedSubgroup( K ) );;
gap> ForAll( K, g -> g ^ chi = SignPerm( g ) );
true
gap> K := IdLast( AlternatingGroup( 4 ) );;
gap> H := DerivedSubgroup( K );;
gap> chi := LinCharByKernel( K, H );;
gap> IsCharacter( chi ) and chi[1] = 1 and KernelOfCharacter( chi ) = H;
true
gap> Set( LinCharByKernel( K, K ) );
[ 1 ]

# GaloisRepsOfCharacters, for a copy of A4 with an irrational character
# first in Irr
gap> G := AlternatingGroup( 4 );;
gap> H := GroupWithGenerators( GeneratorsOfGroup( G ) );;
gap> tbl := CharacterTable( H );;
gap> irr := List( Irr( G ), chi -> Character( tbl,
>      List( ConjugacyClasses( tbl ), c -> Representative( c ) ^ chi ) ) );;
gap> SetIrr( tbl, Concatenation(
>      Filtered( irr, chi -> not ForAll( chi, IsRat ) ),
>      Filtered( irr, chi -> ForAll( chi, IsRat ) ) ) );
gap> ForAll( Irr( H )[1], IsRat );
false
gap> Length( GaloisRepsOfCharacters( Rationals, H ) );
3
gap> Length( GaloisRepsOfCharacters( CF(3), H ) );
4
gap> wd := WedderburnDecompositionByCharacterDescent( Rationals, H );;
gap> Sum( wd, x -> x[1]^2 * Dimension( x[2] ) );
12

# DefectGroupsOfPBlock: orders of the defect groups of the blocks of S5
gap> G := IdLast( SymmetricGroup( 5 ) );;
gap> List( [ 2, 3, 5 ], p -> SortedList( List( Irr( G ), chi ->
>      Size( Representative( DefectGroupsOfPBlock( G, chi, p ) ) ) ) ) );
[ [ 2, 2, 8, 8, 8, 8, 8 ], [ 1, 3, 3, 3, 3, 3, 3 ], [ 1, 1, 5, 5, 5, 5, 5 ] ]
gap> SortedList( List( Irr( G ), chi -> DefectOfCharacterAtP( G, chi, 3 ) ) );
[ 0, 1, 1, 1, 1, 1, 1 ]

#
gap> STOP_TEST( "classorder.tst", 1 );
