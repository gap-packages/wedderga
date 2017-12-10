gap> START_TEST( "ExtremeSSPs.tst");

# PrimitiveCentralIdempotentsByESSP(QG) and PrimitiveCentralIdempotentsByStrongSP(QG)
# All groups of order 2^5 being metabelian are normally monomial and hence two sets 
# must be equal
gap> for n in [1..NrSmallGroups(32)] do
>  G := SmallGroup(32,n);
>  QG:= GroupRing(Rationals,G);
>  if not IsEqualSet( PrimitiveCentralIdempotentsByESSP(QG),
>                     PrimitiveCentralIdempotentsByStrongSP(QG) ) then
>    Print("Error! Different PCIs for [32,", n, "]\n");
>  fi;
> od;

#PrimitiveCentralIdempotentsByESSP(QG) on non-normally monomial group
gap> G:=SmallGroup(486,38);     
<pc group of size 486 with 6 generators>
gap> QG:=GroupRing(Rationals,G);
<algebra-with-one over Rationals, with 6 generators>
gap> PrimitiveCentralIdempotentsByESSP(QG);;
Wedderga: Warning!!!
The output is a NON-COMPLETE list of prim. central idemp.s of the input! 

# All groups of order 2^5 being metabelian are normally monomial and hence sizes of two 
# sets  must be equal
gap> for n in [1..NrSmallGroups(32)] do
>  G := SmallGroup(32,n);
>  QG:= GroupRing(Rationals,G);
>  if not Size(ExtremelyStrongShodaPairs(G))= Size(StrongShodaPairs(G)) then
>    Print("Error! Different PCIs for [32,", n, "]\n");
>  fi;
> od;

# IsNormallyMonomial
gap> ForAll( [1..NrSmallGroups(32)], n -> IsNormallyMonomial(SmallGroup(32,n) ) );
true

# Among the groups of odd order up to 2000, the only groups which are not 
# normally monomial are below (https://doi.org/10.1016/j.jsc.2015.12.002)
gap> ForAll( [ [375,2], [1029,12], [1053,51], [1125,3], [1125,7], [1215,68],
> [1875,18], [1875,19] ], id-> not IsNormallyMonomial(SmallGroup(id) ) );
true

# non-metabelian but normally monomial group.
gap> ForAll([[1,1],[72,41],[192,1023],[192,1025]],id-> IsNormallyMonomial(SmallGroup(id)));
true

# nilpotent and hence strongly monomial but not normally monomial
gap> IsNormallyMonomial(SmallGroup([128,134]));
false
gap> IsStronglyMonomial(SmallGroup([128,134]));
true
gap> IsStronglyMonomial(SmallGroup([24,3]));
false

#strongly monomial but not normally monomial group.
gap> IsStronglyMonomial(SmallGroup([24,12]));
true
gap> IsNormallyMonomial(SmallGroup([24,12]));
false

#simple group, not normally monomial.
gap> IsNormallyMonomial(SmallGroup([60,5]));
false

#Normally monomial group, Uses is maximal function non-trivially
gap> IsNormallyMonomial(SmallGroup([486,36]));
false

# to check StrongShodaPairs(G);
gap> IdSample:=[[1,1],[24,12],[40,3],[60,5],[128,134],[256,52],[1000,86]];;
gap> for id in IdSample do
> G:=SmallGroup(id);
> S:=Size(StrongShodaPairs(G));
> Print(id," ",S,"\n");
> od;
[ 1, 1 ] 1
[ 24, 12 ] 5
[ 40, 3 ] 6
[ 60, 5 ] 1
[ 128, 134 ] 13
[ 256, 52 ] 21
[ 1000, 86 ] 7

# to check PrimitiveCentralIdempotentsByStrongSP(QG)
gap> for id in IdSample do
> G:=SmallGroup(id);
> QG:=GroupRing(Rationals,G);
> S:=Size(PrimitiveCentralIdempotentsByStrongSP(QG));
> Print(id," ",S,"\n");
> od;
[ 1, 1 ] 1
[ 24, 12 ] 5
[ 40, 3 ] 6
Wedderga: Warning!!!
The output is a NON-COMPLETE list of prim. central idemp.s of the input! 
[ 60, 5 ] 1
[ 128, 134 ] 13
[ 256, 52 ] 21
Wedderga: Warning!!!
The output is a NON-COMPLETE list of prim. central idemp.s of the input! 
[ 1000, 86 ] 7

# too long - leave for extended test:
# [128,2328], [256,56090], [640,21540], [640,21541], [768,7667],
# [896,19348],[896,19349],[972,875],[1000,86],[2000,399], [5^6,643]
gap> STOP_TEST( "ExtremeSSPs.tst", 1 );
