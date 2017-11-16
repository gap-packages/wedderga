gap> START_TEST( "ExtremeSSPs.tst");
gap> for n in [1..NrSmallGroups(32)] do
>  G := SmallGroup(32,n);
>  QG:= GroupRing(Rationals,G);
>  if not IsEqualSet( PrimitiveCentralIdempotentsByStSP(QG),
>                     PrimitiveCentralIdempotentsByStrongSP(QG) ) or
>     not IsEqualSet( PrimitiveCentralIdempotentsByESSP(QG),
>                     PrimitiveCentralIdempotentsByStrongSP(QG) ) then
>    Print("Error! Different PCIs for [32,", n, "]\n");
>  fi;
> od;

# too long - leave for extended test:
# [128,2328], [256,56090], [640,21540], [640,21541], [768,7667],
# [896,19348],[896,19349],[972,875],[1000,86],[2000,399], [5^6,643]
gap> ids:=[[24,3],[24,12],[40,3],[54,8],[60,5],[72,41],[128,134],
> [192,1023],[192,1025],[256,52], [486,36],[486,38]];;
gap> for id in ids do
> G := SmallGroup(id);
> QG:= GroupRing(Rationals,G);
> if not IsEqualSet( PrimitiveCentralIdempotentsByStSP(QG),
>                    PrimitiveCentralIdempotentsByStrongSP(QG) ) then
>   Print("different PCIs for n=", n, "\n");
> fi;
> od;
Warning! The output is not complete list of pcis of the input! 
Wedderga: Warning!!!
The output is a NON-COMPLETE list of prim. central idemp.s of the input! 
Warning! The output is not complete list of pcis of the input! 
Wedderga: Warning!!!
The output is a NON-COMPLETE list of prim. central idemp.s of the input! 

# Checking IsNormallyMonomial
gap> ForAll( [1..NrSmallGroups(32)], n -> IsNormallyMonomial(SmallGroup(32,n) ) );
true

# Among the groups of odd order up to 2000, the only groups which are not 
# normally monomial are below (https://doi.org/10.1016/j.jsc.2015.12.002)
gap> ForAll( [ [375,2], [1029,12], [1053,51], [1125,3], [1125,7], [1215,68],
> [1875,18], [1875,19] ], id-> not IsNormallyMonomial(SmallGroup(id) ) );
true
gap> STOP_TEST( "ExtremeSSPs.tst");