#############################################################################
##
#R  IsCrossedProductObjDefaultRep( <obj> )
##
##  The default representation of an element object is a list of length 2,
##  at first position the zero coefficient, at second position a list with
##  the coefficients at the even positions, and the magma elements at the
##  odd positions, with the ordering as defined for the magma elements.
##
##  It is assumed that the arithmetic operations of $G$ produce only
##  normalized elements.
##
DeclareRepresentation( "IsCrossedProductObjDefaultRep", 
    IsPositionalObjectRep,
    [ 1, 2 ] );


#############################################################################
##
#M  ElementOfCrossedProduct( <Fam>, <zerocoeff>, <coeff>, <words> )
##
##  check whether <coeff> and <words> lie in the correct domains,
##  and remove zeroes.
##
InstallMethod( ElementOfCrossedProduct,
    "for family, ring element, and two homogeneous lists",
    [ IsFamily, IsRingElement, IsHomogeneousList, IsHomogeneousList ],
    function( Fam, zerocoeff, coeff, words )
    local rep, i, j;

    # Check that the data is admissible.
    if not IsBound( Fam!.defaultType ) then
      TryNextMethod();
    elif IsEmpty( coeff ) and IsEmpty( words ) then
      return Objectify( Fam!.defaultType, [ zerocoeff, [] ] );
    elif not IsIdenticalObj( FamilyObj( coeff ), Fam!.familyRing ) then
      Error( "<coeff> are not all in the correct domain" );
    elif not IsIdenticalObj( FamilyObj( words ), Fam!.familyGroup ) then
      Error( "<words> are not all in the correct domain" );
    elif Length( coeff ) <> Length( words ) then
      Error( "<coeff> and <words> must have same length" );
    fi;

    # Make sure that the list of words is strictly sorted.
    if not IsSSortedList( words ) then
      words:= ShallowCopy( words );
      coeff:= ShallowCopy( coeff );
      SortParallel( words, coeff );
      if not IsSSortedList( words ) then
        j:= 1;
        for i in [ 2 .. Length( coeff ) ] do
          if words[i] = words[j] then
            coeff[j]:= coeff[j] + coeff[i];
          else
            j:= j+1;
            words[j]:= words[i];
            coeff[j]:= coeff[i];
          fi;
        od;
        for i in [ j+1 .. Length( coeff ) ] do
          Unbind( words[i] );
          Unbind( coeff[i] );
        od;
      fi;
    fi;

    # Create the default representation, and remove zeros.
    rep:= [];
    j:= 1;
    for i in [ 1 .. Length( coeff ) ] do
      if coeff[i] <> zerocoeff then
        rep[  j  ]:= words[i];
        rep[ j+1 ]:= coeff[i];
        j:= j+2;
      fi;
    od;

    # Return the result
    return Objectify( Fam!.defaultType, [ zerocoeff, rep ] );
    end );


#############################################################################
##
#M  ZeroCoefficient( <elm> )
##
InstallMethod( ZeroCoefficient,
    "for crossed product element in default repr.",
    [ IsElementOfCrossedProduct and IsCrossedProductObjDefaultRep ],
    elm -> FamilyObj( elm )!.zeroRing );


#############################################################################
##
#M  CoefficientsAndMagmaElements( <elm> )
##
InstallMethod( CoefficientsAndMagmaElements,
    "for crossed product element in default repr.",
    [ IsElementOfCrossedProduct and IsCrossedProductObjDefaultRep ],
    elm -> elm![2] );


#############################################################################
##
#M  PrintObj( <elm> ) . . . . . for crossed product element in default repr.
##
InstallMethod( PrintObj,
    "for crossed product element",
    [ IsElementOfCrossedProduct ],
    function( elm )

    local coeffs_and_words,
          i;

    coeffs_and_words:= CoefficientsAndMagmaElements( elm );
    for i in [ 1, 3 .. Length( coeffs_and_words ) - 3 ] do
      Print( "(", coeffs_and_words[i], ")*", coeffs_and_words[i+1], "+" );
    od;
    i:= Length( coeffs_and_words );
    if i = 0 then
      Print( "<zero> of ..." );
    else
      Print( "(", coeffs_and_words[i-1], ")*", coeffs_and_words[i] );
    fi;
    end );


#############################################################################
##
#M  \=( <x>, <y> )  . . . . for two crossed product elements in default repr.
##
InstallMethod( \=,
    "for two crossed product elements",
    IsIdenticalObj,
    [ IsElementOfCrossedProduct,
      IsElementOfCrossedProduct ],
    function( x, y )
    return   CoefficientsAndMagmaElements( x )
           = CoefficientsAndMagmaElements( y );
    end );


#############################################################################
##
#M  \<( <x>, <y> )  . . . . for two crossed product elements in default repr.
##
InstallMethod( \<,
    "for two crossed product elements",
    IsIdenticalObj,
    [ IsElementOfCrossedProduct,
      IsElementOfCrossedProduct ],
    function( x, y )
    local i;
    x:= CoefficientsAndMagmaElements( x );
    y:= CoefficientsAndMagmaElements( y );
    for i in [ 1 .. Minimum( Length( x ), Length( y ) ) ] do
      if   x[i] < y[i] then
        return true;
      elif y[i] < x[i] then
        return false;
      fi;
    od;
    return Length( x ) < Length( y );
    end );


#############################################################################
##
#M  \+( <x>, <y> )  . . . . for two crossed product elements in default repr.
##
InstallMethod( \+,
    "for two crossed product elements",
    IsIdenticalObj,
    [ IsElementOfCrossedProduct,
      IsElementOfCrossedProduct ],
    function( x, y )
    local F, sum, z;
    F := FamilyObj( x );
    z := ZeroCoefficient( x );
    x := CoefficientsAndMagmaElements( x );
    y := CoefficientsAndMagmaElements( y );
    sum:= ZippedSum( x, y, z, [ \<, \+ ] );
    return Objectify( F!.defaultType, [ z, sum ] );
    end );


#############################################################################
##
#M  AdditiveInverseOp( <x> ) . . for crossed product element in default repr.
##
InstallMethod( AdditiveInverseOp,
    "for crossed product element",
    [ IsElementOfCrossedProduct ],
    function( x )
    local ext, i;
    ext:= ShallowCopy( CoefficientsAndMagmaElements( x ) );
    for i in [ 2, 4 .. Length( ext ) ] do
      ext[i]:= AdditiveInverse( ext[i] );
    od;
    return Objectify( FamilyObj( x )!.defaultType, [ ZeroCoefficient(x), ext] );
    end );


#############################################################################
##
#M  \*( <x>, <y> )  . . . . for two crossed product elements in default repr.
##
InstallMethod( \*,
    "for two crossed product elements",
    IsIdenticalObj,
    [ IsElementOfCrossedProduct,
      IsElementOfCrossedProduct ],
    function( x, y )
    local F, prod, z,  mons,  cofs,  i,  j,  c, Twisting, Action;
    F := FamilyObj( x );
    Twisting := F!.twisting;
    Action := F!.action;
    z := ZeroCoefficient( x );
    x := CoefficientsAndMagmaElements( x );
    y := CoefficientsAndMagmaElements( y );

    # fold the product
    mons := [];
    cofs := [];
    for i  in [ 1, 3 .. Length(x)-1 ]  do
      for j  in [ 1, 3 .. Length(y)-1 ]  do
	      # we compute product of the coefficients as follows 
	      # (x * a1) * (y * a2) = (x) * (y) * a1^action(y) * a2 = 
        # (x * y) * twisting(x,y) *a1^action(y) * a2 
        c := Twisting(x[i],y[j]) * ( x[i+1]^Action(y[i]) * y[j+1] );
	      if c <> z  then
	        ##  add the product of the monomials
	        Add( mons, x[i] * y[j] );
	        ##  and the coefficient
	        Add( cofs, c );
	      fi;
      od;
    od;

    # sort monomials
    SortParallel( mons, cofs );

    # sum coeffs
    prod := [];
    i := 1;
    while i <= Length(mons)  do
      c := cofs[i];
      while i < Length(mons) and mons[i] = mons[i+1]  do
	      i := i+1;
	      c := c + cofs[i];    ##  add coefficients
      od;
      if c <> z  then
	      ## add the term to the product
	      Add( prod, mons[i] );
	      Add( prod, c );
      fi;
      i := i+1;
    od;

    return Objectify( F!.defaultType, [ z, prod ] );
    end );


#############################################################################
##
#M  OneOp( <elm> )
##
InstallMethod( OneOp,
    "for crossed product element",
    [ IsElementOfCrossedProduct ],
    function( elm )
    local F, z;
    F:= FamilyObj( elm );
    if not IsBound( F!.oneGroup ) then
      return fail;
    fi;
    z:= ZeroCoefficient( elm );
    return Objectify( F!.defaultType, [ z, [ F!.oneGroup, One( z ) ] ] );
    end );


#############################################################################
##
#M  ZeroOp( <elm> )
##
InstallMethod( ZeroOp,
    "for crossed product element",
    [ IsElementOfCrossedProduct ],
    x -> Objectify( FamilyObj(x)!.defaultType,
             [ ZeroCoefficient( x ), [] ] ) );


#############################################################################
##
#F  CrossedProduct( <R>, <G>, act, twist )
##
## An example of trivial action and twisting:
## action should return a mapping F-->F that can be applied via "^" operation
##
##   function(a)
##     return IdentityMapping(F);
##   end,
##
## twisting should return an (invertible) element of F
## 
##   function( g, h)
##     return One(F);
##   end );
##
## to be used in the following way:
##
##   g * h = g * h * twisting(g,h)    for g,h in G
##   a * g = g * a^action(g)          for a in R and g in G
##
InstallGlobalFunction( CrossedProduct, 
function( R, G, act, twist )
    local filter,  # implied filter of all elements in the new domain
          F,       # family of crossed product elements
          one,     # identity of `R'
          zero,    # zero of `R'
          m,       # one element of `G'
          RG,      # free magma ring, result
          gens;    # generators of the magma ring

    # Check the arguments.
    if not IsRing( R ) or One( R ) = fail then
      Error( "<R> must be a ring with identity" );
    fi;
    
    if not IsGroup( G ) then
      Error( "<G> must be a group" );
    fi;

    F:= NewFamily( "CrossedProductObjFamily",
                   IsElementOfCrossedProduct,
                   IsMultiplicativeElementWithInverse and
                   IsAssociativeElement );

    one:= One( R );
    zero:= Zero( R );

    F!.defaultType := NewType( F, IsCrossedProductObjDefaultRep );
    F!.familyRing  := FamilyObj( R );
    F!.familyGroup := FamilyObj( G );
    F!.zeroRing    := zero;
    F!.oneGroup    := One( G );
    F!.action      := act;
    F!.twisting    := twist;

    # Set the characteristic.
    if HasCharacteristic( R ) or HasCharacteristic( FamilyObj( R ) ) then
      SetCharacteristic( F, Characteristic( R ) );
    fi;

    RG:= Objectify( NewType( CollectionsFamily( F ),
                            IsCrossedProduct and IsAttributeStoringRep ),
                      rec() );

    # Set the necessary attributes.
    SetLeftActingDomain( RG, R );
    SetUnderlyingGroup(  RG, G );
    SetIsAssociative( RG, true );
       
    # Deduce other useful information.
    if HasIsFinite( G ) then
      SetIsFiniteDimensional( RG, IsFinite( G ) );
    fi;
    
    # What about IsCommutative ? In MagmaRings it is as below:   
    # if HasIsCommutative( R ) and HasIsCommutative( G ) then
    #   SetIsCommutative( RG, IsCommutative( R ) and IsCommutative( G ) );
    # fi;
    
    if HasIsWholeFamily( R ) and HasIsWholeFamily( G ) then
      SetIsWholeFamily( RG, IsWholeFamily( R ) and IsWholeFamily( G ) );
    fi;

    # Construct the generators. To get meaningful generators, 
    # we have to handle the case that the groups is trivial.
    
    gens:= GeneratorsOfGroup( G );
    if IsEmpty( gens ) then
      SetGeneratorsOfLeftOperatorRingWithOne( RG,
              [ ElementOfCrossedProduct( F, zero, [ one ], [ One( G ) ] ) ] );
    else
      SetGeneratorsOfLeftOperatorRingWithOne( RG,
        List( gens,
                x -> ElementOfCrossedProduct( F, zero, [ one ], [ x ] ) ) );
    fi;

    # Return the crossed product
    return RG;
end );


#############################################################################
##
#M  ViewObj( <RG> ) . . . . . . . . . . . . . . . . . . for a crossed product
##
InstallMethod( ViewObj,
    "for a crossed product",
    [ IsCrossedProduct ],
    10,
    function( RG )
    Print( "<crossed product over ", LeftActingDomain( RG ), ", with ", 
           Length(GeneratorsOfGroup(UnderlyingGroup(RG))), " generators>" );
    end );


#############################################################################
##
#M  PrintObj( <RG> )  . . . . . . . . . . . . . . . . . for a crossed product
##
InstallMethod( PrintObj,
    "for a crossed product",
    [ IsCrossedProduct ],
    10,
    function( RG )
    Print( "CrossedProduct( ", LeftActingDomain( RG ), ", ",
                               UnderlyingGroup(  RG ), " )" );
    end );


#############################################################################
##
#R  IsCanonicalBasisCrossedProductRep( <B> )
##
DeclareRepresentation( "IsCanonicalBasisCrossedProductRep",
    IsCanonicalBasis and IsAttributeStoringRep,
    [ "zerovector" ] );

        
#############################################################################
##
#M  Coefficients( <B>, <v> )  . . . . . . for canon. basis of crossed product
##
InstallMethod( Coefficients,
    "for canon. basis of a crossed product, and an element of a crossed product",
    IsCollsElms,
    [ IsCanonicalBasisCrossedProductRep, IsElementOfCrossedProduct ],
    function( B, v )

    local coeffs,
          data,
          elms,
          i;

    data:= CoefficientsAndMagmaElements( v );
    coeffs:= ShallowCopy( B!.zerovector );
    elms:= EnumeratorSorted( UnderlyingGroup( UnderlyingLeftModule( B ) ) );
    for i in [ 1, 3 .. Length( data )-1 ] do
      coeffs[ Position( elms, data[i] ) ]:= data[i+1];
    od;
    return coeffs;
    end );
    
    
#############################################################################
##
#M  Basis( <RG> ) . . . . . . . . . . . . . . . . . . . for a crossed product
##
InstallMethod( Basis,
    "for a crossed product (delegate to `CanonicalBasis')",
    [ IsCrossedProduct ], CANONICAL_BASIS_FLAGS,
    CanonicalBasis );


#############################################################################
##
#M  CanonicalBasis( <RG> )  . . . . . . . . . . . . . . for a crossed product
##
InstallMethod( CanonicalBasis,
    "for a crossed product",
    [ IsCrossedProduct ],
    function( RG )

    local B, one, zero, F;

    F:= ElementsFamily( FamilyObj( RG ) );
    if not IsBound( F!.defaultType ) then
      TryNextMethod();
    fi;

    one  := One(  LeftActingDomain( RG ) );
    zero := Zero( LeftActingDomain( RG ) );

    B:= Objectify( NewType( FamilyObj( RG ),
                                IsFiniteBasisDefault
                            and IsCanonicalBasisCrossedProductRep ),
                   rec() );

    SetUnderlyingLeftModule( B, RG );
    if IsFiniteDimensional( RG ) then
      SetBasisVectors( B,
          List( EnumeratorSorted( UnderlyingGroup( RG ) ),
                x -> ElementOfCrossedProduct( F, zero, [ one ], [ x ] ) ) );
      B!.zerovector:= List( BasisVectors( B ), x -> zero );
    fi;

    return B;
    end );


#############################################################################
##
#M  IsFinite( <RG> )  . . . . . . . . . . . . . . . . . for a crossed product
##
InstallMethod( IsFinite,
    "for a crossed product",
    [ IsCrossedProduct ],
    RG -> IsFinite( LeftActingDomain( RG ) ) and 
          IsFinite( UnderlyingGroup( RG ) ) );


#############################################################################
##
#M  Representative( <RG> )  . . . . . . . . . . . . . . for a crossed product
##
##  this is a quick-hack solution, should be replaced
##  
InstallMethod( Representative,
    "for a crossed product",
    [ IsCrossedProduct ],
    RG -> GeneratorsOfLeftOperatorRingWithOne(RG)[1] );
   

#############################################################################
##
#M  IsFiniteDimensional( <RG> ) . . . . . . . . . . . . for a crossed product
##
# InstallMethod( IsFiniteDimensional,
#     "for a crossed product",
#     [ IsCrossedProduct ],
#     RG -> IsFinite( UnderlyingGroup( RG ) ) );
