#![doc="Defines the monoid algebraic structure.

A monoid is an algebraic structure with a single associative 
binary operation and an identity element. 
Monoids are studied in semigroup theory as they are 
semigroups with identity. 

A commutative monoid is a monoid whose binary operation is
commutative.

We define four kinds of monoids.

* Monoid with an addition operation with partial equivalence
* Monoid with an addition operation with full equivalence
* Monoid with a multiplication operation with partial equivalence
* Monoid with a multiplication operation with full equivalence


References:

* http://en.wikipedia.org/wiki/Algebraic_structure
* http://en.wikipedia.org/wiki/Monoid

"]

// std imports


// local imports
use num::traits::{Zero, One};
use algebra::structure::semigroup::{SemiGroupAddPartial, SemiGroupAdd, 
    SemiGroupMulPartial, SemiGroupMul};

///////////////////////////////////////////////////////////

/// Monoid with an addition operation with partial equivalence
pub trait MonoidAddPartial
    : SemiGroupAddPartial
    + Zero
{

}


impl<T> MonoidAddPartial for T where
    T: SemiGroupAddPartial + Zero 
{

}



///////////////////////////////////////////////////////////


/// Monoid with an addition operation with full equivalence
pub trait MonoidAdd 
    : SemiGroupAdd
    + MonoidAddPartial
{

}


impl<T> MonoidAdd for T where
    T: SemiGroupAdd + MonoidAddPartial
{}


///////////////////////////////////////////////////////////

/// Monoid with a multiplication operation with partial equivalence
pub trait MonoidMulPartial
    : SemiGroupMulPartial
    + One
{

}


impl<T> MonoidMulPartial for T where
    T: SemiGroupMulPartial + One 
{

}

///////////////////////////////////////////////////////////


/// Monoid with a multiplication operation with full equivalence
pub trait MonoidMul 
    : SemiGroupMul
    + MonoidMulPartial
{

}


impl<T> MonoidMul for T where
    T: SemiGroupMul + MonoidMulPartial
{}

///////////////////////////////////////////////////////////

/// Commutative monoid with an addition operation with partial equivalence
pub trait CommutativeMonoidAddPartial 
: MonoidAddPartial
{

    /// Returns `true` if the addition operator is approximately commutative for
    /// the given argument tuple.
    fn prop_is_commutative(a : Self, b : Self) -> bool {
        let ab = a.clone() + b.clone();
        let ba = b.clone() + a.clone();
        ab == ba
    }

}

impl CommutativeMonoidAddPartial for i8   {}
impl CommutativeMonoidAddPartial for i16  {}
impl CommutativeMonoidAddPartial for i32  {}
impl CommutativeMonoidAddPartial for i64  {}
impl CommutativeMonoidAddPartial for u8   {}
impl CommutativeMonoidAddPartial for u16  {}
impl CommutativeMonoidAddPartial for u32  {}
impl CommutativeMonoidAddPartial for u64  {}
impl CommutativeMonoidAddPartial for f32  {}
impl CommutativeMonoidAddPartial for f64  {}

///////////////////////////////////////////////////////////

/// Commutative monoid with an addition operation with full equivalence
pub trait CommutativeMonoidAdd
: CommutativeMonoidAddPartial + MonoidAdd
{

    /// Returns `true` if the addition operator is approximately commutative for
    /// the given argument tuple.
    fn prop_is_commutative(a : Self, b : Self) -> bool {
        let ab = a.clone() + b.clone();
        let ba = b.clone() + a.clone();
        ab == ba
    }

}

impl CommutativeMonoidAdd for i8   {}
impl CommutativeMonoidAdd for i16  {}
impl CommutativeMonoidAdd for i32  {}
impl CommutativeMonoidAdd for i64  {}
impl CommutativeMonoidAdd for u8   {}
impl CommutativeMonoidAdd for u16  {}
impl CommutativeMonoidAdd for u32  {}
impl CommutativeMonoidAdd for u64  {}

///////////////////////////////////////////////////////////


/// Commutative monoid with a multiplication operation with partial equivalence
pub trait CommutativeMonoidMulPartial 
: MonoidMulPartial
{

    /// Returns `true` if the multiplication operator is approximately commutative for
    /// the given argument tuple.
    fn prop_is_commutative(a : Self, b : Self) -> bool {
        let ab = a.clone() * b.clone();
        let ba = b.clone() * a.clone();
        ab == ba
    }

}

impl CommutativeMonoidMulPartial for i8   {}
impl CommutativeMonoidMulPartial for i16  {}
impl CommutativeMonoidMulPartial for i32  {}
impl CommutativeMonoidMulPartial for i64  {}
impl CommutativeMonoidMulPartial for u8   {}
impl CommutativeMonoidMulPartial for u16  {}
impl CommutativeMonoidMulPartial for u32  {}
impl CommutativeMonoidMulPartial for u64  {}
impl CommutativeMonoidMulPartial for f32  {}
impl CommutativeMonoidMulPartial for f64  {}

///////////////////////////////////////////////////////////


/// Commutative monoid with a multiplication operation with full equivalence
pub trait CommutativeMonoidMul
: CommutativeMonoidMulPartial + MonoidMul
{

    /// Returns `true` if the multiplication operator is approximately commutative for
    /// the given argument tuple.
    fn prop_is_commutative(a : Self, b : Self) -> bool {
        let ab = a.clone() * b.clone();
        let ba = b.clone() * a.clone();
        ab == ba
    }

}

impl CommutativeMonoidMul for i8   {}
impl CommutativeMonoidMul for i16  {}
impl CommutativeMonoidMul for i32  {}
impl CommutativeMonoidMul for i64  {}
impl CommutativeMonoidMul for u8   {}
impl CommutativeMonoidMul for u16  {}
impl CommutativeMonoidMul for u32  {}
impl CommutativeMonoidMul for u64  {}

///////////////////////////////////////////////////////////


#[cfg(test)]
mod tests {
    use super::*;
    use num::traits::{Zero, One};


    fn check_monoid_add_partial<T: MonoidAddPartial>(a: T, b: T)->T{
        let d = a + b;
        d + Zero::zero()
    }


    #[test]
    fn test_monoid_add_partial() {
        assert_eq!(check_monoid_add_partial(2i8, 3i8 ), 5);
        assert_eq!(check_monoid_add_partial(2i16, 3i16), 5);
        assert_eq!(check_monoid_add_partial(2i32, 3i32), 5);
        assert_eq!(check_monoid_add_partial(2i64, 3i64), 5);
        assert_eq!(check_monoid_add_partial(2u8, 3u8), 5);
        assert_eq!(check_monoid_add_partial(2u16, 3u16), 5);
        assert_eq!(check_monoid_add_partial(2u32, 3u32), 5);
        assert_eq!(check_monoid_add_partial(2u64, 3u64), 5);
        assert_eq!(check_monoid_add_partial(2f32, 3f32), 5f32);
        assert_eq!(check_monoid_add_partial(2f64, 3f64), 5f64);
        // The following will not compile
    }


    fn check_monoid_add<T: MonoidAdd>(a: T, b: T)->T{
        let d = a + b;
        d + Zero::zero()
    }

    #[test]
    fn test_monoid_add() {
        assert_eq!(check_monoid_add(2i8, 3i8), 5);
        assert_eq!(check_monoid_add(2i16, 3i16), 5);
        assert_eq!(check_monoid_add(2i32, 3i32), 5);
        assert_eq!(check_monoid_add(2i64, 3i64), 5);
        assert_eq!(check_monoid_add(2u8, 3u8), 5);
        assert_eq!(check_monoid_add(2u16, 3u16), 5);
        assert_eq!(check_monoid_add(2u32, 3u32), 5);
        assert_eq!(check_monoid_add(2u64, 3u64), 5);
        // The following will not compile
        // assert_eq!(check_monoid_add(2f32, 3f32), 5f32);
        // assert_eq!(check_monoid_add(2f64, 3f64, 1f64), -4f64);
    }


    fn check_monoid_mul_partial<T: MonoidMulPartial>(a: T, b: T)->T{
        let d = a * b;
        d * One::one()
    }


    #[test]
    fn test_monoid_mul_partial() {
        assert_eq!(check_monoid_mul_partial(2i8, 3i8), 6);
        assert_eq!(check_monoid_mul_partial(2i16, 3i16), 6);
        assert_eq!(check_monoid_mul_partial(2i32, 3i32), 6);
        assert_eq!(check_monoid_mul_partial(2i64, 3i64), 6);
        assert_eq!(check_monoid_mul_partial(2u8, 3u8), 6);
        assert_eq!(check_monoid_mul_partial(2u16, 3u16), 6);
        assert_eq!(check_monoid_mul_partial(2u32, 3u32), 6);
        assert_eq!(check_monoid_mul_partial(2u64, 3u64), 6);
        assert_eq!(check_monoid_mul_partial(2f32, 3f32), 6.0);
        assert_eq!(check_monoid_mul_partial(2f64, 3f64), 6.0);
        // The following will not compile
    }



}

