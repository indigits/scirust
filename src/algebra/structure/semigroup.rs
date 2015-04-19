#![doc="Defines the semigroup algebraic structure.

A semigroup  is an algebraic structure consisting of a set 
together with an associative binary operation. 

Semigroup builds on top of magma and provides associativity.


We define four kinds of semigroups.

* Semigroup with an addition operation with partial equivalence
* Semigroup with an addition operation with full equivalence
* Semigroup with a multiplication operation with partial equivalence
* Semigroup with a multiplication operation with full equivalence

It is not possible to check the associativity of the
group operation at the compile time. We do provide 
a function (with each semigroup trait) to validate
the associativity of the operation.



References:

* http://en.wikipedia.org/wiki/Algebraic_structure
* http://en.wikipedia.org/wiki/Semigroup

"]


// local imports
use algebra::structure::magma::{MagmaAddPartial, MagmaAdd,
    MagmaMulPartial, MagmaMul};


///////////////////////////////////////////////////////////

/// Semigroup with an addition operation with partial equivalence
pub trait SemiGroupAddPartial
    : MagmaAddPartial
{

    fn prop_is_associative(a : Self, b : Self, c : Self) -> bool {
        let ab = a.clone() + b.clone();
        let bc = b.clone() + c.clone();
        (ab) + c == a + (bc)
    }


}

impl SemiGroupAddPartial for u8 {}
impl SemiGroupAddPartial for u16  {}
impl SemiGroupAddPartial for u32  {}
impl SemiGroupAddPartial for u64  {}
impl SemiGroupAddPartial for i8   {}
impl SemiGroupAddPartial for i16  {}
impl SemiGroupAddPartial for i32  {}
impl SemiGroupAddPartial for i64  {}
impl SemiGroupAddPartial for f32  {}
impl SemiGroupAddPartial for f64  {}

///////////////////////////////////////////////////////////


/// Semigroup with an addition operation with full equivalence
pub trait SemiGroupAdd 
    : MagmaAdd
    + SemiGroupAddPartial
{
    fn prop_is_associative(a : Self, b : Self, c : Self) -> bool {
        let ab = a.clone() + b.clone();
        let bc = b.clone() + c.clone();
        (ab) + c == a + (bc)
    }
}


impl SemiGroupAdd for u8 {}
impl SemiGroupAdd for u16  {}
impl SemiGroupAdd for u32  {}
impl SemiGroupAdd for u64  {}
impl SemiGroupAdd for i8   {}
impl SemiGroupAdd for i16  {}
impl SemiGroupAdd for i32  {}
impl SemiGroupAdd for i64  {}
// The following won't compile.
// impl SemiGroupAdd for f32  {}
// impl SemiGroupAdd for f64  {}

///////////////////////////////////////////////////////////

/// Semigroup with a multiplication operation with partial equivalence
pub trait SemiGroupMulPartial
    : MagmaMulPartial
{

    fn prop_is_associative(a : Self, b : Self, c : Self) -> bool {
        let ab = a.clone() * b.clone();
        let bc = b.clone() * c.clone();
        (ab) * c == a * (bc)
    }


}

impl SemiGroupMulPartial for u8 {}
impl SemiGroupMulPartial for u16  {}
impl SemiGroupMulPartial for u32  {}
impl SemiGroupMulPartial for u64  {}
impl SemiGroupMulPartial for i8   {}
impl SemiGroupMulPartial for i16  {}
impl SemiGroupMulPartial for i32  {}
impl SemiGroupMulPartial for i64  {}
impl SemiGroupMulPartial for f32  {}
impl SemiGroupMulPartial for f64  {}

///////////////////////////////////////////////////////////

/// Semigroup with a multiplication operation with full equivalence
pub trait SemiGroupMul 
    : SemiGroupMulPartial + MagmaMul
{

}

impl<T> SemiGroupMul for T where
    T: MagmaMul + SemiGroupMulPartial
{}

