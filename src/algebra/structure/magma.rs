#![doc="Defines the magma algebraic structure.




A magma is a basic kind of algebraic structure.
A magma is a set with a single binary operation.
The binary operation must be closed. No other
properties are imposed.



We define four kinds of magmas.

* Magma with an addition operation with partial equivalence
* Magma with an addition operation with full equivalence
* Magma with a multiplication operation with partial equivalence
* Magma with a multiplication operation with full equivalence

Types like float don't support full equivalence. That's why
separate magma types are needed for them.

We define separate magma types for addition and multiplication
operations. This gives us freedom in combining them to form
algebraic structures with two operations. 

In a sense, the traits defined in this module are marker traits.
It is not really possible to verify the closure property for a 
type at the compile time. It is assumed that the implementation
will ensure that (to a reasonable degree).

We also define a base trait ``MagmaBase``. This trait imposes
some basic requirements for all types taking advantage
of the algebraic traits defined in ``SciRust``. They all
are required to support ``Debug``, ``Clone`` and ``Sized`` traits.


References:

* http://en.wikipedia.org/wiki/Algebraic_structure

* http://en.wikipedia.org/wiki/Magma_(algebra)


"]



// std imports
use std::fmt::Debug;
// use std::fmt::Display;
// use std::marker::MarkerTrait;
use std::ops::{Add, Mul};


/// Defines basic requirements for all types implementing
/// the algebraic traits defined  in SciRust
pub trait MagmaBase : Debug + Copy + Clone + Sized + PartialEq{

}

impl<T> MagmaBase for T where
    T : Debug + Clone + Sized + Copy + PartialEq
{

}

/// Magma with an addition operation with partial equivalence
pub trait MagmaAddPartial 
    : MagmaBase 
    + Add<Output=Self> {

}

impl<T> MagmaAddPartial for T where
    T: MagmaBase + Add<Output=T>,
{}

///////////////////////////////////////////////////////////


/// Magma with an addition operation with full equivalence
pub trait MagmaAdd 
    : MagmaAddPartial
    + Eq {

}

impl<T> MagmaAdd for T where
    T: MagmaAddPartial + Eq,
{}

///////////////////////////////////////////////////////////

/// Magma with a multiplication operation with partial equivalence
pub trait MagmaMulPartial
    : MagmaBase 
    + Mul<Output=Self>{

}

impl<T> MagmaMulPartial for T where
    T: MagmaBase + Mul<Output=T>,
{}


///////////////////////////////////////////////////////////

/// Magma with a multiplication operation with full equivalence
pub trait MagmaMul 
    : MagmaMulPartial
    + Eq {

}

impl<T> MagmaMul for T where
    T: MagmaMulPartial + Eq,
{}


///////////////////////////////////////////////////////////


#[cfg(test)]
mod tests {
    use super::*;

    fn check_magma_add_partial<T: MagmaAddPartial>(a: T, b: T)->T{
        a + b
    }

    #[test]
    fn test_magma_add_partial() {
        assert_eq!(check_magma_add_partial(2i8, 3i8), 5);
        assert_eq!(check_magma_add_partial(2u8, 3u8), 5);
        assert_eq!(check_magma_add_partial(2i16, 3i16), 5);
        assert_eq!(check_magma_add_partial(2u16, 3u16), 5);
        assert_eq!(check_magma_add_partial(2i32, 3i32), 5);
        assert_eq!(check_magma_add_partial(2u32, 3u32), 5);
        assert_eq!(check_magma_add_partial(2f32, 3f32), 5f32);
        assert_eq!(check_magma_add_partial(2f64, 3f64), 5f64);
    }

    fn check_magma_add<T: MagmaAdd>(a: T, b: T)->T{
        a + b
    }

    #[test]
    fn test_magma_add() {
        assert_eq!(check_magma_add(2i8, 3i8), 5);
        assert_eq!(check_magma_add(2u8, 3u8), 5);
        assert_eq!(check_magma_add(2i16, 3i16), 5);
        assert_eq!(check_magma_add(2u16, 3u16), 5);
        assert_eq!(check_magma_add(2i32, 3i32), 5);
        assert_eq!(check_magma_add(2u32, 3u32), 5);
        // The following lines won't compile
        //assert_eq!(check_magma_add(2f32, 3f32), 5f32);
        //assert_eq!(check_magma_add(2f64, 3f64), 5f64);
    }

    fn check_magma_mul_partial<T: MagmaMulPartial>(a: T, b: T)->T{
        a * b
    }

    #[test]
    fn test_magma_mul_partial() {
        assert_eq!(check_magma_mul_partial(2i8, 3i8), 6);
        assert_eq!(check_magma_mul_partial(2u8, 3u8), 6);
        assert_eq!(check_magma_mul_partial(2i16, 3i16), 6);
        assert_eq!(check_magma_mul_partial(2u16, 3u16), 6);
        assert_eq!(check_magma_mul_partial(2i32, 3i32), 6);
        assert_eq!(check_magma_mul_partial(2u32, 3u32), 6);
        assert_eq!(check_magma_mul_partial(2f32, 3f32), 6f32);
        assert_eq!(check_magma_mul_partial(2f64, 3f64), 6f64);
    }

    fn check_magma_mul<T: MagmaMul>(a: T, b: T)->T{
        a * b
    }


    #[test]
    fn test_magma_mul() {
        assert_eq!(check_magma_mul(2i8, 3i8), 6);
        assert_eq!(check_magma_mul(2u8, 3u8), 6);
        assert_eq!(check_magma_mul(2i16, 3i16), 6);
        assert_eq!(check_magma_mul(2u16, 3u16), 6);
        assert_eq!(check_magma_mul(2i32, 3i32), 6);
        assert_eq!(check_magma_mul(2u32, 3u32), 6);
        // The following lines won't compile
        //assert_eq!(check_magma_mul(2f32, 3f32), 6f32);
        //assert_eq!(check_magma_mul(2f64, 3f64), 6f64);
    }

}
