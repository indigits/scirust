#![doc="Defines the quasi-group algebraic structure.

 A quasigroup is an algebraic structure resembling a group 
 in the sense that division is always possible. 
 Quasigroups differ from groups mainly in that they need not 
 be associative.


We define four kinds of quasigroups.

* Quasigroup with an addition operation with partial equivalence
* Quasigroup with an addition operation with full equivalence
* Quasigroup with a multiplication operation with partial equivalence
* Quasigroup with a multiplication operation with full equivalence



References:

* http://en.wikipedia.org/wiki/Algebraic_structure
* http://en.wikipedia.org/wiki/Quasigroup


"]


// std imports
use std::ops::{Sub};
use ops::Division;


// local imports
use magma::{MagmaAddPartial, MagmaAdd,
    MagmaMulPartial, MagmaMul};

///////////////////////////////////////////////////////////

/// Quasigroup with an addition operation with partial equivalence
pub trait QuasiGroupAddPartial
    : MagmaAddPartial
    + Sub<Output=Self>
{

}


impl QuasiGroupAddPartial for u8   {}
impl QuasiGroupAddPartial for u16  {}
impl QuasiGroupAddPartial for u32  {}
impl QuasiGroupAddPartial for u64  {}
impl QuasiGroupAddPartial for i8   {}
impl QuasiGroupAddPartial for i16  {}
impl QuasiGroupAddPartial for i32  {}
impl QuasiGroupAddPartial for i64  {}
impl QuasiGroupAddPartial for f32  {}
impl QuasiGroupAddPartial for f64  {}


///////////////////////////////////////////////////////////


/// Quasigroup with an addition operation with full equivalence
pub trait QuasiGroupAdd 
    : MagmaAdd
    + QuasiGroupAddPartial
{

}


impl QuasiGroupAdd for u8   {}
impl QuasiGroupAdd for u16  {}
impl QuasiGroupAdd for u32  {}
impl QuasiGroupAdd for u64  {}
impl QuasiGroupAdd for i8   {}
impl QuasiGroupAdd for i16  {}
impl QuasiGroupAdd for i32  {}
impl QuasiGroupAdd for i64  {}

///////////////////////////////////////////////////////////

/// Quasigroup with a multiplication operation with partial equivalence
pub trait QuasiGroupMulPartial
    : MagmaMulPartial
    + Division<Output=Self>
{

}




///////////////////////////////////////////////////////////

/// Quasigroup with a multiplication operation with full equivalence
pub trait QuasiGroupMul 
    : MagmaMul
    + QuasiGroupMulPartial
{

}


///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;

    fn check_quasigroup_add_partial<T: QuasiGroupAddPartial>(a: T, b: T, c : T)->T{
        let d = a + b;
        let e = d - c;
        e
    }


    #[test]
    fn test_quasigroup_add_partial() {
        assert_eq!(check_quasigroup_add_partial(2i8, 3i8, 1i8), 4);
        assert_eq!(check_quasigroup_add_partial(2i16, 3i16, 1i16), 4);
        assert_eq!(check_quasigroup_add_partial(2i32, 3i32, 1i32), 4);
        assert_eq!(check_quasigroup_add_partial(2i64, 3i64, 1i64), 4);
        assert_eq!(check_quasigroup_add_partial(2f32, 3f32, 1f32), 4f32);
        assert_eq!(check_quasigroup_add_partial(2f64, 3f64, 1f64), 4f64);
        // The following will not compile
        // assert_eq!(check_quasigroup_add_partial(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_quasigroup_add_partial(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_quasigroup_add_partial(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_quasigroup_add_partial(2u64, 3u64, 1u64), -4);
    }


    fn check_quasigroup_add<T: QuasiGroupAdd>(a: T, b: T, c : T)->T{
        let d = a + b;
        let e = d - c;
        e
    }


    #[test]
    fn test_quasigroup_add() {
        assert_eq!(check_quasigroup_add(2i8, 3i8, 1i8), 4);
        assert_eq!(check_quasigroup_add(2i16, 3i16, 1i16), 4);
        assert_eq!(check_quasigroup_add(2i32, 3i32, 1i32), 4);
        assert_eq!(check_quasigroup_add(2i64, 3i64, 1i64), 4);
        // The following will not compile
        // assert_eq!(check_quasigroup_add(2u8, 3u8, 1u8), 4);
        // assert_eq!(check_quasigroup_add(2u16, 3u16, 1u16), 4);
        // assert_eq!(check_quasigroup_add(2u32, 3u32, 1u32), 4);
        // assert_eq!(check_quasigroup_add(2u64, 3u64, 1u64), 4);
        //assert_eq!(check_quasigroup_add(2f32, 3f32, 1f32), 4f32);
        //assert_eq!(check_quasigroup_add(2f64, 3f64, 1f64), 4f64);
    }


    #[allow(dead_code)]
    fn check_quasigroup_mul_partial<T: QuasiGroupMulPartial>(a: T, b: T, c : T)->T{
        let d = a * b;
        let e = d / c;
        e
    }


    #[test]
    fn test_quasigroup_mul_partial() {
        //assert_eq!(check_quasigroup_mul_partial(2f32, 3f32, 3f32), 2_f32);
        //assert_eq!(check_quasigroup_mul_partial(2f64, 3f64, 3f64), 2_f64);
        // The following will not compile
        // assert_eq!(check_quasigroup_mul_partial(2i8, 3i8, 1i8), -4);
        // assert_eq!(check_quasigroup_mul_partial(2i16, 3i16, 1i16), -4);
        // assert_eq!(check_quasigroup_mul_partial(2i32, 3i32, 1i32), -4);
        // assert_eq!(check_quasigroup_mul_partial(2i64, 3i64, 1i64), -4);
        // assert_eq!(check_quasigroup_mul_partial(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_quasigroup_mul_partial(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_quasigroup_mul_partial(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_quasigroup_mul_partial(2u64, 3u64, 1u64), -4);
    }



}
