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
use std::ops::{Sub, Div};


// local imports
use algebra::magma::{MagmaAddPartial, MagmaAdd,
    MagmaMulPartial, MagmaMul};

///////////////////////////////////////////////////////////

/// Quasigroup with an addition operation with partial equivalence
pub trait QuasiGroupAddPartial
    : MagmaAddPartial
    + Sub<Output=Self>
{

}


impl<T> QuasiGroupAddPartial for T where
    T: MagmaAddPartial
      + Sub<Output=T>,
{}


///////////////////////////////////////////////////////////


/// Quasigroup with an addition operation with full equivalence
pub trait QuasiGroupAdd 
    : MagmaAdd
    + QuasiGroupAddPartial
{

}


impl<T> QuasiGroupAdd for T where
    T: MagmaAdd + QuasiGroupAddPartial,
{}


///////////////////////////////////////////////////////////

/// Quasigroup with a multiplication operation with partial equivalence
pub trait QuasiGroupMulPartial
    : MagmaMulPartial
    + Div<Output=Self>
{

}


impl<T> QuasiGroupMulPartial for T where
    T: MagmaMulPartial 
      + Div<Output=T>,
{}

///////////////////////////////////////////////////////////

/// Quasigroup with a multiplication operation with full equivalence
pub trait QuasiGroupMul 
    : MagmaMul
    + QuasiGroupMulPartial
{

}

impl<T> QuasiGroupMul for T where
    T: MagmaMul + QuasiGroupMulPartial
{}



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



    fn check_quasigroup_mul_partial<T: QuasiGroupMulPartial>(a: T, b: T, c : T)->T{
        let d = a * b;
        let e = d / c;
        e
    }


    #[test]
    fn test_quasigroup_mul_partial() {
        assert_eq!(check_quasigroup_mul_partial(2f32, 3f32, 3f32), 2_f32);
        assert_eq!(check_quasigroup_mul_partial(2f64, 3f64, 3f64), 2_f64);
        // The following will not compile
        // assert_eq!(check_quasigroup_mul_partial(2i8, 3i8, 1i8), -4);
        // assert_eq!(check_quasigroup_mul_partial(2i16, 3i16, 1i16), -4);
        // assert_eq!(check_quasigroup_mul_partial(2i32, 3i32, 1i32), -4);
        // assert_eq!(check_quasigroup_mul_partial(2i64, 3i64, 1i64), -4);
        // assert_eq!(check_quasigroup_add_partial(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_quasigroup_add_partial(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_quasigroup_add_partial(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_quasigroup_add_partial(2u64, 3u64, 1u64), -4);
    }



}
