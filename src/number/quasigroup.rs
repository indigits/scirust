

// std imports
use std::ops::{Sub, Neg, Div};


// local imports
use number::ops::{Recip};
use number::magma::{MagmaAddPartial, MagmaAdd,
    MagmaMulPartial, MagmaMul};

///////////////////////////////////////////////////////////

pub trait QuasiGroupAddPartial
    : MagmaAddPartial
    + Sub<Self, Output=Self>
    + Neg<Output=Self>
{

}


impl<T> QuasiGroupAddPartial for T where
    T: MagmaAddPartial
      + Sub<Self, Output=Self> + Neg<Output=Self>,
{}


///////////////////////////////////////////////////////////


pub trait QuasiGroupAdd 
    : MagmaAdd
    + QuasiGroupAddPartial
{

}


impl<T> QuasiGroupAdd for T where
    T: MagmaAdd + QuasiGroupAddPartial,
{}


///////////////////////////////////////////////////////////

pub trait QuasiGroupMulPartial
    : MagmaMulPartial
    + Div<Self, Output=Self>
    + Recip<Output=Self>
{

}


impl<T> QuasiGroupMulPartial for T where
    T: MagmaMulPartial 
      + Div<Self, Output=Self> + Recip<Output=Self>,
{}

///////////////////////////////////////////////////////////

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
        -e
    }


    #[test]
    fn test_quasigroup_add_partial() {
        assert_eq!(check_quasigroup_add_partial(2i8, 3i8, 1i8), -4);
        assert_eq!(check_quasigroup_add_partial(2i16, 3i16, 1i16), -4);
        assert_eq!(check_quasigroup_add_partial(2i32, 3i32, 1i32), -4);
        assert_eq!(check_quasigroup_add_partial(2i64, 3i64, 1i64), -4);
        assert_eq!(check_quasigroup_add_partial(2f32, 3f32, 1f32), -4f32);
        assert_eq!(check_quasigroup_add_partial(2f64, 3f64, 1f64), -4f64);
        // The following will not compile
        // assert_eq!(check_quasigroup_add_partial(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_quasigroup_add_partial(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_quasigroup_add_partial(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_quasigroup_add_partial(2u64, 3u64, 1u64), -4);
    }


    fn check_quasigroup_add<T: QuasiGroupAdd>(a: T, b: T, c : T)->T{
        let d = a + b;
        let e = d - c;
        -e
    }


    #[test]
    fn test_quasigroup_add() {
        assert_eq!(check_quasigroup_add(2i8, 3i8, 1i8), -4);
        assert_eq!(check_quasigroup_add(2i16, 3i16, 1i16), -4);
        assert_eq!(check_quasigroup_add(2i32, 3i32, 1i32), -4);
        assert_eq!(check_quasigroup_add(2i64, 3i64, 1i64), -4);
        // The following will not compile
        // assert_eq!(check_quasigroup_add(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_quasigroup_add(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_quasigroup_add(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_quasigroup_add(2u64, 3u64, 1u64), -4);
        //assert_eq!(check_quasigroup_add(2f32, 3f32, 1f32), -4f32);
        //assert_eq!(check_quasigroup_add(2f64, 3f64, 1f64), -4f64);
    }



    fn check_quasigroup_mul_partial<T: QuasiGroupMulPartial>(a: T, b: T, c : T)->T{
        let d = a * b;
        let e = d / c;
        e.recip()
    }


    #[test]
    fn test_quasigroup_mul_partial() {
        assert_eq!(check_quasigroup_mul_partial(2f32, 3f32, 3f32), 0.5f32);
        assert_eq!(check_quasigroup_mul_partial(2f64, 3f64, 3f64), 0.5f64);
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
