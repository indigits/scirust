
// std imports
use std::fmt::Debug;
// use std::fmt::Display;
// use std::marker::MarkerTrait;
use std::ops::{Add, Mul};

pub trait MagmaBase : Debug + Clone + Sized {

}

impl<T> MagmaBase for T where
    T : Debug + Clone + Sized
{

}

pub trait MagmaAddPartial 
    : MagmaBase 
    + Add<Self, Output=Self> 
    + PartialEq{

}

impl<T> MagmaAddPartial for T where
    T: MagmaBase + Add<Self, Output=Self> + PartialEq,
{}

///////////////////////////////////////////////////////////


pub trait MagmaAdd 
    : MagmaAddPartial
    + Eq {

}

impl<T> MagmaAdd for T where
    T: MagmaAddPartial + Eq,
{}

///////////////////////////////////////////////////////////


pub trait MagmaMulPartial
    : MagmaBase 
    + Mul<Self, Output=Self> 
    + PartialEq{

}

impl<T> MagmaMulPartial for T where
    T: MagmaBase + Mul<Self, Output=Self> + PartialEq,
{}


///////////////////////////////////////////////////////////


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
