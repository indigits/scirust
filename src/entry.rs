#![doc="Defines trait for an entry in a matrix. 
"]

// std imports
use std::fmt::Show;

// complex numbers
pub use external::complex::{Complex, Complex32, Complex64};


// local imports

/// Defines an additive identity element for Self
pub trait Zero {
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
}


/// Defines a multiplicative identity element for Self
pub trait One {
    fn one() -> Self;
}

/// Defines basic requirements for a matrix entry
pub trait Entry : Show + Clone + Zero {

}

/******************************************************
 *
 *   Zero implementations.
 *
 *******************************************************/

impl Zero for i8{

    #[inline]
    fn zero() -> i8{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for i16{
    #[inline]
    fn zero() -> i16{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for i32{
    #[inline]
    fn zero() -> i32{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for i64{
    #[inline]
    fn zero() -> i64{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for int{
    #[inline]
    fn zero() -> int{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for u8{
    #[inline]
    fn zero() -> u8{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for u16{
    #[inline]
    fn zero() -> u16{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for u32{
    #[inline]
    fn zero() -> u32{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for u64{
    #[inline]
    fn zero() -> u64{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for uint{
    #[inline]
    fn zero() -> uint{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for f32{
    #[inline]
    fn zero() -> f32{
        0.
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for f64{
    #[inline]
    fn zero() -> f64{
        0.
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for Complex32{
    #[inline]
    fn zero() -> Complex32{
        let v : Complex32 = Complex::zero();
        v
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

impl Zero for Complex64{
    #[inline]
    fn zero() -> Complex64{
        Complex::zero()
    }

    #[inline]
    fn is_zero(&self) -> bool{
        self.is_zero()
    }
}

/******************************************************
 *
 *   One implementations.
 *
 *******************************************************/



impl One for i8{

    #[inline]
    fn one() -> i8 {
        1
    }
}

impl One for i16{

    #[inline]
    fn one() -> i16 {
        1
    }
}

impl One for i32{

    #[inline]
    fn one() -> i32 {
        1
    }
}

impl One for i64{

    #[inline]
    fn one() -> i64 {
        1
    }
}

impl One for int{

    #[inline]
    fn one() -> int {
        1
    }
}

impl One for u8{

    #[inline]
    fn one() -> u8 {
        1
    }
}

impl One for u16{

    #[inline]
    fn one() -> u16 {
        1
    }
}

impl One for u32{

    #[inline]
    fn one() -> u32 {
        1
    }
}

impl One for u64{

    #[inline]
    fn one() -> u64 {
        1
    }
}

impl One for uint{

    #[inline]
    fn one() -> uint {
        1
    }
}

impl One for f32{

    #[inline]
    fn one() -> f32 {
        1.
    }
}

impl One for f64{

    #[inline]
    fn one() -> f64 {
        1.
    }
}

impl One for Complex32{

    #[inline]
    fn one() -> Complex32 {
        Complex::one()
    }
}

impl One for Complex64{

    #[inline]
    fn one() -> Complex64 {
        Complex::one()
    }
}


/******************************************************
 *
 *   Entry implementations.
 *
 *******************************************************/



impl Entry for i8{

}

impl Entry for i16{

}

impl Entry for i32{

}

impl Entry for i64{

}

impl Entry for int{

}

impl Entry for u8{

}

impl Entry for u16{

}

impl Entry for u32{

}

impl Entry for u64{

}

impl Entry for uint{

}

impl Entry for f32{

}

impl Entry for f64{

}

impl Entry for Complex32{

}

impl Entry for Complex64{

}




