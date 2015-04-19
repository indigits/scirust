#![doc="Defines generic traits for signed  numbers.
"]


// std imports
use std::ops::{Neg};
use std::num::{SignedInt, Float};

// local imports
use algebra::number::Number;

pub trait Signed: Number + Neg<Output=Self> {
    fn abs_val(&self) -> Self;
    fn signum(&self) -> Self;
    fn is_positive(&self) -> bool;
    fn is_negative(&self) -> bool;
}


/******************************************************
 *      Signed implementation
 ******************************************************/


impl Signed for i8 {

    #[inline]
    fn abs_val(&self) -> i8{
        self.abs()
    }

    #[inline]
    fn signum(&self) -> i8{
        SignedInt::signum(*self)
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        SignedInt::is_positive(*self)
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        SignedInt::is_negative(*self)
    }
}


impl Signed for i16 {

    #[inline]
    fn abs_val(&self) -> i16{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i16{
        SignedInt::signum(*self)
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        SignedInt::is_positive(*self)
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        SignedInt::is_negative(*self)
    }
}

impl Signed for i32 {

    #[inline]
    fn abs_val(&self) -> i32{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i32{
        SignedInt::signum(*self)
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        SignedInt::is_positive(*self)
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        SignedInt::is_negative(*self)
    }
}
    
impl Signed for i64 {

    #[inline]
    fn abs_val(&self) -> i64{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i64{
        SignedInt::signum(*self)
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        SignedInt::is_positive(*self)
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        SignedInt::is_negative(*self)
    }
}  

impl Signed for f32 {

    #[inline]
    fn abs_val(&self) -> f32{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> f32{
        Float::signum(*self)
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        Float::is_positive(*self)
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        Float::is_negative(*self)
    }
}

impl Signed for f64 {

    #[inline]
    fn abs_val(&self) -> f64{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> f64{
        Float::signum(*self)
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        Float::is_positive(*self)
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        Float::is_negative(*self)
    }
}




