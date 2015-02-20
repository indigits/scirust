#![doc="Defines generic traits for signed  numbers.
"]


// std imports
use std::num::{SignedInt, Float};

// local imports
use number::number::Number;

pub trait Signed: Number + Neg<Self> {
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
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}


impl Signed for i16 {

    #[inline]
    fn abs_val(&self) -> i16{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i16{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}

impl Signed for i32 {

    #[inline]
    fn abs_val(&self) -> i32{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i32{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}
    
impl Signed for i64 {

    #[inline]
    fn abs_val(&self) -> i64{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i64{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}
   

impl Signed for isize {

    #[inline]
    fn abs_val(&self) -> isize{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> isize{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}
    

impl Signed for f32 {

    #[inline]
    fn abs_val(&self) -> f32{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> f32{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}

impl Signed for f64 {

    #[inline]
    fn abs_val(&self) -> f64{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> f64{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}




