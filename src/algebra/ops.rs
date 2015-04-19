#![doc="Defines some operators.
"]



pub trait Recip {
    /// The resulting type after applying the `-` operator
    type Output;

    /// The method for the unary `-` operator
    fn recip(self) -> Self::Output;
}


impl Recip for f32 {
    type Output = Self;
    #[inline]
    fn recip(self) -> f32{
        1.0 / self
    }    
}

impl Recip for f64 {
    type Output = Self;
    #[inline]
    fn recip(self) -> f64{
        1.0 / self
    }    
}

