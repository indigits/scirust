#![doc="Defines some operators.
"]

// std imports
use std::ops::Div;

/// A type for which computing the reciprocal is supported.
/// Supported types include f32, f64, etc.. 
pub trait Recip {
    /// The resulting type after computing 1 / x.
    type Output = Self;

    /// The method for the computing 1 / x.
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


/// Marker interface for restricting division to only those types which
/// support exact division
pub trait Division : Div
{

}

impl Division for f32 {}
impl Division for f64 {}
