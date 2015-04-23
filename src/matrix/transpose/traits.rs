
// local imports
use algebra::structure::{MagmaBase, FieldPartial};
use matrix::traits::Shape;

#[doc="Implemented by matrix types
which support transpose operations.
"]
pub trait Transpose<T:FieldPartial> : Shape<T>{

    /// Returns a new matrix holding the transpose
    fn transpose(&self) -> Self;

    /// Returns the gram matrix : A' * A
    fn gram(&self) -> Self;

    // Performs transpose operation within the matrix itself
    //fn transpose_self(&mut self)->&mut Self;
}


