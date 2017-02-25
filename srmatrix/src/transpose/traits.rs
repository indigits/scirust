
// local imports
use sralgebra::{MagmaBase, CommutativeMonoidAddPartial, CommutativeMonoidMulPartial};
use traits::Shape;

#[doc="Implemented by matrix types
which support transpose operations.
"]
pub trait Transpose<T:MagmaBase> : Shape<T>{

    /// Returns a new matrix holding the transpose
    fn transpose(&self) -> Self;

    // Performs transpose operation within the matrix itself
    //fn transpose_self(&mut self)->&mut Self;
}


#[doc="Frame and Grammian operators for a matrix.
"]
pub trait Frame<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial> : Transpose<T>{

    /// Returns the gram matrix : A' * A
    fn gram(&self) -> Self;

}