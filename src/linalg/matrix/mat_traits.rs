
// std imports
use std::marker::MarkerTrait;



// local imports
use algebra::{Number, Signed};

use error::SRError;


/// Linear algebra methods for matrix of numbers
pub trait LANumberMatrix<T:Number+Signed> {
    /// Returns the determinant of the matrix
    fn det(&self) -> Result<T,SRError>;
}


/// Linear algebra methods for integer matrix
pub trait LAIntMatrix : MarkerTrait{

}

/// Linear algebra methods for signed integer matrix
pub trait LASignedMatrix : MarkerTrait{

}


/// Linear algebra methods for float matrix
pub trait LAFloatMatrix : MarkerTrait{

}

