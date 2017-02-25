
// std imports



// local imports
use sralgebra::{CommutativeRingPartial};

use error::SRError;


/// Linear algebra methods for matrix of numbers
pub trait LANumberMatrix<T:CommutativeRingPartial> {
    /// Returns the determinant of the matrix
    fn det(&self) -> Result<T,SRError>;
}


/// Linear algebra methods for integer matrix
pub trait LAIntMatrix {

}

/// Linear algebra methods for signed integer matrix
pub trait LASignedMatrix {

}


/// Linear algebra methods for float matrix
pub trait LAFloatMatrix {

}

