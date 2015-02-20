


// local imports
use number::{Number, Signed};

use error::SRError;


/// Linear algebra methods for matrix of numbers
pub trait LANumberMatrix<T:Number+Signed> {
    /// Returns the determinant of the matrix
    fn det(&self) -> Result<T,SRError>;
}


/// Linear algebra methods for isizeeger matrix
pub trait LAIntMatrix{

}

/// Linear algebra methods for signed isizeeger matrix
pub trait LASignedMatrix{

}


/// Linear algebra methods for float matrix
pub trait LAFloatMatrix{

}

