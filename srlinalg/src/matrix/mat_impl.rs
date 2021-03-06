
// std imports

// external imports
use num::Signed;

// local imports
use sralgebra::{CommutativeRingPartial};
use super::mat_traits::LANumberMatrix;
use srmatrix::api::*;
use det;

impl<T:CommutativeRingPartial+Signed> LANumberMatrix<T> for Matrix<T>{
    /// Returns determinant of the matrix
    fn det(&self) -> Result<T,SRError>{
        det::det(self)
    }
}




/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/


#[cfg(test)]
mod test{
    //use super::*;
}


/******************************************************
 *
 *   Bench marks
 *
 *******************************************************/


#[cfg(test)]
mod bench{
    //extern crate test;
    //use self::test::Bencher;
    //use super::*;
}


