

// std imports

// local imports
use number::Number;
use super::eo_traits::{ERO, ECO};
use matrix::matrix::Matrix;

/// Implementation of Elementary row operations.
impl<T:Number> ERO<T> for Matrix<T> {
}

/// Implementation of Elementary column operations.
impl<T:Number> ECO<T> for Matrix<T> {
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


