#![doc="Functions for vectors

For many operations below, 
it doesn't matter whether the vector is a column vector
or a row vector. The stride is always one in both cases.
"]


// std imports

// local imports
use number::{One, Zero};
use number::Number;
use matrix::matrix::Matrix;
use matrix::traits::{Shape, MatrixBuffer};

pub struct VecIterator<T:Number> {
    ptr : *const T,
    pos : usize,
    len : usize
}

impl <T:Number> VecIterator<T> {
    pub fn new (ptr: *const T, len : usize)->VecIterator<T>{
        VecIterator{ptr: ptr, pos : 0, len : len}
    } 
}

/// Implementation of iterator trait
impl<T:Number> Iterator for VecIterator<T>{
    type Item = T;
    /// Next element in the vector
    fn next(&mut self)->Option<T> {
        if self.pos == self.len {
            return None;
        }
        let offset = self.pos;
        self.pos += 1;
        Some(unsafe{*self.ptr.offset(offset as isize)})
    }

    /// Returns the upper and lower bound on the remaining length
    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) { 
        let n = self.len  - self.pos;
        (n, Some(n)) 
    }    
} 

/// Constructs a vector iterator for a matrix with the
/// assumption that the matrix is indeed a vector
pub fn vec_iter<T:Number>(m : &Matrix<T>) -> VecIterator<T> {
    assert!(m.is_vector());
    VecIterator::new(m.as_ptr(), m.num_cells())
}

/// Computes the sum of entries in vector v
pub fn vec_reduce_sum<T:Number>(v : &Matrix<T>) -> T{
    let mut result : T = Zero::zero();
    for entry in vec_iter(v){
        result = result + entry;
    }
    result
}


/// Computes the product of entries in vector v
pub fn vec_reduce_product<T:Number>(v : &Matrix<T>) -> T{
    let mut result : T = One::one();
    for entry in vec_iter(v){
        result = result * entry;
    }
    result
}


/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/
#[cfg(test)]
mod test{

    use super::*;
    use matrix::traits::*;
    use matrix::constructors::*;

    #[test]
    fn test_vec_iter(){
        // A column vector
        let v = vector_i64([1, 2, 3, 4].as_slice());
        assert!(v.is_col());
        let mut i = vec_iter(&v);
        assert_eq!(i.next(), Some(1));
        assert_eq!(i.next(), Some(2));
        assert_eq!(i.next(), Some(3));
        assert_eq!(i.next(), Some(4));
        assert_eq!(i.next(), None);
        // A row vector
        let v = v.transpose();
        assert!(v.is_row());
        let mut i = vec_iter(&v);
        assert_eq!(i.next(), Some(1));
        assert_eq!(i.next(), Some(2));
        assert_eq!(i.next(), Some(3));
        assert_eq!(i.next(), Some(4));
        assert_eq!(i.next(), None);
    }

    #[test]
    fn test_vec_reduce_sum(){
        // A column vector
        let v = vector_i64([1, 2, 3, 4].as_slice());
        assert_eq!(vec_reduce_sum(&v), 10);
        // A row vector
        let v = v.transpose();
        assert_eq!(vec_reduce_sum(&v), 10);
    }

    #[test]
    fn test_vec_reduce_prod(){
        // A column vector
        let v = vector_i64([1, 2, 3, 4].as_slice());
        assert_eq!(vec_reduce_product(&v), 24);
        // A row vector
        let v = v.transpose();
        assert_eq!(vec_reduce_product(&v), 24);
    }

}
