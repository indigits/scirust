#![doc="Routines related to permutations
"]

// std imports
use std::num;
use std::num::Int;

// local imports
use number::Number;
use matrix::matrix::{Matrix, MatrixU8};
use matrix::traits::Shape;

/// Tells whether a vector is a permutation or not.
pub fn is_permutation<T:Number+Int>(vector : &Matrix<T>)-> bool{
    assert!(vector.is_vector());
    let n = vector.num_cells();
    let mut flags : MatrixU8 = Matrix::zeros(n, 1);
    for i in range(0, n){
        let v = vector[i];
        let v = v.to_uint().unwrap();
        if v >= n{
            return false;
        }
        flags.set(v, 0, 1);
    }
    for i in range(0, n){
        if flags[i] == 0 {
            return false;
        }
    }
    true
}

/// Finds the inverse permutation
pub fn inverse_permutation<T:Number+Int+FromPrimitive>(vector : &Matrix<T>)-> Matrix<T>{
    assert!(vector.is_vector());
    debug_assert!(is_permutation(vector));
    let n = vector.num_cells();
    let mut result : Matrix<T> = Matrix::zeros(n, 1);
    // Reverse the order
    for i in range(0, n){
        let index = vector[i];
        let index = index.to_uint().unwrap();
        let value : T = num::from_uint(i).unwrap();
        result.set(index, 0, value);
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
    use matrix::constructors::*;

    #[test]
    fn test_permutation(){
        let v = vector_i64([0, 1, 2].as_slice());
        assert!(is_permutation(&v));
        let v = vector_i64([2, 1, 0].as_slice());
        assert!(is_permutation(&v));
        let v = vector_i64([3, 1, 0].as_slice());
        assert!(!is_permutation(&v));
        let v = vector_i64([2, 0, 1].as_slice());
        assert!(is_permutation(&v));
        let v = vector_i64([2, 1, 1].as_slice());
        assert!(!is_permutation(&v));
    }

    #[test]
    fn test_inverse_permutation(){
        let v = vector_i64([0, 3, 2, 1].as_slice());
        let v2 = inverse_permutation(&v);
        assert_eq!(v2, vector_i64([0, 3, 2, 1].as_slice()));
        let v = vector_i64([0, 2, 3, 1, 4].as_slice());
        let v2 = inverse_permutation(&v);
        assert_eq!(v2, vector_i64([0, 3, 1, 2, 4].as_slice()));
    }
}

