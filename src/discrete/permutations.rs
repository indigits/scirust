#![doc="Routines related to permutations
"]

// std imports
use num::{FromPrimitive, ToPrimitive};

// local imports
use algebra::structure::CommutativeMonoidAddPartial;
use matrix::matrix::{Matrix, MatrixU8};
use matrix::traits::Shape;

/// Tells whether a vector is a permutation or not.
pub fn is_permutation<T:CommutativeMonoidAddPartial+ToPrimitive>(vector : &Matrix<T>)-> bool{
    assert!(vector.is_vector());
    let n: usize = vector.num_cells();
    let mut flags : MatrixU8 = Matrix::zeros(n, 1);
    for i in 0..n{
        let v = vector[i];
        let v2 : u64 = v.to_u64().unwrap();
        let v3 = v2 as usize;
        if v2 >= (n  as u64) {
            return false;
        }
        flags.set(v3, 0, 1);
    }
    for v in flags.cell_iter(){
        if v == 0u8 {
            return false;
        }
    }
    true
}

/// Finds the inverse permutation
pub fn inverse_permutation<T:CommutativeMonoidAddPartial+ToPrimitive+FromPrimitive>(vector : &Matrix<T>)-> Matrix<T>{
    assert!(vector.is_vector());
    debug_assert!(is_permutation(vector));
    let n = vector.num_cells();
    let mut result : Matrix<T> = Matrix::zeros(n, 1);
    // Reverse the order
    for i in 0..n{
        let index = vector[i];
        let index = index.to_usize().unwrap();
        let value : T = FromPrimitive::from_usize(i).unwrap();
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

