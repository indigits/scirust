#![doc="Routines related to permutations
"]

// local imports
use number::Number;
use matrix::{Matrix, MatrixU8};
use matrix::traits::MatrixType;

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
    println!("flags {}", flags);
    for i in range(0, n){
        if flags[i] == 0 {
            return false;
        }
    }
    true
}



/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/

#[cfg(test)]
mod test{
    use super::*;
    use matrix::*;

    #[test]
    fn test_permutation(){
        let v = vector_i64([0, 1, 2]);
        assert!(is_permutation(&v));
        let v = vector_i64([2, 1, 0]);
        assert!(is_permutation(&v));
        let v = vector_i64([3, 1, 0]);
        assert!(!is_permutation(&v));
        let v = vector_i64([2, 0, 1]);
        assert!(is_permutation(&v));
        let v = vector_i64([2, 1, 1]);
        assert!(!is_permutation(&v));
    }
}

