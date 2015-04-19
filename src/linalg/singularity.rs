#![doc="Methods for checking singularity or invertibility of matrices.
"]


// std imports
use std::cmp;

// local imports
use algebra::Zero;
use algebra::Number;
//use error::SRError;
use matrix::matrix::Matrix;
use matrix::traits::{Shape, NumberMatrix};

/// Indicates if a lower triangular matrix is singular or not.
pub fn is_singular_lt<T:Number>(m : &Matrix<T>) -> bool {
    if ! m.is_square() {
        return false;
    }
    debug_assert!(m.is_lt());
    has_zero_on_diagonal(m)
}

/// Indicates if a upper triangular matrix is singular or not.
pub fn is_singular_ut<T:Number>(m : &Matrix<T>) -> bool {
    if ! m.is_square() {
        return false;
    }
    debug_assert!(m.is_ut());
    has_zero_on_diagonal(m)
}

/// Indicates if a triangular matrix is singular or not.
pub fn is_singular_triangular<T:Number>(m : &Matrix<T>) -> bool {
    if ! m.is_square() {
        return false;
    }
    debug_assert!(m.is_triangular());
    has_zero_on_diagonal(m)
}

/// Indicates if a diagonal matrix is singular or not.
pub fn is_singular_diagonal<T:Number>(m : &Matrix<T>) -> bool {
    if ! m.is_square() {
        return false;
    }
    debug_assert!(m.is_diagonal());
    has_zero_on_diagonal(m)
}

/// Checks if any entry on the main diagonal is zero
pub fn has_zero_on_diagonal<T:Number>(m : &Matrix<T>) -> bool {
    let n = cmp::min(m.num_rows(), m.num_cols());
    let z : T = Zero::zero();
    for i in range(0, n){
        if m.get(i, i) == z {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod test{

    use matrix::constructors::*;
    use matrix::traits::*;
    use super::*;
    use linalg::matrix::mat_traits::*;

    #[test]
    fn test_triangular_singularity(){
        let  m = from_range_rw_f64(6, 6, 1., 500.);
        println!{"{}", m};
        let mut l = m.lt();
        assert!(!is_singular_lt(&l));
        assert!(l.det().unwrap() != 0.);
        l.set(4, 4, 0.);
        assert!(is_singular_lt(&l));
        assert_eq!(l.det().unwrap(), 0.);

        let mut u = m.ut();
        assert!(!is_singular_ut(&u));
        assert!(u.det().unwrap() != 0.);
        u.set(4, 4, 0.);
        assert!(is_singular_ut(&u));
        assert_eq!(u.det().unwrap(), 0.);


        assert!(is_singular_triangular(&l));
        assert!(is_singular_triangular(&u));

        let mut d = m.diagonal_matrix();
        assert!(!is_singular_diagonal(&d));
        assert!(d.det().unwrap() != 0.);
        d.set(3, 3, 0.);
        assert!(is_singular_diagonal(&d));
        assert_eq!(d.det().unwrap(), 0.);
    }

}