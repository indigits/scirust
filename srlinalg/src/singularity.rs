#![doc="Methods for checking singularity or invertibility of matrices.
"]


// std imports

// external imports

// local imports
use sralgebra::{CommutativeMonoidAddPartial,
CommutativeMonoidMulPartial};
use srmatrix::api::*;

/// Indicates if a lower triangular matrix is singular or not.
pub fn is_singular_lt<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial>(m : &Matrix<T>) -> bool {
    if ! m.is_square() {
        return false;
    }
    debug_assert!(m.is_lt());
    has_zero_on_diagonal(m)
}

/// Indicates if a upper triangular matrix is singular or not.
pub fn is_singular_ut<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial>(m : &Matrix<T>) -> bool {
    if ! m.is_square() {
        return false;
    }
    debug_assert!(m.is_ut());
    has_zero_on_diagonal(m)
}

/// Indicates if a triangular matrix is singular or not.
pub fn is_singular_triangular<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial>(m : &Matrix<T>) -> bool {
    if ! m.is_square() {
        return false;
    }
    debug_assert!(m.is_triangular());
    has_zero_on_diagonal(m)
}

/// Indicates if a diagonal matrix is singular or not.
pub fn is_singular_diagonal<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial>(m : &Matrix<T>) -> bool {
    if ! m.is_square() {
        return false;
    }
    debug_assert!(m.is_diagonal());
    has_zero_on_diagonal(m)
}

/// Checks if any entry on the main diagonal is zero
pub fn has_zero_on_diagonal<T:CommutativeMonoidAddPartial>(m : &Matrix<T>) -> bool {
    m.diagonal_iter().any(|x| x == T::zero())
}

#[cfg(test)]
mod test{

    use srmatrix::api::*;
    use super::*;
    use matrix::mat_traits::*;

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
