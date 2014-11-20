#![doc="Statistical moments
"]

// local imports
use entry::Zero;
use number::Number;
use matrix::matrix::Matrix;
use matrix::traits::{Shape, MatrixBuffer, Strided, 
 StridedNumberMatrix};


/// Computes sum over columns and returns a row vector
pub fn sum_cw<T:Number>(m : & StridedNumberMatrix<T>) -> Matrix<T> {
    let cols = m.num_cols();
    let rows = m.num_rows();
    let mut result = Matrix::new(1, cols);
    let ptr = m.as_ptr();
    let stride = m.stride() as int;
    let mut offset = m.start_offset();
    for c in range(0, cols) {
        let mut sum : T = Zero::zero(); 
        for r in range(0, rows){
            sum = sum + unsafe{*ptr.offset(offset + r as int)};
        }
        offset += stride;
        result.set(0, c, sum);
    }
    result
}

/// Computes sum over rows and returns a column vector
pub fn sum_rw<T:Number>(m : & StridedNumberMatrix<T>) -> Matrix<T> {
    let cols = m.num_cols();
    let rows = m.num_rows();
    let mut result = Matrix::new(rows, 1);
    let ptr = m.as_ptr();
    let stride = m.stride() as int;
    let mut offset = m.start_offset();
    for r in range(0, rows) {
        let mut sum : T = Zero::zero();
        let mut src_offset  = offset; 
        for _ in range(0, cols){
            sum = sum + unsafe{*ptr.offset(src_offset)};
            src_offset += stride;
        }
        offset += 1;
        result.set(r, 0, sum);
    }
    result
}

/// Computes sum over columns and returns a row vector
pub fn sum_sqr_cw<T:Number>(m : & StridedNumberMatrix<T>) -> Matrix<T> {
    let cols = m.num_cols();
    let rows = m.num_rows();
    let mut result = Matrix::new(1, cols);
    let ptr = m.as_ptr();
    let stride = m.stride() as int;
    let mut offset = m.start_offset();
    for c in range(0, cols) {
        let mut sum : T = Zero::zero(); 
        for r in range(0, rows){
            let v  = unsafe{*ptr.offset(offset + r as int)};
            sum = sum + v * v;
        }
        offset += stride;
        result.set(0, c, sum);
    }
    result
}

/// Computes sum over rows and returns a column vector
pub fn sum_sqr_rw<T:Number>(m : & StridedNumberMatrix<T>) -> Matrix<T> {
    let cols = m.num_cols();
    let rows = m.num_rows();
    let mut result = Matrix::new(rows, 1);
    let ptr = m.as_ptr();
    let stride = m.stride() as int;
    let mut offset = m.start_offset();
    for r in range(0, rows) {
        let mut sum : T = Zero::zero();
        let mut src_offset  = offset; 
        for _ in range(0, cols){
            let v = unsafe{*ptr.offset(src_offset)};
            sum = sum + v * v;
            src_offset += stride;
        }
        offset += 1;
        result.set(r, 0, sum);
    }
    result
}


/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/
#[cfg(test)]
mod test{

    use super::*;
    use api::*;


    #[test]
    fn test_moment_sum_cw_1(){
        let m = matrix_rw_i32(3, 3, &[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9]);
        let s = sum_cw(&m);
        assert_eq!(s, matrix_cw_i32(1,3, &[12, 15, 18]));
    }


    #[test]
    fn test_moment_sum_cw_2(){
        let m = matrix_rw_i32(3, 3, &[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9]);
        let v = m.view(1, 1, 2, 2);
        let s = sum_cw(&v);
        assert_eq!(s, matrix_cw_i32(1,2, &[13, 15]));
    }

    #[test]
    fn test_moment_sum_rw_1(){
        let m = matrix_rw_i32(3, 3, &[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9]);
        let s = sum_rw(&m);
        assert_eq!(s, matrix_rw_i32(3,1, &[6, 15, 24]));
    }

    #[test]
    fn test_moment_sum_rw_2(){
        let m = matrix_rw_i32(3, 3, &[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9]);
        let v = m.view(1, 1, 2, 2);
        let s = sum_rw(&v);
        assert_eq!(s, matrix_rw_i32(2,1, &[11, 17]));
    }

    #[test]
    fn test_moment_sum_sqr_cw_1(){
        let m = matrix_rw_i32(3, 3, &[
            1, 1, 2,
            2, 2, 1,
            3, 2, 2]);
        let s = sum_sqr_cw(&m);
        assert_eq!(s, matrix_cw_i32(1,3, &[14, 9, 9]));
    }

    #[test]
    fn test_moment_sum_sqr_cw_2(){
        let m = matrix_rw_i32(3, 3, &[
            1, 1, 2,
            2, 2, 1,
            3, 2, 2]);
        let v = m.view(1, 1, 2, 2);
        let s = sum_sqr_cw(&v);
        assert_eq!(s, matrix_cw_i32(1,2, &[8, 5]));
    }

    #[test]
    fn test_moment_sum_sqr_rw_1(){
        let m = matrix_rw_i32(3, 3, &[
            1, 1, 2,
            2, 2, 1,
            3, 2, 2]);
        let s = sum_sqr_rw(&m);
        assert_eq!(s, matrix_cw_i32(3,1, &[6, 9, 17]));
    }

    #[test]
    fn test_moment_sum_sqr_rw_2(){
        let m = matrix_rw_i32(3, 3, &[
            1, 1, 2,
            2, 2, 1,
            3, 2, 2]);
        let v = m.view(1, 1, 2, 2);
        let s = sum_sqr_rw(&v);
        assert_eq!(s, matrix_cw_i32(2,1, &[5, 8]));
    }
}