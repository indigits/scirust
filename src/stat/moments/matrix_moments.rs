#![doc="Implementation of statistical moments for matrices
"]

// std imports
use std::num::Float;

// local imports
use number::{Zero, One};
use number::Number;
use matrix::matrix::Matrix;
use matrix::traits::{Shape, 
    MatrixBuffer, Strided,
    CopyUpdates, InPlaceUpdates, Transpose, 
    ERO, ECO};
use stat::moments::traits::{Sums, Moments};

impl <T:Number> Sums<T> for Matrix<T> {

    /// Computes sum over columns and returns a row vector
    fn sum_cw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(1, cols);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
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

    fn sum_rw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(rows, 1);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
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
    fn sum_sqr_cw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(1, cols);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
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
    fn sum_sqr_rw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(rows, 1);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
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

}

impl <T:Number + Float + FromPrimitive> Moments<T> for Matrix<T> {

    fn mean_cw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let rows_t : T = FromPrimitive::from_uint(rows).unwrap();
        let mut result = Matrix::new(1, cols);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for c in range(0, cols) {
            let mut sum : T = Zero::zero(); 
            for r in range(0, rows){
                sum = sum + unsafe{*ptr.offset(offset + r as int)};
            }
            offset += stride;
            result.set(0, c, sum / rows_t);
        }
        result
    }

    /// Computes mean over rows and returns a column vector
    fn mean_rw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let cols_t : T = FromPrimitive::from_uint(cols).unwrap();
        let mut result = Matrix::new(rows, 1);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for r in range(0, rows) {
            let mut sum : T = Zero::zero();
            let mut src_offset  = offset; 
            for _ in range(0, cols){
                sum = sum + unsafe{*ptr.offset(src_offset)};
                src_offset += stride;
            }
            offset += 1;
            result.set(r, 0, sum / cols_t);
        }
        result
    }

    fn mean_sqr_cw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(1, cols);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for c in range(0, cols) {
            let mut sum : T = Zero::zero(); 
            for r in range(0, rows){
                let v  = unsafe{*ptr.offset(offset + r as int)};
                sum = sum + v * v;
            }
            offset += stride;
            result.set(0, c, sum);
        }
        let rows_t : T = FromPrimitive::from_uint(rows).unwrap();
        result.ero_scale(0, rows_t.powi(-1));
        result
    }

    /// Computes mean square over rows and returns a column vector
    fn mean_sqr_rw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(rows, 1);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for r in range(0, rows) {
            let mut sum : T = Zero::zero();
            let mut src_offset  = offset; 
            for _ in range(0, cols){
                let v = unsafe{*ptr.offset(src_offset)};
                // sum = sum + sum + v * v;
                sum = v.mul_add(v, sum);
                src_offset += stride;
            }
            offset += 1;
            result.set(r, 0, sum);
        }
        let cols_t : T = FromPrimitive::from_uint(cols).unwrap();
        result.eco_scale(0, cols_t.powi(-1));
        result
    }

    /// Computes sample variance over columns and returns a row vector
    /// sum((m - repmat(mean(m), r, 1)).^2 ) / (r - 1)
    fn var_cw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut sum_vec = Matrix::new(1, cols);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for c in range(0, cols) {
            let mut sum : T = Zero::zero(); 
            for r in range(0, rows){
                let v  = unsafe{*ptr.offset(offset + r as int)};
                sum = sum + v;
            }
            offset += stride;
            sum_vec.set(0, c, sum);
        }
        let rows_t : T = FromPrimitive::from_uint(self.num_rows()).unwrap();
        // get the mean.
        sum_vec.ero_scale(0, rows_t.powi(-1));
        // now subtract and square.
        let mut sum_sqr_vec = Matrix::new(1, cols);
        offset = self.start_offset();
        for c in range(0, cols){
            let mut sum_sqr : T = Zero::zero(); 
            let mean = sum_vec.get(0, c);
            for r in range(0, rows){
                let v  = unsafe{*ptr.offset(offset + r as int)} - mean;
                sum_sqr = sum_sqr + v * v;
            }
            offset += stride;
            sum_sqr_vec.set(0, c, sum_sqr);
        }
        let denom = rows_t - One::one();
        sum_sqr_vec.ero_scale(0, denom.powi(-1));
        sum_sqr_vec
    }


    /// Computes sample variance over rows and returns a column vector
    fn var_rw(&self) -> Matrix<T> {
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(rows, 1);
        let ptr = self.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        // convert to type T
        let cols_t : T = FromPrimitive::from_uint(cols).unwrap();
        for r in range(0, rows) {
            let mut sum : T = Zero::zero();
            let mut src_offset  = offset; 
            for _ in range(0, cols){
                sum = sum + unsafe{*ptr.offset(src_offset)};
                src_offset += stride;
            }
            offset += 1;
            result.set(r, 0, sum / cols_t);
        }
        // now subtract and square.
        let mut var_vec = Matrix::new(rows, 1);
        offset = self.start_offset();
        let denom = cols_t - One::one();
        for r in range(0, rows) {
            let mut sum_sqr : T = Zero::zero(); 
            let mut src_offset  = offset; 
            let mean = result.get(r, 0);
            for _ in range(0, cols){
                let v  = unsafe{*ptr.offset(src_offset)} - mean;
                sum_sqr = sum_sqr + v * v;
                src_offset += stride;
            }
            offset += 1;
            var_vec.set(r, 0, sum_sqr / denom);
        }
        var_vec
    }



#[doc="Computes the sample covariance from a set of example vectors.
Each column in x is a sample vector from the population.

If x is a vector, then we return the variance of the vector.
"]
    fn cov(&self) -> Matrix<T> {
        if self.is_row(){
            return self.var_rw();
        }
        if self.is_col(){
            return self.var_cw();
        }
        // We have got a full matrix of rows of observations
        // The dimension of covariance matrix
        // compute the mean over observation rows.
        let mean = self.mean_cw();
        // subtract the mean from each row
        let y = self.copy_sub_vec_from_rows(&mean).unwrap();
        // Compute the gram matrix of y
        let mut g = y.gram();
        // scale the entries in g
        let rows_t : T = FromPrimitive::from_uint(self.num_rows()).unwrap();
        assert!(g.div_scalar(rows_t - One::one()).is_ok());
        g
    }


}











/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/
#[cfg(test)]
mod test{

    use api::*;
    use stat::moments::traits::Sums;
    use stat::moments::traits::Moments;


    #[test]
    fn test_moment_sum_cw_1(){
        let m = matrix_rw_i32(3, 3, &[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9]);
        let s = m.sum_cw();
        assert_eq!(s, matrix_cw_i32(1,3, &[12, 15, 18]));
    }

    #[test]
    fn test_moment_sum_rw_1(){
        let m = matrix_rw_i32(3, 3, &[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9]);
        let s = m.sum_rw();
        assert_eq!(s, matrix_rw_i32(3,1, &[6, 15, 24]));
    }

    #[test]
    fn test_moment_sum_sqr_cw_1(){
        let m = matrix_rw_i32(3, 3, &[
            1, 1, 2,
            2, 2, 1,
            3, 2, 2]);
        let s = m.sum_sqr_cw();
        assert_eq!(s, matrix_cw_i32(1,3, &[14, 9, 9]));
    }

    #[test]
    fn test_moment_sum_sqr_rw_1(){
        let m = matrix_rw_i32(3, 3, &[
            1, 1, 2,
            2, 2, 1,
            3, 2, 2]);
        let s = m.sum_sqr_rw();
        assert_eq!(s, matrix_cw_i32(3,1, &[6, 9, 17]));
    }

    #[test]
    fn test_moment_mean_cw_1(){
        let m = matrix_rw_f32(3, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            7., 8., 9.]);
        let s = m.mean_cw();
        assert_eq!(s, matrix_cw_f32(1,3, &[12./3., 15./3., 18./3.]));
    }


    #[test]
    fn test_moment_mean_rw_1(){
        let m = matrix_rw_f32(3, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            7., 8., 9.]);
        let s = m.mean_rw();
        assert_eq!(s, matrix_rw_f32(3,1, &[6./3., 15./3., 24./3.]));
    }

    #[test]
    fn test_moment_mean_sqr_cw_1(){
        let m = matrix_rw_f32(3, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            7., 8., 9.]);
        let s = m.mean_sqr_cw();
        assert_eq!(s, matrix_cw_f32(1,3, &[22., 31., 42.]));
    }

    #[test]
    fn test_moment_mean_sqr_rw_1(){
        let m = matrix_rw_f32(4, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            4., 5., 6.,
            7., 8., 9.]);
        let s = m.mean_sqr_rw();
        let d = s - matrix_cw_f32(4,1, &[14./3., 
            77. / 3., 
            77. / 3., 
            194.0 / 3.0
            ]);
        println!("{:e}", d.max_abs_scalar_value());
        assert!(d.max_abs_scalar_value() < 1e-5);
        // for 64-bit floating point, we can be more accurate.
        //assert!(d.max_abs_scalar_value() < 1e-13);
    }

    #[test]
    fn test_moment_var_cw_1(){
        let m = matrix_rw_f32(4, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            4., 5., 6.,
            7., 8., 9.]);
        let s = m.var_cw();
        let e = matrix_cw_f32(1,3, &[6., 
            6. , 
            6.
            ]);
        let d = s - e;
        println!("{}", s);
        println!("{:e}", d.max_abs_scalar_value());
        assert!(d.max_abs_scalar_value() < 1e-12);
        // for 64-bit floating point, we can be more accurate.
        //assert!(d.max_abs_scalar_value() < 1e-13);
    }

    #[test]
    fn test_moment_var_rw_1(){
        let m = matrix_rw_f32(4, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            4., 5., 6.,
            7., 8., 9.]);
        let s = m.var_rw();
        let e = matrix_cw_f32(4,1, &[1., 1., 1., 1.]);
        println!("{}", s);
        let d = s - e;
        println!("{:e}", d.max_abs_scalar_value());
        assert!(d.max_abs_scalar_value() < 1e-12);
        // for 64-bit floating point, we can be more accurate.
        //assert!(d.max_abs_scalar_value() < 1e-13);
    }

    #[test]
    fn test_moment_cov_1(){
        let m = matrix_rw_f32(4, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            4., 5., 6.,
            7., 8., 9.]);
        let s = m.cov();
        let e = matrix_cw_f32(3, 3, &[
            6., 6., 6., 
            6., 6., 6.,
            6., 6., 6.]);
        println!("{}", s);
        let d = s - e;
        println!("{:e}", d.max_abs_scalar_value());
        assert!(d.max_abs_scalar_value() < 1e-12);
        // for 64-bit floating point, we can be more accurate.
        //assert!(d.max_abs_scalar_value() < 1e-13);
    }
}