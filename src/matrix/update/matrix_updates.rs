/// std imports

/// local imports
use number::Number;
use matrix::traits::{Shape, MatrixBuffer, Strided};
use matrix::update::traits::Updates;
use matrix::eo::eo_traits::{ERO, ECO};
use matrix::matrix::Matrix;

/// Implementation of Matrix general update operations.
impl<T:Number> Updates<T> for Matrix<T> {
    fn scale_row_lt(&mut self, r :  uint, scale_factor : T)-> &mut Matrix<T>{
        debug_assert!(r < self.num_rows());
        let stride = self.stride() as int;
        let ptr = self.as_mut_ptr();
        let mut offset = self.cell_to_offset(r, 0);
        for _ in range(0, r + 1){
            unsafe{
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale_factor * v;
                offset += stride;
                debug_assert!(offset < self.capacity() as int);
            }
        }
        self
    }
    fn scale_col_lt(&mut self, c :  uint, scale_factor : T)-> &mut Matrix<T>{
        debug_assert!(c < self.num_cols());
        let ptr = self.as_mut_ptr();
        let mut offset = self.cell_to_offset(c, c);
        for _ in range(c, self.num_cols()){
            unsafe{
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale_factor * v;
                offset += 1;
            }
        }
        self
    }
    fn scale_row_ut(&mut self, r :  uint, scale_factor : T)-> &mut Matrix<T>{
        debug_assert!(r < self.num_rows());
        let stride = self.stride() as int;
        let ptr = self.as_mut_ptr();
        let mut offset = self.cell_to_offset(r, r);
        for _ in range(r, self.num_cols()){
            unsafe{
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale_factor * v;
                offset += stride;
            }
        }
        self
    }
    fn scale_col_ut(&mut self, c :  uint, scale_factor : T)-> &mut Matrix<T>{
        debug_assert!(c < self.num_cols());
        let ptr = self.as_mut_ptr();
        let mut offset = self.cell_to_offset(0, c);
        for _ in range(0, c + 1){
            unsafe{
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale_factor * v;
                offset += 1;
            }
        }
        self
    }
    fn scale_rows(&mut self, scale_factors : &Matrix<T>)-> &mut Matrix<T>{
        assert!(scale_factors.is_col());
        assert_eq!(scale_factors.num_cells(), self.num_rows());
        for r in range(0, self.num_rows()){
            let factor = scale_factors[r];
            self.ero_scale(r, factor);
        }
        self
    }
    fn scale_cols(&mut self, scale_factors : &Matrix<T>)-> &mut Matrix<T>{
        assert!(scale_factors.is_col());
        assert_eq!(scale_factors.num_cells(), self.num_cols());
        for c in range(0, self.num_cols()){
            let factor = scale_factors[c];
            self.eco_scale(c, factor);
        }
        self
    }
}


/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/


#[cfg(test)]
mod test{

    use matrix::constructors::*;
    use matrix::update::traits::Updates;

    #[test]
    fn test_set_diag(){
        let mut m = from_range_rw_f64(4, 4, 1., 100.);
        m.set_diagonal(0.);
        let m2 = matrix_rw_f64(4, 4, [
            0., 2., 3., 4., 
            5., 0., 7., 8.,
            9., 10., 0., 12.,
            13., 14., 15., 0.,
            ].as_slice());
        assert_eq!(m, m2);
        m.set_row(0, 1.);
        let m2 = matrix_rw_f64(4, 4, [
            1., 1., 1., 1., 
            5., 0., 7., 8.,
            9., 10., 0., 12.,
            13., 14., 15., 0.,
            ].as_slice());
        assert_eq!(m, m2);
        m.set_col(2, 20.);
        let m2 = matrix_rw_f64(4, 4, [
            1., 1., 20., 1., 
            5., 0., 20., 8.,
            9., 10., 20., 12.,
            13., 14., 20., 0.,
            ].as_slice());
        assert_eq!(m, m2);
        m.set_block(2, 2, 2, 2, 30.);
        let m2 = matrix_rw_f64(4, 4, [
            1., 1., 20., 1., 
            5., 0., 20., 8.,
            9., 10., 30., 30.,
            13., 14., 30., 30.,
            ].as_slice());
        assert_eq!(m, m2);
    }

    #[test]
    fn test_scale_row_lt(){
        let mut m = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9
            ].as_slice());
        let m2 = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 5, 6,
            21, 24, 27
            ].as_slice());
        m.scale_row_lt(2, 3);
        assert_eq!(m, m2);
    }
    #[test]
    fn test_scale_col_lt(){
        let mut m = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9
            ].as_slice());
        let m2 = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 5, 6,
            7, 8, 27
            ].as_slice());
        m.scale_col_lt(2, 3);
        assert_eq!(m, m2);
        let m2 = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 10, 6,
            7, 16, 27
            ].as_slice());
        m.scale_col_lt(1, 2);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_scale_row_ut(){
        let mut m = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9
            ].as_slice());
        let m2 = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 5, 6,
            7, 8, 27
            ].as_slice());
        m.scale_row_ut(2, 3);
        assert_eq!(m, m2);
        m.scale_row_ut(1, 2);
        let m2 = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 10, 12,
            7, 8, 27
            ].as_slice());
        assert_eq!(m, m2);
    }
    #[test]
    fn test_scale_col_ut(){
        let mut m = matrix_rw_i64(3, 3, [
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9
            ].as_slice());
        let m2 = matrix_rw_i64(3, 3, [
            1, 2, 9, 
            4, 5, 18,
            7, 8, 27
            ].as_slice());
        m.scale_col_ut(2, 3);
        assert_eq!(m, m2);
        let m2 = matrix_rw_i64(3, 3, [
            1, 4, 9, 
            4, 10, 18,
            7, 8, 27
            ].as_slice());
        m.scale_col_ut(1, 2);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_scale_rows(){
        let mut m = matrix_rw_i64(4, 3, [
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9,
            10, 11, 12
            ].as_slice());
        let factors = vector_i64([1, 2, 3, 4].as_slice());
        m.scale_rows(&factors);
        let m2 = matrix_rw_i64(4, 3, [
            1, 2, 3, 
            8, 10, 12,
            21, 24, 27,
            40, 44, 48
            ].as_slice());
        assert_eq!(m, m2);
    }
    #[test]
    fn test_scale_cols(){
        let mut m = matrix_rw_i64(2, 3, [
            1, 2, 3, 
            4, 5, 6,
            ].as_slice());
        let factors = vector_i64([1, 2, 3].as_slice());
        m.scale_cols(&factors);
        let m2 = matrix_rw_i64(2, 3, [
            1, 4, 9, 
            4, 10, 18,
            ].as_slice());
        assert_eq!(m, m2);
    }

}


/******************************************************
 *
 *   Bench marks
 *
 *******************************************************/


#[cfg(test)]
mod bench{

}



