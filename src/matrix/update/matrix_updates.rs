/// std imports

/// local imports
use number::Number;
use error::SRResult;
use error::SRError;
use matrix::traits::{Shape, MatrixBuffer, Strided};
use matrix::update::traits::{InPlaceUpdates, CopyUpdates};
use matrix::eo::eo_traits::{ERO, ECO};
use matrix::matrix::Matrix;

/// Implementation of Matrix general update operations.
impl<T:Number> InPlaceUpdates<T> for Matrix<T> {
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

    /// Subtract a vector from each column
    fn sub_vec_from_cols(&mut self, vec: &Matrix<T>)->SRResult<()>{
        if ! vec.is_col() {
            return Err(SRError::IsNotAColVector);
        }
        let rows = self.num_rows();
        if vec.num_rows() != rows {
            return Err(SRError::RowsMismatch);
        }
        let cols = self.num_cols();
        let pd = self.as_mut_ptr();
        let ps = vec.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for _ in range(0, cols){
            for r in range(0i, rows as int){
                unsafe{
                    let v1 = *pd.offset(offset + r);
                    let v2 = *ps.offset(r);
                    *pd.offset(offset + r) = v1  - v2;
                }
            }
            offset += stride;
        }
        Ok(())
    }
    /// Subtract a vector from each row
    fn sub_vec_from_rows(&mut self, vec: &Matrix<T>)->SRResult<()>{
        if ! vec.is_row() {
            return Err(SRError::IsNotARowVector);
        }
        let cols = self.num_cols();
        if vec.num_cols() != cols {
            return Err(SRError::ColsMismatch);
        }
        let rows = self.num_rows();
        let pd = self.as_mut_ptr();
        let ps = vec.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for c in range(0, cols){
            let v2 = unsafe{*ps.offset(c as int)};
            for r in range(0i, rows as int){
                unsafe{
                    let v1 = *pd.offset(offset + r);
                    *pd.offset(offset + r) = v1  - v2;
                }
            }
            offset += stride;
        }
        Ok(())
    }
    /// Subtract a vector from each column
    fn add_vec_to_cols(&mut self, vec: &Matrix<T>)->SRResult<()>{
        if ! vec.is_col() {
            return Err(SRError::IsNotAColVector);
        }
        let rows = self.num_rows();
        if vec.num_rows() != rows {
            return Err(SRError::RowsMismatch);
        }
        let cols = self.num_cols();
        let pd = self.as_mut_ptr();
        let ps = vec.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for _ in range(0, cols){
            for r in range(0i, rows as int){
                unsafe{
                    let v1 = *pd.offset(offset + r);
                    let v2 = *ps.offset(r);
                    *pd.offset(offset + r) = v1  + v2;
                }
            }
            offset += stride;
        }
        Ok(())
    }
    /// Subtract a vector from each row
    fn add_vec_to_rows(&mut self, vec: &Matrix<T>)->SRResult<()>{
        if ! vec.is_row() {
            return Err(SRError::IsNotARowVector);
        }
        let cols = self.num_cols();
        if vec.num_cols() != cols {
            return Err(SRError::ColsMismatch);
        }
        let rows = self.num_rows();
        let pd = self.as_mut_ptr();
        let ps = vec.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for c in range(0, cols){
            let v2 = unsafe{*ps.offset(c as int)};
            for r in range(0i, rows as int){
                unsafe{
                    let v1 = *pd.offset(offset + r);
                    *pd.offset(offset + r) = v1  + v2;
                }
            }
            offset += stride;
        }
        Ok(())
    }
    /// Subtract a vector from each column
    fn mul_vec_to_cols(&mut self, vec: &Matrix<T>)->SRResult<()>{
        if ! vec.is_col() {
            return Err(SRError::IsNotAColVector);
        }
        let rows = self.num_rows();
        if vec.num_rows() != rows {
            return Err(SRError::RowsMismatch);
        }
        let cols = self.num_cols();
        let pd = self.as_mut_ptr();
        let ps = vec.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for _ in range(0, cols){
            for r in range(0i, rows as int){
                unsafe{
                    let v1 = *pd.offset(offset + r);
                    let v2 = *ps.offset(r);
                    *pd.offset(offset + r) = v1  * v2;
                }
            }
            offset += stride;
        }
        Ok(())
    }
    /// Subtract a vector from each row
    fn mul_vec_to_rows(&mut self, vec: &Matrix<T>)->SRResult<()>{
        if ! vec.is_row() {
            return Err(SRError::IsNotARowVector);
        }
        let cols = self.num_cols();
        if vec.num_cols() != cols {
            return Err(SRError::ColsMismatch);
        }
        let rows = self.num_rows();
        let pd = self.as_mut_ptr();
        let ps = vec.as_ptr();
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        for c in range(0, cols){
            let v2 = unsafe{*ps.offset(c as int)};
            for r in range(0i, rows as int){
                unsafe{
                    let v1 = *pd.offset(offset + r);
                    *pd.offset(offset + r) = v1  * v2;
                }
            }
            offset += stride;
        }
        Ok(())
    }

}


/// Implementation of Matrix general copy and update operations.
/// TODO Optimize implementations.
impl<T:Number> CopyUpdates<T> for Matrix<T> {

    /// Subtract a vector from each column
    fn copy_sub_vec_from_cols(&self, vec: &Matrix<T>)->SRResult<Matrix<T>>{
        let mut m = self.clone();
        let result = m.sub_vec_from_cols(vec);
        match result {
            Err(code) => Err(code),
            Ok(_) => Ok(m)
        }
    }
    /// Subtract a vector from each row
    fn copy_sub_vec_from_rows(&self, vec: &Matrix<T>)->SRResult<Matrix<T>>{
        let mut m = self.clone();
        let result = m.sub_vec_from_rows(vec);
        match result {
            Err(code) => Err(code),
            Ok(_) => Ok(m)
        }
    }
    /// Subtract a vector from each column
    fn copy_add_vec_to_cols(&self, vec: &Matrix<T>)->SRResult<Matrix<T>>{
        let mut m = self.clone();
        let result = m.add_vec_to_cols(vec);
        match result {
            Err(code) => Err(code),
            Ok(_) => Ok(m)
        }
    }
    /// Subtract a vector from each row
    fn copy_add_vec_to_rows(&self, vec: &Matrix<T>)->SRResult<Matrix<T>>{
        let mut m = self.clone();
        let result = m.add_vec_to_rows(vec);
        match result {
            Err(code) => Err(code),
            Ok(_) => Ok(m)
        }
    }
    /// Subtract a vector from each column
    fn copy_mul_vec_to_cols(&self, vec: &Matrix<T>)->SRResult<Matrix<T>>{
        let mut m = self.clone();
        let result = m.mul_vec_to_cols(vec);
        match result {
            Err(code) => Err(code),
            Ok(_) => Ok(m)
        }
    }
    /// Subtract a vector from each row
    fn copy_mul_vec_to_rows(&self, vec: &Matrix<T>)->SRResult<Matrix<T>>{
        let mut m = self.clone();
        let result = m.mul_vec_to_rows(vec);
        match result {
            Err(code) => Err(code),
            Ok(_) => Ok(m)
        }
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
    use matrix::update::traits::*;
    use matrix::traits::Transpose;

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

    #[test]
    fn test_sub_vec_from_cols(){
        let mut m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[1, 2]);
        assert!(m.sub_vec_from_cols(&v).is_ok());
        let m2 = matrix_rw_i64(2, 3, &[
            0, 1, 2, 
            2, 3, 4,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_copy_sub_vec_from_cols(){
        let m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[1, 2]);
        let m = m.copy_sub_vec_from_cols(&v).unwrap();
        let m2 = matrix_rw_i64(2, 3, &[
            0, 1, 2, 
            2, 3, 4,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_sub_vec_from_rows(){
        let mut m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2, 3]);
        assert!(m.sub_vec_from_rows(&v.transpose()).is_ok());
        let m2 = matrix_rw_i64(2, 3, &[
            2, 4, 0, 
            5, 7, 3,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_copy_sub_vec_from_rows(){
        let m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2, 3]);
        let m = m.copy_sub_vec_from_rows(&v.transpose()).unwrap();
        let m2 = matrix_rw_i64(2, 3, &[
            2, 4, 0, 
            5, 7, 3,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_add_vec_to_cols(){
        let mut m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2]);
        assert!(m.add_vec_to_cols(&v).is_ok());
        let m2 = matrix_rw_i64(2, 3, &[
            0, 1, 2, 
            2, 3, 4,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_copy_add_vec_to_cols(){
        let m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2]);
        let m = m.copy_add_vec_to_cols(&v).unwrap();
        let m2 = matrix_rw_i64(2, 3, &[
            0, 1, 2, 
            2, 3, 4,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_add_vec_to_rows(){
        let mut m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2, 3]);
        assert!(m.add_vec_to_rows(&v.transpose()).is_ok());
        let m2 = matrix_rw_i64(2, 3, &[
            0, 0, 6, 
            3, 3, 9,
            ]);
        assert_eq!(m, m2);
    }


    #[test]
    fn test_copy_add_vec_to_rows(){
        let m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2, 3]);
        let m = m.copy_add_vec_to_rows(&v.transpose()).unwrap();
        let m2 = matrix_rw_i64(2, 3, &[
            0, 0, 6, 
            3, 3, 9,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_mul_vec_to_cols(){
        let mut m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2]);
        assert!(m.mul_vec_to_cols(&v).is_ok());
        let m2 = matrix_rw_i64(2, 3, &[
            -1, -2, -3, 
            -8, -10, -12,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_copy_mul_vec_to_cols(){
        let m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2]);
        let m = m.copy_mul_vec_to_cols(&v).unwrap();
        let m2 = matrix_rw_i64(2, 3, &[
            -1, -2, -3, 
            -8, -10, -12,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_mul_vec_to_rows(){
        let mut m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2, 3]);
        assert!(m.mul_vec_to_rows(&v.transpose()).is_ok());
        let m2 = matrix_rw_i64(2, 3, &[
            -1, -4, 9,
            -4, -10, 18
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_copy_mul_vec_to_rows(){
        let m = matrix_rw_i64(2, 3, &[
            1, 2, 3, 
            4, 5, 6,
            ]);
        let v = vector_i64(&[-1, -2, 3]);
        let m = m.copy_mul_vec_to_rows(&v.transpose()).unwrap();
        let m2 = matrix_rw_i64(2, 3, &[
            -1, -4, 9,
            -4, -10, 18
            ]);
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



