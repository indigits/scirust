#![doc="Implementation of extraction API for matrix type
"]

// local imports
use mod_n;
use sralgebra::{CommutativeMonoidAddPartial};
use num::traits::Zero;
use matrix::Matrix;
use traits::{Shape, MatrixBuffer, Strided};
use extract::traits::Extraction;


/// Implements matrix extraction API
impl <T:CommutativeMonoidAddPartial> Extraction<T> for Matrix<T> {

    /// Returns the r'th row vector
    fn row(&self, r : isize) -> Matrix<T> {
        // Lets ensure that the row value is mapped to
        // a value in the range [0, rows - 1]
        let r = mod_n(r, self.num_rows() as isize);        
        let mut result : Matrix<T> = Matrix::new_uninitialized(1, self.num_cols());
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        for c in 0..self.num_cols(){
            let src_offset = self.cell_to_offset(r, c);
            let dst_offset = result.cell_to_offset(0, c);
            unsafe{
                *pd.offset(dst_offset) = *ps.offset(src_offset);
            }
        }
        result
    }

    /// Returns the c'th column vector
    fn col(&self, c : isize) -> Matrix<T>{
        // Lets ensure that the col value is mapped to
        // a value in the range [0, cols - 1]
        let c = mod_n(c, self.num_cols() as isize);        
        let mut result : Matrix<T> = Matrix::new_uninitialized(self.num_rows(), 1);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        for r in 0..self.num_rows(){
            let src_offset = self.cell_to_offset(r, c);
            let dst_offset = result.cell_to_offset(r, 0);
            unsafe{
                *pd.offset(dst_offset) = *ps.offset(src_offset);
            }
        }
        result
    }

    /// Extract a submatrix from the matrix
    /// rows can easily repeat if the number of requested rows is higher than actual rows
    /// cols can easily repeat if the number of requested cols is higher than actual cols
    fn sub_matrix(&self, start_row : isize, 
        start_col : isize , 
        num_rows: usize, 
        num_cols : usize) -> Matrix<T>{
        let r = mod_n(start_row, self.num_rows() as isize);        
        let c = mod_n(start_col, self.num_cols() as isize);
        let mut result : Matrix<T> = Matrix::new_uninitialized(num_rows, num_cols);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        let mut dc = 0;
        for c in (c..(c + num_cols)).map(|x | x % self.num_cols()) {
            let mut dr = 0;
            for r in (r..(r + num_rows)).map(|x|  x % self.num_rows()) {
                let src_offset = self.cell_to_offset(r, c);
                let dst_offset = result.cell_to_offset(dr, dc);
                unsafe{
                    *pd.offset(dst_offset) = *ps.offset(src_offset);
                }
                dr += 1;
            }
            dc += 1;
        }
        result
    }


    /// Returns the lower triangular part of the matrix
    fn lt_matrix(&self)->Matrix<T>{
        let n = self.smaller_dim();
        let mut result = Matrix::new_uninitialized(n, n);
        let z  : T = Zero::zero();
        let mut ps = self.as_ptr();
        let mut pd = result.as_mut_ptr();
        let src_stride = self.stride() as isize;
        let dst_stride = result.stride() as isize;
        for c in 0..n{
            for r in 0..c{
                unsafe{
                    *pd.offset(r as isize) = z;
                }
            }
            for r in c..n{
                let r  = r as isize;
                unsafe {
                    *pd.offset(r) = *ps.offset(r);
                }
            }
            unsafe
            {
                ps = ps.offset(src_stride);
                pd = pd.offset(dst_stride);
            }
        } 
        result
    }

    /// Returns the upper triangular part of the matrix
    fn ut_matrix(&self)->Matrix<T>{
        let n = self.smaller_dim();
        let mut result = Matrix::new_uninitialized(n, n);
        let z  : T = Zero::zero();
        let mut ps = self.as_ptr();
        let mut pd = result.as_mut_ptr();
        let src_stride = self.stride() as isize;
        let dst_stride = result.stride() as isize;
        for c in 0..n{
            for r in 0..(c + 1){
                let r  = r as isize;
                unsafe {
                    *pd.offset(r) = *ps.offset(r);
                }
            }
            for r in (c + 1)..n{
                unsafe{
                    *pd.offset(r as isize) = z;
                }
            }
            unsafe
            {
                ps = ps.offset(src_stride);
                pd = pd.offset(dst_stride);
            }
        } 
        result
    }

}


#[cfg(test)]
mod test{
    use matrix::*;
    use constructors::*;
    use traits::*;


    #[test]
    fn test_extract_row(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(4, 4, (0..16));
        let m2  = m1.row(0);
        assert_eq!(m2.to_std_vec(), vec![0, 4, 8, 12]);
        assert_eq!(m2.num_rows() , 1);
        assert!(m2.is_row());
        assert_eq!(m2.num_cols() , m1.num_cols());
        assert_eq!(m2.num_cells() , m1.num_cols());
        assert_eq!(m1.row(1).to_std_vec(), vec![1, 5, 9, 13]);
        assert_eq!(m1.row(2).to_std_vec(), vec![2, 6, 10, 14]);
        assert_eq!(m1.row(3).to_std_vec(), vec![3, 7, 11, 15]);
        assert_eq!(m1.row(-1).to_std_vec(), vec![3, 7, 11, 15]);
        assert_eq!(m1.row(-2).to_std_vec(), vec![2, 6, 10, 14]);
        assert_eq!(m1.row(-6).to_std_vec(), vec![2, 6, 10, 14]);
    }


    #[test]
    fn test_extract_col(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(4, 4, (0..16));
        let m2  = m1.col(0);
        assert_eq!(m2.to_std_vec(), vec![0, 1, 2, 3]);
        assert!(!m2.is_row());
        assert!(m2.is_col());
        assert!(!m2.is_scalar());
        assert_eq!(m2.num_rows() , m1.num_rows());
        assert_eq!(m2.num_cells() , m1.num_rows());
        assert_eq!(m1.col(1).to_std_vec(), vec![4, 5, 6, 7]);
        assert_eq!(m1.col(2).to_std_vec(), vec![8, 9, 10, 11]);
        assert_eq!(m1.col(3).to_std_vec(), vec![12, 13, 14, 15]);
        assert_eq!(m1.col(-1).to_std_vec(), vec![12, 13, 14, 15]);
        assert_eq!(m1.col(-2).to_std_vec(), vec![8, 9, 10, 11]);
        assert_eq!(m1.col(-6).to_std_vec(), vec![8, 9, 10, 11]);
    }


    #[test]
    fn test_sub_matrix(){
        let m  : MatrixI64 = Matrix::from_iter_cw(4, 4, (0..16));
        let m1 = m.sub_matrix(0, 0, 2, 2);
        assert_eq!(m1.num_cells(), 4);
        assert_eq!(m1.num_rows(), 2);
        assert_eq!(m1.num_cols(), 2);
        assert_eq!(m1.to_std_vec(), vec![0, 1, 4, 5]);
        assert_eq!(m.sub_matrix(1, 0, 2, 2).to_std_vec(), vec![1, 2, 5, 6]);
        assert_eq!(m.sub_matrix(4, 4, 2, 2).to_std_vec(), vec![0, 1, 4, 5]);
        assert_eq!(m.sub_matrix(-1, -1, 2, 2).to_std_vec(), vec![15, 12, 3, 0]);
        assert_eq!(m.sub_matrix(-4, -4, 2, 2).to_std_vec(), vec![0, 1, 4, 5]);
    }

    #[test]
    fn test_ut(){
        let m = matrix_rw_i32(4, 2, &[
            1, 2,
            3, 4,
            5, 6, 
            7, 8
            ]);
        assert_eq!(m.ut_matrix(), 
            matrix_rw_i32(2,2, &[
            1, 2,
            0, 4]));
        assert_eq!(m.lt_matrix(), 
            matrix_rw_i32(2,2, &[
            1, 0,
            3, 4]));
    }

}


/******************************************************
 *
 *   Bench marks
 *
 *******************************************************/


#[cfg(test)]
mod bench{
    //extern crate test;
    //use self::test::Bencher;
    //use super::*;
}


