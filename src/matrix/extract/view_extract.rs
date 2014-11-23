#![doc="Implementation of extraction API for matrix view type
"]


// local imports
use number::Number;
use matrix::matrix::Matrix;
use matrix::view::MatrixView;
use matrix::extract::traits::Extraction;
use matrix::traits::{Shape, MatrixBuffer};

use discrete::{mod_n};


/// Implement extraction API for matrix view 
impl <'a, T:Number> Extraction<T> for MatrixView<'a, T> {

    /// Returns the r'th row vector
    fn row(&self, r : int) -> Matrix<T> {
        // Lets ensure that the row value is mapped to
        // a value in the range [0, rows - 1]
        let r = mod_n(r, self.num_rows() as int);        
        let mut result : Matrix<T> = Matrix::new(1, self.num_cols());
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        for c in range(0, self.num_cols()){
            let src_offset = self.cell_to_offset(r, c);
            let dst_offset = result.cell_to_offset(0, c);
            unsafe{
                *pd.offset(dst_offset) = *ps.offset(src_offset);
            }
        }
        result
    }

    /// Returns the c'th column vector
    fn col(&self, c : int) -> Matrix<T>{
        // Lets ensure that the col value is mapped to
        // a value in the range [0, cols - 1]
        let c = mod_n(c, self.num_cols() as int);        
        let mut result : Matrix<T> = Matrix::new(self.num_rows(), 1);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        for r in range(0, self.num_rows()){
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
    fn sub_matrix(&self, start_row : int, 
        start_col : int , 
        num_rows: uint, 
        num_cols : uint) -> Matrix<T>{
        let r = mod_n(start_row, self.num_rows() as int);        
        let c = mod_n(start_col, self.num_cols() as int);
        let mut result : Matrix<T> = Matrix::new(num_rows, num_cols);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        let mut dc = 0;
        for c in range(c, c + num_cols).map(|x | x % self.num_cols()) {
            let mut dr = 0;
            for r in range(r , r + num_rows).map(|x|  x % self.num_rows()) {
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

    /// Returns the upper triangular part of the matrix
    fn ut_matrix(&self)->Matrix<T>{
        unimplemented!();
    }

    /// Returns the lower triangular part of the matrix
    fn lt_matrix(&self)->Matrix<T>{
        unimplemented!();
    }

}


#[cfg(test)]
mod test{
    use api::*;
    #[test]
    fn test_extract_row(){
        let m :  MatrixI64 = Matrix::from_iter_cw(20, 20, range(-100, 400));
        let v1   = m.view(2, 2, 6, 6);
        println!("v1 : {}", v1);
        let r1 = v1.row(0);
        let v2 = m.view(2,2, 1, 6);
        assert_eq!(v2.to_matrix(), r1);
    }

    #[test]
    fn test_extract_col(){
        let m :  MatrixI64 = from_range_rw_i64(20, 20, -100, 400);
        let v1   = m.view(2, 2, 6, 6);
        println!("v1 : {}", v1);
        let c1 = v1.col(0);
        let v2 = m.view(2,2, 6, 1);
        assert_eq!(v2.to_matrix(), c1);
    }

    #[test]
    fn test_extract_sub_matrix(){
        let m :  MatrixI64 = from_range_rw_i64(20, 20, -100, 400);
        let v1   = m.view(2, 2, 6, 6);
        let m2  = v1.sub_matrix(2, 2, 4, 4);
        let m3 = m.sub_matrix(4,4, 4, 4);
        println!("m2 : {}", m2);
        assert_eq!(m2, m3);
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



