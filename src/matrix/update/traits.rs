#![doc="Various kind of updates to the matrix
"]

/// std imports
use std::cmp;


/// local imports
use error::SRResult;
use discrete::mod_n;
use algebra::Number;
use matrix::matrix::Matrix;
use matrix::traits::{Shape, MatrixBuffer};


#[doc="Matrix In Place Updates API

This trait defines various updates to the contents
of a matrix in-place. Such updates include: 

- Setting all the entries on a diagonal
- Setting all the entries on a row
- Setting all the entries on a column 
- etc.

These are not standard linear algebra operations, but since they
are needed once in a while, hence the helper functions
are provided to achieve them quickly. Further,
optimizing the low level code can help a lot in improving performance
of higher level functions.
"]
pub trait InPlaceUpdates<T:Number> : Shape<T>+MatrixBuffer<T> {

    // Sets all the entries on the main diagonal to a particular value
    fn set_diagonal(&mut self, value : T)-> &mut Self{
        let n = cmp::min(self.num_rows(), self.num_cols());
        let p = self.as_mut_ptr();
        for i in range(0, n){
            unsafe{
                *p.offset(self.cell_to_offset(i, i)) = value;
            }
        }
        self
    }

   /// Sets all the entries on a row to a particular value
    fn set_row(&mut self, 
        r :  usize, 
        value : T
        )-> &mut Self {
        debug_assert! (r  < self.num_rows());
        let ptr = self.as_mut_ptr();
        for c in range(0, self.num_cols()){
            unsafe{
                *ptr.offset(self.cell_to_offset(r, c)) = value;
            }
        }
        self
    }


   /// Sets all the entries on a column to a particular value
    fn set_col(&mut self, 
        c :  usize, 
        value : T
        )-> &mut Self {
        debug_assert! (c  < self.num_cols());
        let ptr = self.as_mut_ptr();
        for r in range(0, self.num_rows()){
            unsafe{
                *ptr.offset(self.cell_to_offset(r, c)) = value;
            }
        }
        self
    }

   /// Sets all the entries on a block to a particular value
    fn set_block(&mut self, 
        start_row : isize, 
        start_col : isize , 
        num_rows: usize, 
        num_cols : usize,
        value : T
        )-> &mut Self {
        let r = mod_n(start_row, self.num_rows() as isize);        
        let c = mod_n(start_col, self.num_cols() as isize);
        assert!(num_rows < self.num_rows());
        assert!(num_cols < self.num_cols());
        let p = self.as_mut_ptr();
        for c in range(c, c + num_cols).map(|x | x % self.num_cols()) {
            for r in range(r , r + num_rows).map(|x|  x % self.num_rows()) {
                let offset = self.cell_to_offset(r, c);
                unsafe{
                    *p.offset(offset) = value;
                }
            }
        }
        self
    }
     /// Add the matrix by a scalar
    fn add_scalar(&mut self, rhs: T) -> &mut Self;
     /// Multiply the matrix by a scalar
    fn mul_scalar(&mut self, rhs: T) -> &mut Self;
     /// Divide the matrix by a scalar
    fn div_scalar(&mut self, rhs: T) -> SRResult<&mut Self>;
    /// Scales a specific row in lower triangular part 
    fn scale_row_lt(&mut self, r :  usize, scale_factor : T)-> &mut Self;
    /// Scales a specific column in lower triangular part
    fn scale_col_lt(&mut self, c :  usize, scale_factor : T)-> &mut Self;
    /// Scales a specific row in upper triangular part
    fn scale_row_ut(&mut self, r :  usize, scale_factor : T)-> &mut Self;
    /// Scales a specific column in upper triangular part
    fn scale_col_ut(&mut self, c :  usize, scale_factor : T)-> &mut Self;
    /// Scale all rows as per the scale factors
    fn scale_rows(&mut self, scale_factors : &Matrix<T>)-> &mut Self;
    /// Scale all columns as per the scale factors
    fn scale_cols(&mut self, scale_factors : &Matrix<T>)-> &mut Self;
    /// Subtract a vector from each column
    fn sub_vec_from_cols(&mut self, vec: &Matrix<T>)->SRResult<()>;
    /// Subtract a vector from each row
    fn sub_vec_from_rows(&mut self, vec: &Matrix<T>)->SRResult<()>;
    /// Subtract a vector from each column
    fn add_vec_to_cols(&mut self, vec: &Matrix<T>)->SRResult<()>;
    /// Subtract a vector from each row
    fn add_vec_to_rows(&mut self, vec: &Matrix<T>)->SRResult<()>;
    /// Subtract a vector from each column
    fn mul_vec_to_cols(&mut self, vec: &Matrix<T>)->SRResult<()>;
    /// Subtract a vector from each row
    fn mul_vec_to_rows(&mut self, vec: &Matrix<T>)->SRResult<()>;

    /// Copies data from upper triangular part to lower triangular part
    fn ut_to_lt(&mut self)->&mut Self;
}


#[doc=" Matrix Copy and Update API 

These methods create a new copy of the matrix
and apply changes to it. Care is taken to
avoid extra operations.  The copying process
and update process are merged together to 
save processing cycles.

# Remarks

These methods require immutable reference to the
matrix being updated.
"]
pub trait CopyUpdates<T:Number> : Shape<T>+MatrixBuffer<T> {
     /// Add the matrix by a scalar
    fn copy_add_scalar(&self, rhs: T) -> Self;
     /// Multiply the matrix by a scalar
    fn copy_mul_scalar(&self, rhs: T) -> Self;
     /// Divide the matrix by a scalar
    fn copy_div_scalar(&self, rhs: T) -> Self;
    /// Subtract a vector from each column
    fn copy_sub_vec_from_cols(&self, vec: &Matrix<T>)->SRResult<Self>;
    /// Subtract a vector from each row
    fn copy_sub_vec_from_rows(&self, vec: &Matrix<T>)->SRResult<Self>;
    /// Subtract a vector from each column
    fn copy_add_vec_to_cols(&self, vec: &Matrix<T>)->SRResult<Self>;
    /// Subtract a vector from each row
    fn copy_add_vec_to_rows(&self, vec: &Matrix<T>)->SRResult<Self>;
    /// Subtract a vector from each column
    fn copy_mul_vec_to_cols(&self, vec: &Matrix<T>)->SRResult<Self>;
    /// Subtract a vector from each row
    fn copy_mul_vec_to_rows(&self, vec: &Matrix<T>)->SRResult<Self>;
}
