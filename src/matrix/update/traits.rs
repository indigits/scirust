#![doc="Various kind of updates to the matrix
"]

/// std imports
use std::cmp;


/// local imports
use discrete::mod_n;
use number::Number;
use matrix::matrix::Matrix;
use matrix::traits::{Shape, MatrixBuffer};


#[doc="Various kind of updates to the matrix

Such updates include: 

- Setting all the entries on a diagonal
- Setting all the entries on a row
- Setting all the entries on a column 
- etc.

These are not standard operations, but since they
are needed once in a while, hence the helper functions
are provided to achieve them quickly.

"]
pub trait Updates<T:Number> : Shape<T>+MatrixBuffer<T> {

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
        r :  uint, 
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
        c :  uint, 
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
        start_row : int, 
        start_col : int , 
        num_rows: uint, 
        num_cols : uint,
        value : T
        )-> &mut Self {
        let r = mod_n(start_row, self.num_rows() as int);        
        let c = mod_n(start_col, self.num_cols() as int);
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
    /// Scales a specific row in lower triangular part 
    fn scale_row_lt(&mut self, r :  uint, scale_factor : T)-> &mut Self;
    /// Scales a specific column in lower triangular part
    fn scale_col_lt(&mut self, c :  uint, scale_factor : T)-> &mut Self;
    /// Scales a specific row in upper triangular part
    fn scale_row_ut(&mut self, r :  uint, scale_factor : T)-> &mut Self;
    /// Scales a specific column in upper triangular part
    fn scale_col_ut(&mut self, c :  uint, scale_factor : T)-> &mut Self;
    /// Scale all rows as per the scale factors
    fn scale_rows(&mut self, scale_factors : &Matrix<T>)-> &mut Self;
    /// Scale all columns as per the scale factors
    fn scale_cols(&mut self, scale_factors : &Matrix<T>)-> &mut Self;
}


