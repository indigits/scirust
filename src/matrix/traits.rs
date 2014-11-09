// Standard library imports
use std::ptr;


// local imports
use matrix::element::Number;
use matrix::matrix::Matrix;
use matrix::error::*;

use discrete::{mod_n};


#[doc="Defines the traits which all matrix types must implement.

At the moment there are two matrix types. One is the `Matrix<T>` itself
and another is `MatrixView<T>`. It so happens that they share
a lot of methods in common.  This trait identifies a fundamental
set of methods which every matrix type must implement.
"] 
pub trait MatrixType<T:Number> {
    
    /// Returns the number of rows
    fn num_rows(&self) -> uint;

    /// Returns the number of columns
    fn num_cols(&self) -> uint;


    /// Returns the size in an (r, c) tuple
    fn size (&self)-> (uint, uint);
 
 
     /// Returns the number of cells in matrix
    fn num_cells(&self)->uint;

    /// Indicates if the matrix is a row vector
    fn is_row(&self) -> bool {
        self.num_rows() == 1
    }

    /// Indicates if the matrix is a column vector
    fn is_col(&self) -> bool {
        self.num_cols() == 1
    }

    /// Indicates if the matrix is a scalar actually
    fn is_scalar(&self) -> bool {
        self.num_cells() == 1
    }

    /// Indicates if the matrix is a vector
    fn is_vector(&self) -> bool {
        (self.num_rows() == 1) ^ (self.num_cols() == 1)
    } 

    /// Indicates if the matrix is empty
    fn is_empty(&self) -> bool {
        self.num_rows() * self.num_cols() == 0
    }

    /// Indicates if the matrix is square
    fn is_square(&self) -> bool {
        self.num_rows() == self.num_cols()
    }

    /// Gets an element by its row and column number
    fn get(&self, r : uint, c : uint) -> T;

    /// Sets an element in the view
    fn set(&mut self, r : uint, c : uint, value : T);


    /// Converts an index to cell address (row, column)
    #[inline]
    fn index_to_cell(&self, index : uint) -> (uint, uint){
        debug_assert!(index < self.num_cells());
        let rows = self.num_rows();
        let c = index / rows;
        let r = index - c*rows;
        (r, c)
    }

    /// Converts a cell address to an index (r, c) to index
    #[inline]
    fn cell_to_index(&self, r : uint,  c: uint) -> uint{
        let rows = self.num_rows();
        debug_assert!(r < rows);
        debug_assert!(c < self.num_cols());
        c * rows + r
    }

}


#[doc=" Defines the low level interface to the internal
memory buffer of a matrix implementation.  Use it with
caution. 
"]
pub trait MatrixBuffer<T:Number> {

    /// Returns the number of actual memory elements 
    /// per column
    fn stride (&self)->uint;


    /// Returns a constant pointer to matrix's buffer
    fn as_ptr(&self)-> *const T;

    /// Returns a mutable pointer to matrix's buffer
    fn as_mut_ptr(&mut self) -> *mut T;


    /// Maps a cell index to actual offset in buffer
    fn cell_to_offset(&self, r : uint,  c: uint)-> int;

}

#[doc=" These methods help in run time to query
about the properties of the type which is implementing
the trait. 

As an example, this API helps us decide if an object
of type MatrixType is a real matrix or a view 
extracted from a matrix. 

An implementing type needs to implement
only few of the methods below. For rest, it can
simply use the default implementation.

"]
pub trait Introspection {

    /// Indicates if the matrix is a view
    fn is_view(&self) -> bool {
        false
    }

    /// Indicates if the matrix is a standard matrix
    fn is_matrix(&self) -> bool {
        false
    }

}


/// Matrix conversion API
pub trait Conversion<T:Number> : MatrixType<T> {
    /// Converts the matrix to vector from standard library
    fn to_std_vec(&self) -> Vec<T>;

    /// Converts the matrix to a scalar 
    fn to_scalar(&self) -> T {
        if !self.is_scalar() {
            fail! (DimensionsMismatch.to_string());
        }
        self.get(0, 0)
    }
}


/// Matrix extraction API
pub trait Extraction<T:Number> : MatrixType<T>+MatrixBuffer<T> {

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

}

/// Matrix min-max API
pub trait MinMax<T:Number+PartialOrd> : MatrixType<T> {

    /// Returns a column vector consisting of maximum over each row
    fn max_row_wise(&self) -> Matrix<T>;

    /// Returns a column vector consisting of minimum over each row
    fn min_row_wise(&self) -> Matrix<T>;

    /// Returns a row vector consisting of maximum over each column
    fn max_col_wise(&self) -> Matrix<T>;

    /// Returns a row vector consisting of minimum over each column
    fn min_col_wise(&self) -> Matrix<T>;
}


/// Matrix min-max with absolute values API
pub trait MinMaxAbs<T:Number+PartialOrd+Signed> : MatrixType<T> {

    // Returns the absolute minimum scalar value
    fn min_abs_scalar(&self) -> (T, uint, uint);

    // Returns the maximum scalar value
    fn max_abs_scalar(&self) -> (T, uint, uint);


}


/******************************************************
 *
 *   Elementary row / column operations.
 *
 *******************************************************/

/// Elementary row operations on a matrix
pub trait ERO<T:Number> : MatrixType<T>+MatrixBuffer<T> {

    /// Row switching.
    fn ero_switch(&mut self, 
        i :  uint,
        j : uint
        )-> &mut Self {
        debug_assert! (i  < self.num_rows());
        debug_assert! (j  < self.num_rows());
        let ptr = self.as_mut_ptr();
        for c in range(0, self.num_cols()){
            let offset_a = self.cell_to_offset(i, c);
            let offset_b = self.cell_to_offset(j, c);
            unsafe {
                ptr::swap(ptr.offset(offset_a), ptr.offset(offset_b));
            }
        }
        self
    }

    /// Row scaling by a factor.
    fn ero_scale(&mut self, 
        r :  uint, 
        scale : T
        )-> &mut Self {
        debug_assert! (r  < self.num_rows());
        let ptr = self.as_mut_ptr();
        for c in range(0, self.num_cols()){
            let offset = self.cell_to_offset(r, c);
            unsafe {
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale * v;
            }
        }
        self
    }

    /// Row scaling by a factor and adding to another row.
    /// r_i = r_i + k * r_j
    fn ero_scale_add(&mut self, 
        i :  uint, 
        j :  int, 
        scale : T
        )-> &mut Self {
        debug_assert! (i  < self.num_rows());
        debug_assert! (j  < self.num_rows() as int);
        let j = if j < 0{
            mod_n(j, self.num_rows() as int)
        }
        else {
            j as uint
        };
        let ptr = self.as_mut_ptr();
        for c in range(0, self.num_cols()){
            let offset_a = self.cell_to_offset(i, c);
            let offset_b = self.cell_to_offset(j, c);
            unsafe {
                let va = *ptr.offset(offset_a);
                let vb = *ptr.offset(offset_b);
                *ptr.offset(offset_a) = va + scale * vb;
            }
        }
        self
    }

}