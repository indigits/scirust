#![doc="Traits for matrices
"]





// Standard library imports
use std::cmp;
use std::num::{Int, SignedInt, Float};

// local imports
use number::Entry;
use number::{Number, Signed};
use matrix::matrix::Matrix;
use error::SRError;


// Reexports
pub use matrix::eo::eo_traits::{ERO, ECO};
pub use matrix::update::traits::{
    InPlaceUpdates, CopyUpdates};
pub use matrix::transpose::traits::{Transpose};
pub use matrix::extract::traits::{Extraction};


#[doc="Defines the features which all matrix types must implement.
This API focuses only on the shape of the matrix.
"] 
pub trait Shape<T:Entry> {
    
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

    /// Returns the minimum of rows and columns
    fn smaller_dim(&self) -> uint {
        cmp::min(self.num_rows(), self.num_cols())
    }

    /// Returns the larger of rows and columns
    fn larger_dim(&self)-> uint {
        cmp::max(self.num_rows(), self.num_cols())
    }

}


#[doc=" Defines the low level interface to the internal
memory buffer of a matrix implementation.  Use it with
caution. 
"]
pub trait MatrixBuffer<T:Entry> {

    /// Returns a constant pointer to matrix's buffer
    fn as_ptr(&self)-> *const T;

    /// Returns a mutable pointer to matrix's buffer
    fn as_mut_ptr(&mut self) -> *mut T;


    /// Maps a cell index to actual offset in buffer
    fn cell_to_offset(&self, r : uint,  c: uint)-> int;

    /// Returns the offset of the first cell in the buffer
    #[inline]
    fn start_offset(&self) -> int {
        0
    }
}

#[doc="A matrix structure whose storage is in terms
of columns with a fixed number of storage elements
per column given by its stride. 
"]
pub trait Strided {

    /// Returns the number of actual memory elements 
    /// per column
    fn stride (&self)->uint;


}


/// Defines a set of basic methods implemented by matrices of numbers
pub trait NumberMatrix<T:Number> : Shape<T> + MatrixBuffer<T>{
    
    /// Returns if the matrix is an identity matrix
    fn is_identity(&self) -> bool;

    /// Returns if the matrix is a diagonal matrix
    fn is_diagonal(&self) -> bool;

    /// Returns if the matrix is lower triangular 
    fn is_lt(&self) -> bool;

    /// Returns if the matrix is upper triangular 
    fn is_ut(&self) -> bool;


    /// Returns if the matrix is triangular
    fn is_triangular(&self) -> bool{
        self.is_lt() || self.is_ut()
    }

    /// Returns if the matrix is symmetric
    fn is_symmetric(&self) -> bool;

    /// Returns the trace of the matrix
    fn trace(&self) -> T ;
}

/// A matrix which is both a number matrix and is strided.
pub trait StridedNumberMatrix<T:Number> : 
NumberMatrix<T>+Strided{

}

#[doc=" These methods help in run time to query
about the properties of the type which is implementing
the trait. 

As an example, this API helps us decide if an object
of type Shape is a real matrix or a view 
extracted from a matrix. 

An implementing type needs to implement
only few of the methods below. For rest, it can
simply use the default implementation.

"]
pub trait Introspection {

    /// Indicates if the matrix is a view
    fn is_matrix_view_type(&self) -> bool {
        false
    }

    /// Indicates if the matrix is a standard matrix
    fn is_standard_matrix_type(&self) -> bool {
        false
    }

    /// Indicates if the matrix is a triangular matrix
    fn is_triangular_matrix_type(&self) -> bool {
        false
    }
}


/// Matrix conversion API
pub trait Conversion<T:Number> : Shape<T> {
    /// Converts the matrix to vector from standard library
    fn to_std_vec(&self) -> Vec<T>;

    /// Converts the matrix to a scalar 
    fn to_scalar(&self) -> T {
        if !self.is_scalar() {
            panic! (SRError::DimensionsMismatch.to_string());
        }
        self.get(0, 0)
    }
}



/// Matrix min-max API
pub trait MinMax<T:Number+PartialOrd> : Shape<T> {

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
pub trait MinMaxAbs<T:Signed> : Shape<T> {

    // Returns the absolute minimum scalar value
    fn min_abs_scalar(&self) -> (T, uint, uint);

    // Returns the maximum scalar value
    fn max_abs_scalar(&self) -> (T, uint, uint);


}

#[doc="Features for searching within the matrix
"]
pub trait Search<T:Signed+PartialOrd> : Shape<T>+MatrixBuffer<T> + Strided{

    /// Returns the largest entry (by magnitude) in the row between
    /// [start, end) columns
    fn max_abs_scalar_in_row(&self, 
        row : uint, 
        start_col: uint, 
        end_col : uint)-> (T, uint){
        debug_assert!(row < self.num_rows());
        debug_assert!(end_col <= self.num_cols());
        let p = self.as_ptr();
        let stride = self.stride();
        let mut offset = start_col * stride + row;
        let mut result = unsafe{*p.offset(offset as int)}.abs_val();
        let mut index  = 0;
        for i in range(1, end_col - start_col){
            offset += stride;
            let s = unsafe{*p.offset(offset as int)}.abs_val();
            if s > result {
                index = i;
                result = s;
            }
        }
        (result, index + start_col)
    }


    /// Returns the largest entry (by magnitude) in the column between
    /// [start, end) rows
    fn max_abs_scalar_in_col(&self, 
        col : uint, 
        start_row: uint, 
        end_row : uint)-> (T, uint){
        debug_assert!(end_row <= self.num_rows());
        debug_assert!(col < self.num_cols());
        let p = self.as_ptr();
        let stride = self.stride();
        let mut offset = col * stride + start_row;
        let mut result = unsafe{*p.offset(offset as int)}.abs_val();
        let mut index  = 0;
        for i in range(1, end_row - start_row){
            offset += 1;
            let s = unsafe{*p.offset(offset as int)}.abs_val();
            if s > result {
                index = i;
                result = s;
            }
        }
        (result, index + start_row)
    }
}


/// A matrix of integers
pub trait IntMatrix<T:Number+Int> {

}


/// A matrix of floats
pub trait FloatMatrix<T:Number+Float>{

}


/// A matrix of signed integers
pub trait SignedMatrix<T:Number+SignedInt>{

}


/// A strided buffer matrix of floats
pub trait StridedFloatMatrix <T:Number+Float> 
    : StridedNumberMatrix<T> {

}
