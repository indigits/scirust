
// local imports
use matrix::element::MatrixElt;
use matrix::matrix::Matrix;
use matrix::error::*;

#[doc="Defines the traits which all matrix types must implement.

At the moment there are two matrix types. One is the `Matrix<T>` itself
and another is `MatrixView<T>`. It so happens that they share
a lot of methods in common.  This trait identifies a fundamental
set of methods which every matrix type must implement.
"] 
pub trait MatrixType<T:MatrixElt> {
    
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
pub trait MatrixBuffer<T:MatrixElt> {

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
pub trait Conversion<T:MatrixElt> : MatrixType<T> {
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
pub trait Extraction<T:MatrixElt> : MatrixType<T> {

    /// Returns the r'th row vector
    fn row(&self, r : int) -> Matrix<T>;

    /// Returns the c'th column vector
    fn col(&self, c : int) -> Matrix<T>;

    /// Extract a submatrix from the matrix
    /// rows can easily repeat if the number of requested rows is higher than actual rows
    /// cols can easily repeat if the number of requested cols is higher than actual cols
    fn sub_matrix(&self, start_row : int, 
        start_col : int , 
        num_rows: uint, 
        num_cols : uint) -> Matrix<T>;

}

/// Matrix min-max API
pub trait MinMax<T:MatrixElt+PartialOrd> : MatrixType<T> {

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
pub trait MinMaxAbs<T:MatrixElt+PartialOrd+Signed> : MatrixType<T> {

    // Returns the absolute minimum scalar value
    fn min_abs_scalar(&self) -> (T, uint, uint);

    // Returns the maximum scalar value
    fn max_abs_scalar(&self) -> (T, uint, uint);


}
