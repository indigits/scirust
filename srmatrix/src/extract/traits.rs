
// local imports
use sralgebra::{CommutativeMonoidAddPartial};
use matrix::Matrix;
use traits::{Shape};




/// Matrix extraction API
pub trait Extraction<T:CommutativeMonoidAddPartial> : Shape<T>{

    /// Returns the r'th row vector
    fn row(&self, r : isize) -> Matrix<T>;

    /// Returns the c'th column vector
    fn col(&self, c : isize) -> Matrix<T>;

    /// Extract a submatrix from the matrix
    /// rows can easily repeat if the number of requested rows is higher than actual rows
    /// cols can easily repeat if the number of requested cols is higher than actual cols
    fn sub_matrix(&self, start_row : isize, 
        start_col : isize , 
        num_rows: usize, 
        num_cols : usize) -> Matrix<T>;

    /// Returns the upper triangular part of the matrix
    fn ut_matrix(&self)->Matrix<T>;

    /// Returns the lower triangular part of the matrix
    fn lt_matrix(&self)->Matrix<T>;
}
