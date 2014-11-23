
// local imports
use number::Number;
use matrix::matrix::Matrix;
use matrix::traits::{Shape};




/// Matrix extraction API
pub trait Extraction<T:Number> : Shape<T>{

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

    /// Returns the upper triangular part of the matrix
    fn ut_matrix(&self)->Matrix<T>;

    /// Returns the lower triangular part of the matrix
    fn lt_matrix(&self)->Matrix<T>;
}
