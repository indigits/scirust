#![doc="Traits describing statistical moments of data
"]

// std imports
use num::{Float, FromPrimitive};

// local imports
use matrix::matrix::Matrix;
use algebra::structure::{CommutativeMonoidAddPartial, CommutativeMonoidMulPartial, FieldPartial};


pub trait Sums<T: CommutativeMonoidAddPartial+CommutativeMonoidMulPartial> {

    /// Computes sum over columns and returns a row vector
    fn sum_cw(&self) -> Matrix<T>;

    /// Computes sum over rows and returns a column vector
    fn sum_rw(&self) -> Matrix<T>;

    /// Computes sum over columns and returns a row vector
    fn sum_sqr_cw(&self) -> Matrix<T>;

    /// Computes sum over rows and returns a column vector
    fn sum_sqr_rw(&self) -> Matrix<T>;
}


pub trait Moments <T: FieldPartial + Float + FromPrimitive> {

    /// Computes mean over columns and returns a row vector
    fn mean_cw(&self) -> Matrix<T>;

    /// Computes mean over rows and returns a column vector
    fn mean_rw(&self) -> Matrix<T>;

    /// Computes mean square over columns and returns a row vector
    fn mean_sqr_cw(&self) -> Matrix<T>;

    /// Computes mean square over rows and returns a column vector
    fn mean_sqr_rw(&self) -> Matrix<T>;

    /// Computes sample variance over columns and returns a row vector
    /// sum((m - repmat(mean(m), r, 1)).^2 ) / (r - 1)
    fn var_cw(&self) -> Matrix<T>;

    /// Computes sample variance over rows and returns a column vector
    fn var_rw(&self) -> Matrix<T>;


#[doc="Computes the sample covariance from a set of example vectors.
When self is a proper matrix, then we consider each row as an observation
and each column as a random variable. 

If x is a vector, then we return the variance of the vector.
"]
    fn cov(&self) -> Matrix<T>;
}
