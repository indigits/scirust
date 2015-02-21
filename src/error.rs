#![doc="Defines errors related to SciRust library
"]


// std imports


// local imports



/// Enum for errors related to SciRust library
#[derive(Debug, Copy, Clone)]
pub enum SRError{

    /******************************************************
     *
     *   General matrix operation errors
     *
     *******************************************************/

    /// The matrix is empty
    EmptyMatrix,
    /// The dimensions of two matrices mismatch
    DimensionsMismatch,
    /// Number of rows don't match
    RowsMismatch,
    /// Number of columns don't match
    ColsMismatch,
    /// The matrix is not a square matrix
    IsNotSquareMatrix,
    /// The matrix is not a vector
    IsNotAVector,
    // The matrix is not a column vector
    IsNotAColVector,
    // The matrix is not a row vector
    IsNotARowVector,
    /// Indicates that a matrix is not full rank
    IsNotFullRankMatrix,
    /// Indicates that the matrix is full rank
    IsFullRankMatrix,
    /// The matrix is singular
    IsSingular,
    /// The matrix is invertible (non-singular)
    IsNonSingular,
    /// The matrix is positive definite
    IsPositiveDefinite,
    /// The matrix is positive semi definite
    IsPositiveSemiDefinite,
    /// The matrix is negative definite
    IsNegativeDefinite,
    /// The matrix is negative semi definite
    IsNegativeSemiDefinite,
    /// The matrix is non-definite
    IsNonDefinite,


    /******************************************************
     *
     *   Errors related to solving a system of linear equations
     *
     *******************************************************/

    /// The dimensions of left and right hand side don't match
    LRDimensionMismatch,
    /// There is no solution to the system of equations
    NoSolution,
    /// There are infinite solutions to the system of equations.
    InfiniteSolutions,

    /******************************************************
     *
     *   Arithmetic related stuff
     *
     *******************************************************/
    //  Divide by zero
    DivideByZero,


    /******************************************************
     *
     *   Discrete numbers related stuff
     *
     *******************************************************/

    /// The number is not a power of two
    IsNotPowerOfTwo,

}

///A convenient typedef of the return value of SciRust API
/// whenever applicable.
pub type SRResult<T> = Result<T, SRError>;



impl SRError{

    /// Converts enum values to string representation
    pub fn to_string(&self) -> String {
        match *self{
            //  Matrices
            SRError::EmptyMatrix => format!("Matrix is empty"),
            SRError::DimensionsMismatch => format!("Dimensions don't match"),
            SRError::RowsMismatch => format!("Number of rows don't match"),
            SRError::ColsMismatch => format!("Number of columns don't match"),
            SRError::IsNotSquareMatrix => format!("Matrix is not square"),
            SRError::IsNotAVector => format!("Matrix is not a vector"),
            SRError::IsNotAColVector => format!("Matrix is not a column vector"),
            SRError::IsNotARowVector => format!("Matrix is not a row vector"),
            SRError::IsNotFullRankMatrix => format!("Matrix is not full rank"),
            SRError::IsFullRankMatrix => format!("Matrix is full rank"),
            SRError::IsSingular => format!("Matrix is singular"),
            SRError::IsNonSingular => format!("Matrix is not singular"),
            SRError::IsPositiveDefinite => format!("Matrix is positive definite"),
            SRError::IsPositiveSemiDefinite => format!("Matrix is positive semi-definite"),
            SRError::IsNegativeDefinite => format!("Matrix is negative definite"),
            SRError::IsNegativeSemiDefinite => format!("Matrix is negative semi-definite"),
            SRError::IsNonDefinite => format!("Matrix is non-definite"),
            // Linear systems
            SRError::LRDimensionMismatch => format!("The dimensions of LHS and RHS don't match"),
            SRError::NoSolution => format!("No solution"),
            SRError::InfiniteSolutions => format!("Infinite solutions"),
            // Arithmetic
            SRError::DivideByZero => format!("Attempt to divide by zero"),
            // Discrete numbers
            SRError::IsNotPowerOfTwo => format!("Number is not power of two"),


        }
    }
}

