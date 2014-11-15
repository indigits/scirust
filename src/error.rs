#![doc="Defines errors related to SciRust library
"]


// std imports


// local imports



/// Enum for errors related to SciRust library
#[deriving(Show)]
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
    /// The matrix is not a square matrix
    IsNotSquareMatrix,
    /// The matrix is not a vector
    NotAVector,
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
     *   Discrete numbers related stuff
     *
     *******************************************************/

    /// The number is not a power of two
    IsNotPowerOfTwo,

}


impl SRError{

    /// Converts enum values to string representation
    pub fn to_string(&self) -> String {
        match *self{
            //  Matrices
            EmptyMatrix => format!("Matrix is empty"),
            DimensionsMismatch => format!("Dimensions don't match"),
            IsNotSquareMatrix => format!("Matrix is not square"),
            NotAVector => format!("Matrix is not a vector"),
            IsNotFullRankMatrix => format!("Matrix is not full rank"),
            IsFullRankMatrix => format!("Matrix is full rank"),
            IsSingular => format!("Matrix is singular"),
            IsNonSingular => format!("Matrix is not singular"),
            IsPositiveDefinite => format!("Matrix is positive definite"),
            IsPositiveSemiDefinite => format!("Matrix is positive semi-definite"),
            IsNegativeDefinite => format!("Matrix is negative definite"),
            IsNegativeSemiDefinite => format!("Matrix is negative semi-definite"),
            IsNonDefinite => format!("Matrix is non-definite"),
            // Linear systems
            LRDimensionMismatch => format!("The dimensions of LHS and RHS don't match"),
            NoSolution => format!("No solution"),
            InfiniteSolutions => format!("Infinite solutions"),
            // Discrete numbers
            IsNotPowerOfTwo => format!("Number is not power of two"),


        }
    }
}

