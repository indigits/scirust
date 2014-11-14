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
    NonSquareMatrix,
    /// The matrix is not a vector
    NotAVector,


    /// The number is not a power of two
    IsNotPowerOfTwo,



    /******************************************************
     *
     *   Errors related to solving a system of linear equations
     *
     *******************************************************/

        /// There is no solution to the system of equations
    NoSolution,
    /// There are infinite solutions to the system of equations.
    InfiniteSolutions


}


impl SRError{

    /// Converts enum values to string representation
    pub fn to_string(&self) -> String {
        match *self{
            EmptyMatrix => format!("Matrix is empty"),
            DimensionsMismatch => format!("Dimensions don't match"),
            NonSquareMatrix => format!("Matrix is not square"),
            NotAVector => format!("Matrix is not a vector"),

            IsNotPowerOfTwo => format!("Number is not power of two"),


            NoSolution => format!("No solution"),
            InfiniteSolutions => format!("Infinite solutions"),
        }
    }
}

