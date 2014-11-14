#![doc="Defines errors related to SciRust library
"]


/// Enum for errors related to SciRust library
#[deriving(Show)]
pub enum SRError{
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


}
