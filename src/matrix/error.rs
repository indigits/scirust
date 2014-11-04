/// Errors related to matrix operations
#[deriving(Show)]
pub enum MatrixError{
    /// The matrix is empty
    EmptyMatrix,
    /// The dimensions of two matrices mismatch
    DimensionsMismatch,
    /// The matrix is not a square matrix
    NonSquareMatrix,
    /// The object is not a vector
    NotAVector,
}
