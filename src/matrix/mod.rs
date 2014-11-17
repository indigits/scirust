#![doc="Fundamental matrix structures
"]
pub use super::number::*;
pub use super::error::*;
pub use self::constructors::*;
pub use self::iter::*;
pub use self::matrix::*;
pub use self::matrix_conversion::*;
pub use self::matrix_minmax::*;
pub use self::random::*;
pub use self::traits::*;
pub use self::view::*;
pub use self::view_conversion::*;
pub use self::view_minmax::*;
pub use self::triangular_matrix::*;
pub use self::vector::*;

mod constructors;
mod iter;
pub mod matrix;
mod matrix_conversion;
mod matrix_minmax;
mod random;
pub mod traits;
pub mod vector;
mod view;
mod view_conversion;
mod view_minmax;
mod triangular_matrix;
