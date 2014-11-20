#![doc="SciRust API
"]

pub use error::SRError;

// Matrix library
pub use matrix::*;
pub use matrix::traits::*;
pub use matrix::matrix::*;
pub use matrix::constructors::*;
pub use matrix::view::*;
pub use matrix::eo::eo_traits::*;
// pub use super::number::*;
// pub use super::error::SRError;
// pub use self::constructors::*;
// pub use self::iter::*;
// pub use self::matrix::*;
// pub use self::matrix_conversion::*;
// pub use self::matrix_minmax::*;
// pub use self::random::*;
// pub use self::traits::*;
// pub use self::view::*;
// pub use self::view_conversion::*;
// pub use self::view_minmax::*;
// pub use self::triangular_matrix::*;
// pub use self::vector::*;
//pub use self::eo::*;


// Linear algebra library
pub use linalg::linear_system::*;
pub use linalg::lu::*;
pub use linalg::singularity::*;
pub use linalg::inverse::*;
pub use linalg::rank::*;


// Statistics library
pub use stat::*;
pub use stat::moments::*;

