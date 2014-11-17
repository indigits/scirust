#![doc="Linear algebra algorithms
"]
pub use self::linear_system::*;
pub use self::det::*;
pub use self::lu::*;
pub use super::error::*;
pub use self::singularity::*;
pub use self::inverse::*;
pub use self::rank::*;

mod linear_system;
mod det;
mod lu;
mod singularity;
mod inverse;
mod rank;
