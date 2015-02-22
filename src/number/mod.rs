#![doc="Defines generic traits for numbers on which
SciRust library works.
"]
// Re-exporting symbols
pub use self::number::Number;
pub use self::number::num_range;
pub use self::signed::Signed;
pub use self::zero::Zero;
pub use self::one::One;
pub use self::entry::Entry;


//pub use self::complex::{Complex, Complex32, Complex64};

pub mod number;
pub mod entry;
pub mod zero;
pub mod one;
pub mod signed;
//pub mod properties;

//pub mod complex;
