#![feature(unsafe_destructor, globs)]
pub use  self::matrix::{Mat, MatI64, MatF64};
mod discrete;
mod matelt;
mod materr;
mod matrix;
