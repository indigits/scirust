#![feature(unsafe_destructor, globs)]
pub use  self::matrix::*;
pub use  self::matrand::*;
mod discrete;
mod matelt;
mod materr;
mod matrix;
mod matrand;
