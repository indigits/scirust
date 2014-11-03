#![feature(unsafe_destructor, globs)]
pub use  self::matrix::*;
pub use  self::matrand::*;
pub use  self::matiter::*;
pub use  self::matview::*;
mod discrete;
mod matelt;
mod materr;
mod matrix;
mod matrand;
mod matiter;
mod matview;