#![feature(unsafe_destructor, globs)]
pub use  self::matrix::*;
pub use  self::matrand::*;
pub use  self::matiter::*;
pub use  self::matview::*;

mod matelt;
mod materr;
mod matrix;
mod matrand;
mod matiter;
mod matview;


pub mod discrete{

    pub use self::modular::*;

    mod modular;

}

pub mod linalg {
    pub use self::gauss_elim::*;

    mod gauss_elim;
}