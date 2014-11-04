#![feature(unsafe_destructor, globs)]

pub mod discrete{

    pub use self::modular::*;

    mod modular;

}

pub mod matrix {
    pub use self::element::*;
    pub use self::error::*;
    pub use self::iter::*;
    pub use self::matrix::*;
    pub use self::random::*;
    pub use self::view::*;

    mod element;
    mod error;
    mod iter;
    mod matrix;
    mod random;
    mod view;
}

pub mod linalg {
    pub use self::gauss_elim::*;

    mod gauss_elim;
}