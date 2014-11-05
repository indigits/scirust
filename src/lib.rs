#![feature(unsafe_destructor, globs)]

#![doc="A generics based scientific computing library for Rust


# Introduction

This library aims to provide scientific computing functionality in Rust.
Current focus is to provide a comprehensive API with simple straight-forward
implementations. Different modules cover functionality 
covering matrices, linear algebra, signal processing and statistics.



# Implementation

The implementation's aim to take the best advantage
of language features in Rust, yet avoid getting too much into
heavy degree of performance optimization. We do take 
advantage of Rust features like generic programming, 
immutable by default approach, type traits, iterators,
amongst others. This is a pure Rust implementation. 
There is no integration planned at the moment with
C libraries like BLAS or LAPACK.  


"]




pub mod discrete{
    #![doc="Discrete mathematics
    "]

    pub use self::modular::*;

    mod modular;

}

pub mod matrix {
    #![doc="Fundamental matrix structures
    "]

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
    pub mod special;
    mod view;
}


pub mod linalg {
    #![doc="Linear algebra algorithms
    "]
    pub use self::gauss_elim::*;

    mod gauss_elim;
}