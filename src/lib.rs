
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

#![feature(unsafe_destructor, globs)]
#![feature(macro_rules)]
extern crate log;
extern crate collections;
//extern crate quickcheck_macros;
//extern crate quickcheck;


// Common modules
pub mod number;
pub mod error;

// Main libraries
pub mod external;
pub mod discrete;
pub mod matrix;
pub mod linalg;
pub mod signal;
// pub mod dx;
pub mod stat;
pub mod alg;
// pub mod opt;
pub mod image {
#![doc="Image processing
"]
}

pub mod audio {
#![doc="Audio signal processing
"]
}

// for internal use only
pub mod util;
mod testdata;

// Overall API of SciRust
pub mod api;
