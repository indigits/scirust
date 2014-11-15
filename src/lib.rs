
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
#![feature(phase)]
#[phase(plugin, link)] extern crate log;



pub mod number;
pub mod error;

pub mod discrete{
    #![doc="Discrete mathematics
    "]

    pub use self::modular::*;

    mod modular;

}

pub mod matrix {
    #![doc="Fundamental matrix structures
    "]
    pub use super::number::*;
    pub use super::error::*;
    pub use self::constructors::*;
    pub use self::iter::*;
    pub use self::matrix::*;
    pub use self::matrix_conversion::*;
    pub use self::matrix_minmax::*;
    pub use self::random::*;
    pub use self::traits::*;
    pub use self::view::*;
    pub use self::view_conversion::*;
    pub use self::view_minmax::*;
    pub use self::triangular_matrix::*;

    mod constructors;
    mod iter;
    mod matrix;
    mod matrix_conversion;
    mod matrix_minmax;
    mod random;
    pub mod traits;
    mod view;
    mod view_conversion;
    mod view_minmax;
    mod triangular_matrix;
}


pub mod linalg {
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

}


pub mod stat {
    #![doc="Statistics
    "]
}

pub mod opt {
    #![doc="Optimization
    "]
    pub mod lp {
    #![doc="Linear programming
    "]
    }
    pub mod ls {
    #![doc="Least squares
    "]
    }
    pub mod cvx {
    #![doc="Convex optimization
    "]
    }
}



pub mod signal {
    #![doc="Signal processing
    "]
    pub use self::signal::*;
    mod signal;
}


pub mod image {
#![doc="Image processing
"]
}

pub mod audio {
#![doc="Audio signal processing
"]
}


pub mod dx {
#![doc="Data exchange

This module and its submodules intend to provide
functions for exchanging scientific data in different
formats.

Initial targets are:

* HDF (Hierarchical Data Format)
* MAT (MATLAB format)


# Remarks

At the moment, the graphics capabilities 
do not exist within SciRust. In order
to create graphics, one needs to write
data generated in a SciRust program
and save it in a suitable format 
so that it can be imported in systems like
SciPy etc. The graphics facilities of those
systems can then be used. 

"]

}


pub mod external{
#![doc="External code.

This module holds external code picked up from outside
sources till a more formal solution is found.


# Contents

## complex.rs

Rust code base has a complex.rs file but it is not yet
part of Standard Library. Using the file as it is 
in the code base. When Rust provides built-in support
for complex numbers, then this library will be removed.


"]

pub mod complex;

}

// for internal use only
pub mod util{

    pub mod memory;
}



// for internal use only
mod testdata{

    #[cfg(test)]
    pub mod matrix {

        pub use self::simple::*;

        mod simple;

    }
}
