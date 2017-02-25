#![doc="Linear algebra algorithms
"]

extern crate log;
extern crate num;
extern crate sralgebra;
extern crate srmatrix;
extern crate srdiscrete;

pub mod linear_system;
pub mod det;
pub mod lu;
pub mod singularity;
pub mod inverse;
pub mod rank;

pub mod matrix{
    pub mod mat_impl;
    pub mod mat_traits;
}

pub mod testdata;
