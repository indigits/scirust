#![doc="Linear algebra algorithms
"]

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