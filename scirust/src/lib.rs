#![allow(unused_features)]

//! A generics based scientific computing library for Rust
//!
//!
//! # Introduction
//!
//! This library aims to provide scientific computing functionality in Rust.
//! Current focus is to provide a comprehensive API with simple straight-forward
//! implementations. Different modules cover functionality 
//! covering matrices, linear algebra, signal processing and statistics.
//!
//!
//!
//! # Implementation
//! 
//! The implementation's aim to take the best advantage
//! of language features in Rust, yet avoid getting too much into
//! heavy degree of performance optimization. We do take 
//! advantage of Rust features like generic programming, 
//! immutable by default approach, type traits, iterators,
//! amongst others. This is a pure Rust implementation. 
//! There is no integration planned at the moment with
//! C libraries like BLAS or LAPACK.  
//! 
//! 
//! # Examples
//! 
//! Constructing a simple matrix: 
//! 
//!     use scirust::api::*;
//!     let a = matrix_cw_f64(2,2, &[1., 4., 2., 8.]);
//!     println!("{}", a);
//!
//! 
//! Solving a linear system of equations using Gaussian elimination method
//!
//! ```
//! # // Import scirust library components
//! # use scirust::api::*;  
//! // Construct a 2x2 matrix 
//! let a = matrix_cw_f64(2,2, &[1., 4., 2., 5.]);
//! // Print the contents of the matrix
//! println!("{}", a);
//! // A 2x1 vector.
//! let b = vector_f64(&[3.0, 6.0]);
//! // Solve the linear equation A x = b.
//! let x = GaussElimination::new(&a, &b).solve().unwrap();
//! // Print the solution vector.
//! println!("{}", x);
//! // Verify the solution
//! assert_eq!(x, vector_f64(&[-1., 2.]));
//! // Alternatively use the linear system validation algorithm.
//! let lsv = LinearSystemValidator::new(&a, &x, &b);
//! assert!(lsv.is_max_abs_val_below_threshold(1e-6));
//! ```




// #![feature(core)]
// #![feature(alloc)]
// #![feature(step_by)]
// #![feature(convert)]
// #![feature(test)]
// #![feature(associated_type_defaults)]

#[macro_use(plugin, link)] 
#[macro_use] 
extern crate log;
extern crate rand;
extern crate num;

extern crate sralgebra;
extern crate srdiscrete;
extern crate sralgorithm;
extern crate srmatrix;
extern crate srlinalg;
extern crate srsignal;
extern crate srstat;
extern crate srinterp;

// Common modules


// // Overall API of SciRust
pub mod api;
