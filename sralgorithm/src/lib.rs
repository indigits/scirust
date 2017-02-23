#![doc="Algorithms

Algorithms in this library are customized for
scientific computing needs.  They are specifically
tuned for the underlying buffer architecture
of matrices.

Much of the code also contains unsafe parts
where we deal with raw pointers directly. 

Extensive unit tests are there to validate
the correctness of the implementation.
"]


// external dependencies
extern crate log;
extern crate rand;
extern crate util;

#[macro_use]
extern crate bencher;


// public modules
pub mod sort;


