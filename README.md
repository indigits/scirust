SciRust
==============

Scientific computing library written in 
[Rust](http://www.rust-lang.org/) programming language. 

The objective is to design a generic library which can be used as a backbone for scientific computing.

Current emphasis is less on performance and more on providing a comprehensive API.


[![Build Status](https://travis-ci.org/indigits/scirust.svg?branch=master)](https://travis-ci.org/indigits/scirust)
[![](http://meritbadge.herokuapp.com/scirust)](https://crates.io/crates/scirust)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/indigits/scirust?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)


Current areas of focus

* Fundamental algebraic structures
* Matrices
* Linear algebra
* Statistics 
* Signal processing



A discussion group is setup at
[SciRust Google Group](https://groups.google.com/forum/#!forum/scirust).

# Features

## General 

* Pure Rust implementation
* Focus on generic programming
* Extensive unit tests for all features
* Column major implementation

## Matrices

* Generic matrix class supporting various data-types 
 (u8, i8, u16, i16, ... , f32, f64, Complex32, Complex64)
* Views over parts of matrices
* Comprehensive support for operations on matrices.
* Views over sub-matrices with similar operations.
* Special support for triangular matrices.



## Linear algebra

* Solving systems of linear equations
* LDU factorization
* Rank, Determinant, Inverse



# About Rust and Building the project


If you are unfamiliar with Rust, you are recommended to go through
[The Rust Programming Language Book](https://doc.rust-lang.org/book/README.html).

The library can be built and used using 
[Cargo](http://doc.crates.io/guide.html) which is the official
dependency management and build tool for Rust.


Working with matrices requires a lot of low level code. As
a user of the library, we expect that you won't have to write
the low level code yourself. If you are reading or debugging
through the source code of the library, you would see a lot
of low level code. Good knowledge of Rust 
help you sail through them. Check out topics like: 

* [References and Borrowing](https://doc.rust-lang.org/book/references-and-borrowing.html)
* [Lifetimes](https://doc.rust-lang.org/book/lifetimes.html)
* [Unsafe code](https://doc.rust-lang.org/book/unsafe.html)
* [Raw Pointers](https://doc.rust-lang.org/book/raw-pointers.html)
* [Traits](https://doc.rust-lang.org/book/traits.html)


The library code is full of unit tests. These unit tests serve
multiple purposes

* Making sure that the functions work as advertised.
* Extensively testing those functions which use unsafe and low level
  features of Rust.
* Learning about how to use the library features.


If you haven't read already, please familiarize yourself with
[Unit Testing in Rust](https://doc.rust-lang.org/book/testing.html).
Writing unit tests will help you write better and more reliable code.




