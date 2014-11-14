SciRust
==============

Scientific computing library written in 
[Rust](http://www.rust-lang.org/) programming language. 

The objective is to design a generic library which can be used as a backbone for scientific computing.

Current emphasis is less on performance and more on providing a comprehensive API.

Current areas of focus

* Matrices
* Linear algebra
* Statistics 
* Signal processing


For more details and examples, 
see the [API Guide](http://indigits.github.io/scirust).

If you are unfamiliar with Rust, you are recommended to go through
[The Rust Guide](http://doc.rust-lang.org/guide.html).

The library can be built and used using 
[Cargo](http://doc.crates.io/guide.html) which is the official
dependency management and build tool for Rust.


Working with matrices requires a lot of low level code. As
a user of the library, we expect that you won't have to write
the low level code yourself. If you are reading or debugging
through the source code of the library, you would see a lot
of low level code. Familiarity with following guides will
help you sail through them.

* [References and Lifetime Guide](http://doc.rust-lang.org/guide-lifetimes.html)
* [Writing low level and unsafe code in Rust](http://doc.rust-lang.org/guide-unsafe.html)
* [Rust Pointers Guide](http://doc.rust-lang.org/guide-pointers.html)


The library code is full of unit tests. These unit tests serve
multiple purposes

* Making sure that the functions work as advertised.
* Extensively testing those functions which use unsafe and low level
  features of Rust.
* Learning about how to use the library features.


If you haven't read already, please familiarize yourself with
[Rust Testing Guide](http://doc.rust-lang.org/guide-testing.html).
Writing unit tests will help you write better and more reliable code.




