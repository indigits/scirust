#![doc="Defines trait for an entry in a matrix. 
"]

// std imports
use std::fmt::Show;
use std::num::Zero;

// complex numbers
pub use external::complex::{Complex32, Complex64};


// local imports


/// Defines basic requirements for a matrix entry
pub trait Entry : Show + Zero + Clone {

}


impl Entry for i8{

}

impl Entry for i16{

}

impl Entry for i32{

}

impl Entry for i64{

}

impl Entry for int{

}

impl Entry for u8{

}

impl Entry for u16{

}

impl Entry for u32{

}

impl Entry for u64{

}

impl Entry for uint{

}

impl Entry for f32{

}

impl Entry for f64{

}

impl Entry for Complex32{

}

impl Entry for Complex64{

}




