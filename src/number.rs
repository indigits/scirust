#![doc="Defines generic traits for numbers on which
SciRust library works.
"]


// Imports
use std::fmt::Show;


// complex numbers
use external::complex::{Complex32, Complex64};


/// Defines all the traits which a matrix element must support
pub trait Number : Num+Copy+Show {

}

/// Indicate that i8 fits all requirements for being a matrix element.
impl Number for i8 {
    
}
/// Indicate that i16 fits all requirements for being a matrix element.
impl Number for i16 {
    
}
/// Indicate that i32 fits all requirements for being a matrix element.
impl Number for i32 {
    
}
/// Indicate that i64 fits all requirements for being a matrix element.
impl Number for i64 {
    
}

/// Indicate that int fits all requirements for being a matrix element.
impl Number for int {
    
}


/// Indicate that u8 fits all requirements for being a matrix element.
impl Number for u8 {
    
}
/// Indicate that u16 fits all requirements for being a matrix element.
impl Number for u16 {
    
}
/// Indicate that u32 fits all requirements for being a matrix element.
impl Number for u32 {
    
}

/// Indicate that u64 fits all requirements for being a matrix element.
impl Number for u64 {
    
}

/// Indicate that uint fits all requirements for being a matrix element.
impl Number for uint {
    
}



/// Indicate that f32 fits all requirements for being a matrix element.
impl Number for f32 {
    
}


/// Indicate that f64 fits all requirements for being a matrix element.
impl Number for f64 {
    
}


/// Indicate that Complex32 fits all requirements for being a matrix element.
impl Number for Complex32 {
    
}


/// Indicate that Complex64 fits all requirements for being a matrix element.
impl Number for Complex64 {
    
}
