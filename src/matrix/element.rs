// Imports
use std::fmt::Show;


/// Defines all the traits which a matrix element must support
pub trait MatrixElt : Num+Copy+Show {

}

/// Indicate that i8 fits all requirements for being a matrix element.
impl MatrixElt for i8 {
    
}
/// Indicate that i16 fits all requirements for being a matrix element.
impl MatrixElt for i16 {
    
}
/// Indicate that i32 fits all requirements for being a matrix element.
impl MatrixElt for i32 {
    
}
/// Indicate that i64 fits all requirements for being a matrix element.
impl MatrixElt for i64 {
    
}

/// Indicate that int fits all requirements for being a matrix element.
impl MatrixElt for int {
    
}


/// Indicate that u8 fits all requirements for being a matrix element.
impl MatrixElt for u8 {
    
}
/// Indicate that u16 fits all requirements for being a matrix element.
impl MatrixElt for u16 {
    
}
/// Indicate that u32 fits all requirements for being a matrix element.
impl MatrixElt for u32 {
    
}

/// Indicate that u64 fits all requirements for being a matrix element.
impl MatrixElt for u64 {
    
}

/// Indicate that uint fits all requirements for being a matrix element.
impl MatrixElt for uint {
    
}



/// Indicate that f32 fits all requirements for being a matrix element.
impl MatrixElt for f32 {
    
}


/// Indicate that f64 fits all requirements for being a matrix element.
impl MatrixElt for f64 {
    
}
