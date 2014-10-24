// Imports
use std::fmt::Show;


/// Defines all the traits which a matrix element must support
pub trait MatElt : Num+Copy+Show {

}

/// Indicate that i8 fits all requirements for being a matrix element.
impl MatElt for i8 {
    
}
/// Indicate that i16 fits all requirements for being a matrix element.
impl MatElt for i16 {
    
}
/// Indicate that i32 fits all requirements for being a matrix element.
impl MatElt for i32 {
    
}
/// Indicate that i64 fits all requirements for being a matrix element.
impl MatElt for i64 {
    
}



/// Indicate that u8 fits all requirements for being a matrix element.
impl MatElt for u8 {
    
}
/// Indicate that u16 fits all requirements for being a matrix element.
impl MatElt for u16 {
    
}
/// Indicate that u32 fits all requirements for being a matrix element.
impl MatElt for u32 {
    
}

/// Indicate that u64 fits all requirements for being a matrix element.
impl MatElt for u64 {
    
}



/// Indicate that f32 fits all requirements for being a matrix element.
impl MatElt for f32 {
    
}


/// Indicate that f64 fits all requirements for being a matrix element.
impl MatElt for f64 {
    
}
