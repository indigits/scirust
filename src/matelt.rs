use std::fmt::Show;
/// Defines all the traits which a matrix element must support
pub trait MatElt : Num+Copy+Show {

}

/// Indicate that i64 fits all requirements for being a matrix element.
impl MatElt for i64 {
    
}

/// Indicate that f64 fits all requirements for being a matrix element.
impl MatElt for f64 {
    
}


/// Indicate that u8 fits all requirements for being a matrix element.
impl MatElt for u8 {
    
}

