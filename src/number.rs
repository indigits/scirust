#![doc="Defines generic traits for numbers on which
SciRust library works.
"]


// Imports
use std::fmt::Show;


// complex numbers
use external::complex::{Complex32, Complex64};


/// Defines all the traits which a matrix element must support
pub trait Number : Num+Copy+Show {

    fn is_float() -> bool {
        false
    }

    fn is_int() -> bool {
        false
    }

    fn is_complex() -> bool {
        false
    }

    fn is_signed() -> bool;
}

/// Indicate that i8 fits all requirements for being a matrix element.
impl Number for i8 {

    #[inline]
    fn is_signed() -> bool {
        true
    }
    
}
/// Indicate that i16 fits all requirements for being a matrix element.
impl Number for i16 {

    #[inline]
    fn is_signed() -> bool {
        true
    }
    
    fn is_int() -> bool {
        true
    }

}
/// Indicate that i32 fits all requirements for being a matrix element.
impl Number for i32 {
    
    #[inline]
    fn is_signed() -> bool {
        true
    }

    fn is_int() -> bool {
        true
    }

}
/// Indicate that i64 fits all requirements for being a matrix element.
impl Number for i64 {
    
    #[inline]
    fn is_signed() -> bool {
        true
    }

    
    fn is_int() -> bool {
        true
    }

}

/// Indicate that int fits all requirements for being a matrix element.
impl Number for int {
    
    
    #[inline]
    fn is_signed() -> bool {
        true
    }

    fn is_int() -> bool {
        true
    }

}


/// Indicate that u8 fits all requirements for being a matrix element.
impl Number for u8 {
    
    
    #[inline]
    fn is_signed() -> bool {
        false
    }

    fn is_int() -> bool {
        true
    }

}
/// Indicate that u16 fits all requirements for being a matrix element.
impl Number for u16 {
    
    
    #[inline]
    fn is_signed() -> bool {
        false
    }

    fn is_int() -> bool {
        true
    }

}
/// Indicate that u32 fits all requirements for being a matrix element.
impl Number for u32 {
    
    
    #[inline]
    fn is_signed() -> bool {
        false
    }

    fn is_int() -> bool {
        true
    }

}

/// Indicate that u64 fits all requirements for being a matrix element.
impl Number for u64 {
    
    
    #[inline]
    fn is_signed() -> bool {
        false
    }

    fn is_int() -> bool {
        true
    }

}

/// Indicate that uint fits all requirements for being a matrix element.
impl Number for uint {
    
    
    #[inline]
    fn is_signed() -> bool {
        false
    }

    fn is_int() -> bool {
        true
    }

}



/// Indicate that f32 fits all requirements for being a matrix element.
impl Number for f32 {
    
    
    #[inline]
    fn is_signed() -> bool {
        true
    }

    fn is_float() -> bool {
        true
    }

}


/// Indicate that f64 fits all requirements for being a matrix element.
impl Number for f64 {
    
    
    #[inline]
    fn is_signed() -> bool {
        true
    }

    fn is_float() -> bool {
        true
    }

}


/// Indicate that Complex32 fits all requirements for being a matrix element.
impl Number for Complex32 {
    
    
    #[inline]
    fn is_signed() -> bool {
        true
    }

}


/// Indicate that Complex64 fits all requirements for being a matrix element.
impl Number for Complex64 {
    
    
    #[inline]
    fn is_signed() -> bool {
        true
    }

}


//TODO: figure out how to make this work.
//pub fn describe<T: Number>(){
    //println!("Signed: {}", T.is_signed());
//}
