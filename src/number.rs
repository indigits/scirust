#![doc="Defines generic traits for numbers on which
SciRust library works.
"]


// std imports


// complex numbers
pub use external::complex::{Complex32, Complex64};

// local imports
use entry::Entry;


/// Defines basic requirements for a matrix of numbers
pub trait Number : Num+Copy+Entry {

    #[inline]
    fn is_float(&self) -> bool {
        false
    }

    #[inline]
    fn is_int(&self) -> bool {
        false
    }

    #[inline]
    fn is_complex(&self) -> bool {
        false
    }

    fn is_signed(&self) -> bool;
}

/// Indicate that i8 fits all requirements for being a matrix element.
impl Number for i8 {

    #[inline]
    fn is_signed(&self) -> bool {
        true
    }
    
}
/// Indicate that i16 fits all requirements for being a matrix element.
impl Number for i16 {

    #[inline]
    fn is_signed(&self) -> bool {
        true
    }
    
    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}
/// Indicate that i32 fits all requirements for being a matrix element.
impl Number for i32 {
    
    #[inline]
    fn is_signed(&self) -> bool {
        true
    }

    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}
/// Indicate that i64 fits all requirements for being a matrix element.
impl Number for i64 {
    
    #[inline]
    fn is_signed(&self) -> bool {
        true
    }

    
    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}

/// Indicate that int fits all requirements for being a matrix element.
impl Number for int {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        true
    }

    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}


/// Indicate that u8 fits all requirements for being a matrix element.
impl Number for u8 {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        false
    }

    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}
/// Indicate that u16 fits all requirements for being a matrix element.
impl Number for u16 {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        false
    }

    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}
/// Indicate that u32 fits all requirements for being a matrix element.
impl Number for u32 {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        false
    }

    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}

/// Indicate that u64 fits all requirements for being a matrix element.
impl Number for u64 {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        false
    }

    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}

/// Indicate that uint fits all requirements for being a matrix element.
impl Number for uint {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        false
    }

    #[inline]
    fn is_int(&self) -> bool {
        true
    }

}



/// Indicate that f32 fits all requirements for being a matrix element.
impl Number for f32 {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        true
    }

    #[inline]
    fn is_float(&self) -> bool {
        true
    }

}


/// Indicate that f64 fits all requirements for being a matrix element.
impl Number for f64 {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        true
    }

    #[inline]
    fn is_float(&self) -> bool {
        true
    }

}


/// Indicate that Complex32 fits all requirements for being a matrix element.
impl Number for Complex32 {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        true
    }

    #[inline]
    fn is_complex(&self) -> bool {
        true
    }

}


/// Indicate that Complex64 fits all requirements for being a matrix element.
impl Number for Complex64 {
    
    
    #[inline]
    fn is_signed(&self) -> bool {
        true
    }

    #[inline]
    fn is_complex(&self) -> bool {
        true
    }

}


//TODO: figure out how to make this work.
pub fn describe<T: Number>(z : T){
     println!("Signed: {}", z.is_signed());
}


/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/


#[cfg(test)]
mod tests {

    use super::*;
    use std::num::Zero;



    #[test]
    fn  test_introspection_i64(){
       let z : i64 = Zero::zero();
        describe(z);
        assert!(z.is_signed());
        assert!(!z.is_float());
        assert!(!z.is_complex());
        assert!(z.is_int());
    }

    #[test]
    fn  test_introspection_f64(){
       let z : f64 = Zero::zero();
        describe(z);
        assert!(z.is_signed());
        assert!(z.is_float());
        assert!(!z.is_complex());
        assert!(!z.is_int());
    }


    #[test]
    fn  test_introspection_c64(){
       let z : Complex64 = Zero::zero();
        describe(z);
        assert!(z.is_signed());
        assert!(!z.is_float());
        assert!(z.is_complex());
        assert!(!z.is_int());
    }

}