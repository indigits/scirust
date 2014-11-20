#![doc="Defines generic traits for numbers on which
SciRust library works.
"]


// std imports
use std::num::{Int, SignedInt, Float};

// complex numbers
pub use external::complex::{Complex32, Complex64};

// local imports
use entry::{Entry, Zero, One};


/// Defines basic requirements for a matrix of numbers
pub trait Number : Entry
    + Clone
    + Copy
    + Add<Self,Self>
    + Sub<Self, Self> 
    + Mul<Self, Self> 
    + Div<Self, Self>
    + PartialEq
    + One
    + Zero{

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

    fn power(&self, n : uint)-> Self{
        if n == 0 {
            return One::one();
        }
        let mut result = *self;
        for _ in range(0, n - 1){
            result = result.mul(self);
        }
        result
    }
}


pub trait Signed: Number + Neg<Self> {
    fn abs_val(&self) -> Self;
    fn signum(&self) -> Self;
    fn is_positive(&self) -> bool;
    fn is_negative(&self) -> bool;
}



/******************************************************
 *      i8 implementation
 ******************************************************/


/// Indicate that i8 fits all requirements for being a matrix element.
impl Number for i8 {

    #[inline]
    fn is_signed(&self) -> bool {
        true
    }
    
}

impl Signed for i8 {

    #[inline]
    fn abs_val(&self) -> i8{
        self.abs()
    }

    #[inline]
    fn signum(&self) -> i8{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}


/******************************************************
 *      i16 implementation
 ******************************************************/



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

impl Signed for i16 {

    #[inline]
    fn abs_val(&self) -> i16{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i16{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}
    
/******************************************************
 *      i32 implementation
 ******************************************************/


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

impl Signed for i32 {

    #[inline]
    fn abs_val(&self) -> i32{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i32{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}
    

/******************************************************
 *      i64 implementation
 ******************************************************/



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

impl Signed for i64 {

    #[inline]
    fn abs_val(&self) -> i64{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> i64{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}
    
/******************************************************
 *      int implementation
 ******************************************************/

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

impl Signed for int {

    #[inline]
    fn abs_val(&self) -> int{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> int{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}
    

/******************************************************
 *      u8 implementation
 ******************************************************/


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

/******************************************************
 *      u16 implementation
 ******************************************************/


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

/******************************************************
 *      u32 implementation
 ******************************************************/


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

/******************************************************
 *      u64 implementation
 ******************************************************/



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

/******************************************************
 *      uint implementation
 ******************************************************/


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


/******************************************************
 *      f32 implementation
 ******************************************************/

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

impl Signed for f32 {

    #[inline]
    fn abs_val(&self) -> f32{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> f32{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}


/******************************************************
 *      f64 implementation
 ******************************************************/

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

impl Signed for f64 {

    #[inline]
    fn abs_val(&self) -> f64{
        self.abs()
    }
    
    #[inline]
    fn signum(&self) -> f64{
        self.signum()
    }
    
    #[inline]
    fn is_positive(&self) -> bool{
        self.is_positive()
    }
    
    #[inline]
    fn is_negative(&self) -> bool{
        self.is_negative()
    }
}

/******************************************************
 *      Complex32 implementation
 ******************************************************/

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

/******************************************************
 *      Complex64 implementation
 ******************************************************/

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
 *   Number range
 *
 *******************************************************/

/// An iterator over the range [start, stop)
#[deriving(Clone)]
pub struct NumRange<A> {
    state: A,
    stop: A,
    one: A
}


/// Returns an iterator over the given range [start, stop) (that is, starting
/// at start (inclusive), and ending at stop (exclusive)).
#[inline]
pub fn num_range<A: Number + PartialOrd + One>(start: A, stop: A) -> NumRange<A> {
    NumRange{state: start, stop: stop, one: One::one()}
}

impl<A: Number + PartialOrd + One + ToPrimitive> Iterator<A> for NumRange<A> {
    #[inline]
    fn next(&mut self) -> Option<A> {
        if self.state < self.stop {
            let result = self.state.clone();
            self.state = self.state + self.one;
            Some(result)
        } else {
            None
        }
    }

    #[inline]
    fn size_hint(&self) -> (uint, Option<uint>) {
        // This first checks if the elements are representable as i64. If they aren't, try u64 (to
        // handle cases like range(huge, huger)). We don't use uint/int because the difference of
        // the i64/u64 might lie within their range.
        let bound = match self.state.to_i64() {
            Some(a) => {
                let sz = self.stop.to_i64().map(|b| b.checked_sub(a));
                match sz {
                    Some(Some(bound)) => bound.to_uint(),
                    _ => None,
                }
            },
            None => match self.state.to_u64() {
                Some(a) => {
                    let sz = self.stop.to_u64().map(|b| b.checked_sub(a));
                    match sz {
                        Some(Some(bound)) => bound.to_uint(),
                        _ => None
                    }
                },
                None => None
            }
        };

        match bound {
            Some(b) => (b, Some(b)),
            // Standard fallback for unbounded/unrepresentable bounds
            None => (0, None)
        }
    }
}



/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/


#[cfg(test)]
mod tests {

    use super::*;
    use entry::Zero;



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