/// Defines a multiplicative identity element for Self
pub trait One {
    fn one() -> Self;
}





/******************************************************
 *
 *   One implementations.
 *
 *******************************************************/



impl One for i8{

    #[inline]
    fn one() -> i8 {
        1
    }
}

impl One for i16{

    #[inline]
    fn one() -> i16 {
        1
    }
}

impl One for i32{

    #[inline]
    fn one() -> i32 {
        1
    }
}

impl One for i64{

    #[inline]
    fn one() -> i64 {
        1
    }
}


impl One for u8{

    #[inline]
    fn one() -> u8 {
        1
    }
}

impl One for u16{

    #[inline]
    fn one() -> u16 {
        1
    }
}

impl One for u32{

    #[inline]
    fn one() -> u32 {
        1
    }
}

impl One for u64{

    #[inline]
    fn one() -> u64 {
        1
    }
}

impl One for f32{

    #[inline]
    fn one() -> f32 {
        1.
    }
}

impl One for f64{

    #[inline]
    fn one() -> f64 {
        1.
    }
}

