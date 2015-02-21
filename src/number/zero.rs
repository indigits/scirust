

/// Defines an additive identity element for Self
pub trait Zero {
    fn zero() -> Self;
    fn is_zero(&self) -> bool;
}


/******************************************************
 *
 *   Zero implementations.
 *
 *******************************************************/

impl Zero for i8{

    #[inline]
    fn zero() -> i8{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0
    }
}

impl Zero for i16{
    #[inline]
    fn zero() -> i16{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0
    }
}

impl Zero for i32{
    #[inline]
    fn zero() -> i32{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0
    }
}

impl Zero for i64{
    #[inline]
    fn zero() -> i64{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0
    }
}


impl Zero for u8{
    #[inline]
    fn zero() -> u8{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0
    }
}

impl Zero for u16{
    #[inline]
    fn zero() -> u16{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0
    }
}

impl Zero for u32{
    #[inline]
    fn zero() -> u32{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0
    }
}

impl Zero for u64{
    #[inline]
    fn zero() -> u64{
        0
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0
    }
}

impl Zero for f32{
    #[inline]
    fn zero() -> f32{
        0.
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0.
    }
}

impl Zero for f64{
    #[inline]
    fn zero() -> f64{
        0.
    }

    #[inline]
    fn is_zero(&self) -> bool{
        *self == 0.
    }
}

