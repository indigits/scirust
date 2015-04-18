

// std imports


// local imports
use number::zero::Zero;
use number::one::One;
use number::semigroup::{SemiGroupAddPartial, SemiGroupAdd, 
    SemiGroupMulPartial, SemiGroupMul};

///////////////////////////////////////////////////////////

pub trait MonoidAddPartial
    : SemiGroupAddPartial
    + Zero
{

}


impl<T> MonoidAddPartial for T where
    T: SemiGroupAddPartial + Zero 
{

}



///////////////////////////////////////////////////////////


pub trait MonoidAdd 
    : SemiGroupAdd
    + MonoidAddPartial
{

}


impl<T> MonoidAdd for T where
    T: SemiGroupAdd + MonoidAddPartial
{}


///////////////////////////////////////////////////////////

pub trait MonoidMulPartial
    : SemiGroupMulPartial
    + One
{

}


impl<T> MonoidMulPartial for T where
    T: SemiGroupMulPartial + One 
{

}

///////////////////////////////////////////////////////////


pub trait MonoidMul 
    : SemiGroupMul
    + MonoidMulPartial
{

}


impl<T> MonoidMul for T where
    T: SemiGroupMul + MonoidMulPartial
{}

///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use number::zero::Zero;
    use number::one::One;

    fn check_monoid_add_partial<T: MonoidAddPartial>(a: T, b: T)->T{
        let d = a + b;
        d + Zero::zero()
    }


    #[test]
    fn test_monoid_add_partial() {
        assert_eq!(check_monoid_add_partial(2i8, 3i8 ), 5);
        assert_eq!(check_monoid_add_partial(2i16, 3i16), 5);
        assert_eq!(check_monoid_add_partial(2i32, 3i32), 5);
        assert_eq!(check_monoid_add_partial(2i64, 3i64), 5);
        assert_eq!(check_monoid_add_partial(2f32, 3f32), 5f32);
        assert_eq!(check_monoid_add_partial(2f64, 3f64), 5f64);
        // The following will not compile
        // assert_eq!(check_monoid_add_partial(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_monoid_add_partial(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_monoid_add_partial(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_monoid_add_partial(2u64, 3u64, 1u64), -4);
    }


    fn check_monoid_add<T: MonoidAdd>(a: T, b: T)->T{
        let d = a + b;
        d + Zero::zero()
    }

    #[test]
    fn test_monoid_add() {
        assert_eq!(check_monoid_add(2i8, 3i8), 5);
        assert_eq!(check_monoid_add(2i16, 3i16), 5);
        assert_eq!(check_monoid_add(2i32, 3i32), 5);
        assert_eq!(check_monoid_add(2i64, 3i64), 5);
        // The following will not compile
        // assert_eq!(check_monoid_add(2f32, 3f32, 1f32), -4f32);
        // assert_eq!(check_monoid_add(2f64, 3f64, 1f64), -4f64);
        // assert_eq!(check_monoid_add(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_monoid_add(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_monoid_add(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_monoid_add(2u64, 3u64, 1u64), -4);
    }


    fn check_monoid_mul_partial<T: MonoidMulPartial>(a: T, b: T)->T{
        let d = a * b;
        d * One::one()
    }


    #[test]
    fn test_monoid_mul_partial() {
        assert_eq!(check_monoid_mul_partial(2f32, 3f32), 6.0);
        assert_eq!(check_monoid_mul_partial(2f64, 3f64), 6.0);
        // The following will not compile
        // assert_eq!(check_monoid_mul_partial(2i8, 3i8), -4);
        // assert_eq!(check_monoid_mul_partial(2i16, 3i16, 1i16), -4);
        // assert_eq!(check_monoid_mul_partial(2i32, 3i32, 1i32), -4);
        // assert_eq!(check_monoid_mul_partial(2i64, 3i64, 1i64), -4);
        // assert_eq!(check_monoid_mul_partial(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_monoid_mul_partial(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_monoid_mul_partial(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_monoid_mul_partial(2u64, 3u64, 1u64), -4);
    }



}

