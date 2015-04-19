

// std imports


// local imports
use algebra::zero::Zero;
use algebra::one::One;
use algebra::quasigroup::{QuasiGroupAddPartial, QuasiGroupAdd, 
    QuasiGroupMulPartial, QuasiGroupMul};

///////////////////////////////////////////////////////////

pub trait LoopAddPartial
    : QuasiGroupAddPartial
    + Zero
{

}


impl<T> LoopAddPartial for T where
    T: QuasiGroupAddPartial + Zero 
{

}



///////////////////////////////////////////////////////////


pub trait LoopAdd 
    : QuasiGroupAdd
    + LoopAddPartial
{

}


impl<T> LoopAdd for T where
    T: QuasiGroupAdd + LoopAddPartial
{}


///////////////////////////////////////////////////////////

pub trait LoopMulPartial
    : QuasiGroupMulPartial
    + One
{

}


impl<T> LoopMulPartial for T where
    T: QuasiGroupMulPartial + One 
{

}

///////////////////////////////////////////////////////////


pub trait LoopMul 
    : QuasiGroupMul
    + LoopMulPartial
{

}


impl<T> LoopMul for T where
    T: QuasiGroupMul + LoopMulPartial
{}

///////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::zero::Zero;
    use algebra::one::One;

    fn check_loop_add_partial<T: LoopAddPartial>(a: T, b: T, c : T)->T{
        let d = a + b;
        let e = d - c;
        -e + Zero::zero()
    }


    #[test]
    fn test_loop_add_partial() {
        assert_eq!(check_loop_add_partial(2i8, 3i8, 1i8), -4);
        assert_eq!(check_loop_add_partial(2i16, 3i16, 1i16), -4);
        assert_eq!(check_loop_add_partial(2i32, 3i32, 1i32), -4);
        assert_eq!(check_loop_add_partial(2i64, 3i64, 1i64), -4);
        assert_eq!(check_loop_add_partial(2f32, 3f32, 1f32), -4f32);
        assert_eq!(check_loop_add_partial(2f64, 3f64, 1f64), -4f64);
        // The following will not compile
        // assert_eq!(check_loop_add_partial(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_loop_add_partial(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_loop_add_partial(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_loop_add_partial(2u64, 3u64, 1u64), -4);
    }


    fn check_loop_add<T: LoopAdd>(a: T, b: T, c : T)->T{
        let d = a + b;
        let e = d - c;
        -e + Zero::zero()
    }

    #[test]
    fn test_loop_add() {
        assert_eq!(check_loop_add(2i8, 3i8, 1i8), -4);
        assert_eq!(check_loop_add(2i16, 3i16, 1i16), -4);
        assert_eq!(check_loop_add(2i32, 3i32, 1i32), -4);
        assert_eq!(check_loop_add(2i64, 3i64, 1i64), -4);
        // The following will not compile
        // assert_eq!(check_loop_add(2f32, 3f32, 1f32), -4f32);
        // assert_eq!(check_loop_add(2f64, 3f64, 1f64), -4f64);
        // assert_eq!(check_loop_add(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_loop_add(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_loop_add(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_loop_add(2u64, 3u64, 1u64), -4);
    }


    fn check_loop_mul_partial<T: LoopMulPartial>(a: T, b: T, c : T)->T{
        let d = a * b;
        let e = d / c;
        e.recip() * One::one()
    }


    #[test]
    fn test_loop_mul_partial() {
        assert_eq!(check_loop_mul_partial(2f32, 3f32, 3f32), 0.5);
        assert_eq!(check_loop_mul_partial(2f64, 3f64, 3f64), 0.5);
        // The following will not compile
        // assert_eq!(check_loop_mul_partial(2i8, 3i8, 1i8), -4);
        // assert_eq!(check_loop_mul_partial(2i16, 3i16, 1i16), -4);
        // assert_eq!(check_loop_mul_partial(2i32, 3i32, 1i32), -4);
        // assert_eq!(check_loop_mul_partial(2i64, 3i64, 1i64), -4);
        // assert_eq!(check_loop_mul_partial(2u8, 3u8, 1u8), -4);
        // assert_eq!(check_loop_mul_partial(2u16, 3u16, 1u16), -4);
        // assert_eq!(check_loop_mul_partial(2u32, 3u32, 1u32), -4);
        // assert_eq!(check_loop_mul_partial(2u64, 3u64, 1u64), -4);
    }



}

