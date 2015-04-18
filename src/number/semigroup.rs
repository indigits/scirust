

// local imports
use number::magma::{MagmaAddPartial, MagmaAdd,
    MagmaMulPartial, MagmaMul};


///////////////////////////////////////////////////////////

pub trait SemiGroupAddPartial
    : MagmaAddPartial
{

    fn prop_add_is_associative_partial(a : Self, b : Self, c : Self) -> bool {
        let ab = a.clone() + b.clone();
        let bc = b.clone() + c.clone();
        (ab) + c == a + (bc)
    }


}

impl SemiGroupAddPartial for u8 {}
impl SemiGroupAddPartial for u16  {}
impl SemiGroupAddPartial for u32  {}
impl SemiGroupAddPartial for u64  {}
impl SemiGroupAddPartial for i8   {}
impl SemiGroupAddPartial for i16  {}
impl SemiGroupAddPartial for i32  {}
impl SemiGroupAddPartial for i64  {}
impl SemiGroupAddPartial for f32  {}
impl SemiGroupAddPartial for f64  {}

///////////////////////////////////////////////////////////


pub trait SemiGroupAdd 
    : MagmaAdd
    + SemiGroupAddPartial
{
    fn prop_add_is_associative(a : Self, b : Self, c : Self) -> bool {
        let ab = a.clone() + b.clone();
        let bc = b.clone() + c.clone();
        (ab) + c == a + (bc)
    }
}


impl SemiGroupAdd for u8 {}
impl SemiGroupAdd for u16  {}
impl SemiGroupAdd for u32  {}
impl SemiGroupAdd for u64  {}
impl SemiGroupAdd for i8   {}
impl SemiGroupAdd for i16  {}
impl SemiGroupAdd for i32  {}
impl SemiGroupAdd for i64  {}
// The following won't compile.
// impl SemiGroupAdd for f32  {}
// impl SemiGroupAdd for f64  {}

///////////////////////////////////////////////////////////

pub trait SemiGroupMulPartial
    : MagmaMulPartial
{

    fn prop_mul_is_associative_partial(a : Self, b : Self, c : Self) -> bool {
        let ab = a.clone() * b.clone();
        let bc = b.clone() * c.clone();
        (ab) * c == a * (bc)
    }


}

impl SemiGroupMulPartial for u8 {}
impl SemiGroupMulPartial for u16  {}
impl SemiGroupMulPartial for u32  {}
impl SemiGroupMulPartial for u64  {}
impl SemiGroupMulPartial for i8   {}
impl SemiGroupMulPartial for i16  {}
impl SemiGroupMulPartial for i32  {}
impl SemiGroupMulPartial for i64  {}
impl SemiGroupMulPartial for f32  {}
impl SemiGroupMulPartial for f64  {}

///////////////////////////////////////////////////////////

pub trait SemiGroupMul 
    : MagmaMul
    + SemiGroupMulPartial
{

}

impl<T> SemiGroupMul for T where
    T: MagmaMul + SemiGroupMulPartial
{}

