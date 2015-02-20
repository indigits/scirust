#![doc="Defines trait for an entry in a matrix. 
"]

// std imports
use std::fmt::Debug;



// local imports
use number::zero::Zero;

/// Defines basic requirements for a matrix entry
pub trait Entry : Debug + Clone + Zero {

}




/******************************************************
 *
 *   Entry implementations.
 *
 *******************************************************/



impl Entry for i8{

}

impl Entry for i16{

}

impl Entry for i32{

}

impl Entry for i64{

}

impl Entry for isize{

}

impl Entry for u8{

}

impl Entry for u16{

}

impl Entry for u32{

}

impl Entry for u64{

}

impl Entry for usize{

}

impl Entry for f32{

}

impl Entry for f64{

}





