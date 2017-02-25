#![doc="Fundamental matrix structures


Basic Matrix element: MagmaBase
Support for addition with zero: CommutativeMonoidAddPartial
Support of add, mul, 0, 1: CommutativeMonoidAddPartial+CommutativeMonoidMulPartial
Support for addition and subtraction with zero: CommutativeGroupAddPartial
Support for add, sub, mult: CommutativeRingPartial
Support for add, sub, mult, div: FieldPartial



"]

extern crate num;
extern crate rand;
extern crate sralgebra;
extern crate util;


pub mod error;
pub mod constructors;
pub mod iter;
pub mod matrix;
pub mod random;
pub mod traits;
pub mod vector;


pub mod view;
pub mod view_conversion;
pub mod view_minmax;
// pub mod triangular_matrix;

pub mod eo {
   pub mod eo_traits;
   pub mod eo_matrix;
   pub mod eo_view;
}


pub mod update{

    pub mod traits;
    pub mod matrix_updates;
}


pub mod transpose{
    pub mod traits;
    pub mod matrix_transpose;
}

pub mod extract{
    pub mod traits;
    pub mod matrix_extract;
    pub mod view_extract;
}


pub mod matrix_conversion;
pub mod matrix_minmax;

pub mod api;


#[inline]
fn mod_n (x : isize, n : isize) -> usize {
    let x = x % n;
    if x < 0{
        (x + n.abs()) as usize
    }else{
        x as usize
    }
}


#[inline]
fn cell_to_loc(stride: usize, r : usize,  c: usize)-> usize {
    (c * stride + r)
}
 
#[inline]
fn loc_to_cell(stride: usize, location : usize) -> (usize, usize){
    let c = location / stride;
    let r = location - c*stride;
    (r, c)
}
