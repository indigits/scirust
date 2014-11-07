#![doc="Some example matrices and other objects which are 
useful for testing purposes.
"]

use matrix::*;

#[doc="Matrix properties:
 
determinant = -13
"]  
pub fn square_0() -> MatrixF64 {
    matrix_rw_f64(3,3,[
        0. , 1., 2.,
        3., -1., 0.,
        1., -2., 1.
        ])
} 


#[doc="Matrix properties:
 
determinant = 6
"]pub fn square_1() -> MatrixF64 {
    matrix_rw_f64(3,3,[
        1. , -1., 1.,
        1., 1., 1.,
        1., 2., 4.
        ])
} 
