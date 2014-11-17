
// std imports
use std::mem;

// local imports
use number::Number;
use super::eo_traits::{ERO, ECO};
use matrix::view::MatrixView;
use matrix::traits::{Shape, MatrixBuffer};


/// Implementation of Elementary row operations.
impl<'a, T:Number> ERO<T> for MatrixView<'a, T> {


    /// Row scaling by a factor and adding to another row.
    /// r_i = r_i + k * r_j
    /// The j-th row can be outside the view also.
    /// This is the row relative to the start of the view.
    #[inline]
    fn ero_scale_add(&mut self, 
        i :  uint, 
        j :  int, 
        scale : T
        )-> &mut MatrixView<'a, T> {
        debug_assert! (i  < self.num_rows());
        let m = self.matrix();
        // Compute j-th row in m (by doing offset)
        let j = j + (self.start_row() as int);
        debug_assert! (j  >= 0);
        let j = j as uint;
        debug_assert!(j < m.num_rows());
        let ptr = m.as_ptr();
        // I am allowing modification of the underlying buffer
        let ptr : *mut T = unsafe { mem::transmute(ptr) };
        let sc = self.start_col();
        // Compute initial offsets
        let mut offset_a = self.cell_to_offset(i, 0);
        let mut offset_b = m.cell_to_offset(j, sc);
        let stride_a = self.stride() as int;
        let stride_b = m.stride() as int;
        for _ in range(0, self.num_cols()){
            unsafe {
                let va = *ptr.offset(offset_a);
                let vb = *ptr.offset(offset_b);
                *ptr.offset(offset_a) = va + scale * vb;
            }
            // Update offsets
            offset_a += stride_a; 
            offset_b += stride_b;
        }
        self
    }
}

/// Implementation of Elementary column operations.
impl<'a, T:Number> ECO<T> for MatrixView<'a, T> {
    /// Column scaling by a factor and adding to another column.
    /// c_i = c_i + k * c_j
    /// The j-th column can be outside the view also.
    /// This is the column relative to the start of the view.
    #[inline]
    fn eco_scale_add(&mut self, 
        i :  uint, 
        j :  int, 
        scale : T
        )-> &mut MatrixView<'a, T> {
        debug_assert! (i  < self.num_cols());
        let m = self.matrix();
        // Compute j-th column in m (by doing offset)
        let j = j + (self.start_col() as int);
        debug_assert! (j  >= 0);
        let j = j as uint;
        debug_assert!(j < m.num_cols());
        let ptr = m.as_ptr();
        // I am allowing modification of the underlying buffer
        let ptr : *mut T = unsafe { mem::transmute(ptr) };
        let sr = self.start_row();
        // Compute initial offsets
        let mut offset_a = self.cell_to_offset(0, i);
        let mut offset_b = m.cell_to_offset(sr, j);
        for _ in range(0, self.num_rows()){
            unsafe {
                let va = *ptr.offset(offset_a);
                let vb = *ptr.offset(offset_b);
                *ptr.offset(offset_a) = va + scale * vb;
            }
            // Update offsets
            offset_a += 1; 
            offset_b += 1;
        }
        self
    }
}


/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/


#[cfg(test)]
mod test{
    //use super::*;
}


/******************************************************
 *
 *   Bench marks
 *
 *******************************************************/


#[cfg(test)]
mod bench{
    //extern crate test;
    //use self::test::Bencher;
    //use super::*;
}


