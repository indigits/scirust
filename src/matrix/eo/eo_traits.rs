/******************************************************
 *
 *   Elementary row / column operations.
 *
 *******************************************************/

/// std imports
use std::ptr;

/// local imports
use discrete::mod_n;
use number::Number;
use matrix::traits::{Shape, MatrixBuffer, Strided};



/// Elementary row operations on a matrix
pub trait ERO<T:Number> : Shape<T>+MatrixBuffer<T> + Strided {

    /// Row switching.
    #[inline]
    fn ero_switch(&mut self, 
        i :  usize,
        j : usize
        )-> &mut Self {
        debug_assert! (i  < self.num_rows());
        debug_assert! (j  < self.num_rows());
        let ptr = self.as_mut_ptr();
        let mut offset_a = self.cell_to_offset(i, 0);
        let mut offset_b = self.cell_to_offset(j, 0);
        let stride = self.stride() as isize;
        for _ in range(0, self.num_cols()){
            unsafe {
                ptr::swap(ptr.offset(offset_a), ptr.offset(offset_b));
                offset_a += stride;
                offset_b += stride;
            }
        }
        self
    }

    /// Row scaling by a factor.
    #[inline]
    fn ero_scale(&mut self, 
        r :  usize, 
        scale : T
        )-> &mut Self {
        debug_assert! (r  < self.num_rows());
        let ptr = self.as_mut_ptr();
        for c in range(0, self.num_cols()){
            let offset = self.cell_to_offset(r, c);
            unsafe {
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale * v;
            }
        }
        self
    }

    /// Row scaling by a factor over a slice of the row. [start, end)
    #[inline]
    fn ero_scale_slice(&mut self, 
        r :  usize, 
        scale : T,
        start : usize,
        end : usize,
        )-> &mut Self {
        debug_assert! (r  < self.num_rows());
        debug_assert!(start <= end);
        debug_assert!(end <= self.num_cols());
        let ptr = self.as_mut_ptr();
        for c in range(start, end){
            let offset = self.cell_to_offset(r, c);
            unsafe {
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale * v;
            }
        }
        self
    }

    /// Row scaling by a factor and adding to another row.
    /// r_i = r_i + k * r_j
    #[inline]
    fn ero_scale_add(&mut self, 
        i :  usize, 
        j :  isize, 
        scale : T
        )-> &mut Self {
        debug_assert! (i  < self.num_rows());
        debug_assert! (j  < self.num_rows() as isize);
        let j = if j < 0{
            mod_n(j, self.num_rows() as isize)
        }
        else {
            j as usize
        };
        let ptr = self.as_mut_ptr();
        for c in range(0, self.num_cols()){
            let offset_a = self.cell_to_offset(i, c);
            let offset_b = self.cell_to_offset(j, c);
            unsafe {
                let va = *ptr.offset(offset_a);
                let vb = *ptr.offset(offset_b);
                *ptr.offset(offset_a) = va + scale * vb;
            }
        }
        self
    }

}


/// Elementary column operations on a matrix
pub trait ECO<T:Number> : Shape<T>+MatrixBuffer<T> + Strided {

    /// Column switching.
    #[inline]
    fn eco_switch(&mut self, 
        i :  usize,
        j : usize
        )-> &mut Self {
        debug_assert! (i  < self.num_cols());
        debug_assert! (j  < self.num_cols());
        let ptr = self.as_mut_ptr();
        let mut offset_a = self.cell_to_offset(0, i);
        let mut offset_b = self.cell_to_offset(0, j);
        for _ in range(0, self.num_rows()){
            unsafe {
                ptr::swap(ptr.offset(offset_a), ptr.offset(offset_b));
                offset_a += 1;
                offset_b += 1;
            }
        }
        self
    }

    /// Column scaling by a factor.
    #[inline]
    fn eco_scale(&mut self, 
        c :  usize, 
        scale : T
        )-> &mut Self {
        debug_assert! (c  < self.num_cols());
        let ptr = self.as_mut_ptr();
        let mut offset = self.cell_to_offset(0, c);
        for _ in range(0, self.num_rows()){
            unsafe {
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale * v;
                offset += 1;
            }
        }
        self
    }

    /// Column scaling by a factor over a slice
    #[inline]
    fn eco_scale_slice(&mut self, 
        c :  usize, 
        scale : T,
        start : usize,
        end : usize,
        )-> &mut Self {
        debug_assert! (c  < self.num_cols());
        debug_assert!(start < end);
        debug_assert!(end <= self.num_rows());
        let ptr = self.as_mut_ptr();
        let mut offset = self.cell_to_offset(start, c);
        for _ in range(start, end){
            unsafe {
                let v = *ptr.offset(offset);
                *ptr.offset(offset) = scale * v;
                offset += 1;
            }
        }
        self
    }

    /// Column scaling by a factor and adding to another column.
    /// c_i = c_i + k * c_j
    #[inline]
    fn eco_scale_add(&mut self, 
        i :  usize, 
        j :  isize, 
        scale : T
        )-> &mut Self {
        debug_assert! (i  < self.num_cols());
        debug_assert! (j  < self.num_cols() as isize);
        let j = if j < 0{
            mod_n(j, self.num_cols() as isize)
        }
        else {
            j as usize
        };
        let ptr = self.as_mut_ptr();
        let mut offset_a = self.cell_to_offset(0, i);
        let mut offset_b = self.cell_to_offset(0, j);
        for _ in range(0, self.num_rows()){
            unsafe {
                let va = *ptr.offset(offset_a);
                let vb = *ptr.offset(offset_b);
                *ptr.offset(offset_a) = va + scale * vb;
                offset_a += 1;
                offset_b += 1;
            }
        }
        self
    }

}
