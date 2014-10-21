use std::mem;
use std::ptr;
use std::ops;
use std::cmp;
use std::fmt;
use std::num::One;
use std::iter::Iterator;
use std::rt::heap::{allocate, deallocate};
use std::raw::Slice as RawSlice;

// The following is needed for destroying matrix.




#[doc = "
Represents a matrix of numbers.


The numbers are stored in column-major order.
This is the standard in Fortran and MATLAB.


Currently, the stride is same as the number
of rows. May be changed later on.

For a row vector, the stride  = 1.

"]
pub struct Mat<T:Ord+Copy> {
    rows : uint,
    cols : uint, 
    ptr : *mut T
}

/// A matrix of 64-bit signed integers
pub type MatI64 = Mat<i64>;
/// A matrix of 64-bit floating point numbers.
pub type MatF64 = Mat<f64>;

/// Errors related to matrix operations
#[deriving(Show)]
pub enum MatErr{
    //EmptyMatrix,
    DimensionsMismatch,
}

/// Static functions for creating  a matrix
impl<T:Num+Ord+Copy> Mat<T> {


    pub fn new(rows: uint, cols : uint)-> Mat<T> {
        assert! (mem::size_of::<T>() != 0);
        let size = rows * cols;
        let bytes = size * mem::size_of::<T>();
        let raw = unsafe {
            allocate(bytes, mem::min_align_of::<T>())
        };
        let ptr = raw as *mut T;
        Mat { rows : rows, cols : cols, ptr : ptr}
    }

    /// Constructs a matrix of all zeros
    pub fn zeros(rows: uint, cols : uint)-> Mat<T> {
        let m : Mat<T> = Mat::new(rows, cols);
        // zero out the memory
        unsafe { ptr::zero_memory(m.ptr, m.capacity())};
        m
    }


    pub fn ones(rows: uint, cols : uint)-> Mat<T> {
        let m : Mat<T> = Mat::new(rows, cols);
        // fill with ones
        let ptr = m.ptr;
        let o : T = One::one();
        unsafe {
            for c in range(0, cols){
                for r in range (0, rows){
                    let offset = m.cell_to_offset(r, c);
                    let p  = ptr.offset(offset as int);
                    *p = o;
                }
            } 
        }
        m
    }


    pub fn from_slice(rows: uint, cols : uint, values: &[T]) -> Mat<T>{
        let mut mat : Mat<T> = Mat::new(rows, cols);
        let n_cells = mat.num_cells();
        let stride = mat.stride();
        // get a mutable slice from m
        {
            let dst_slice = mat.as_mut_slice();
            let n_values = values.len();
            // The number of entries we can copy
            let n_fill = cmp::min(n_values, n_cells);
            let mut n = 0;
            let mut offset_src = 0;
            let mut offset_dst = 0;
            for _ in range(0, cols){
                for r in range(0, rows){
                    n+=1;
                    dst_slice[offset_dst + r] = values[offset_src + r];
                    if n == n_fill {
                        break;
                    }
                }
                if n == n_fill {
                    break;
                }
                offset_src += rows;
                offset_dst += stride;
            }
        }
        // return
        mat
    }

    pub fn from_iter< A : Iterator<T>>(rows: uint, cols : uint, mut iter: A) -> Mat<T>{
        let mut mat : Mat<T> = Mat::new(rows, cols);
        let stride = mat.stride();
        // get a mutable slice from m
        {
            let dst_slice = mat.as_mut_slice();
            let mut offset_dst = 0;
            'outer: for _ in range(0, cols){
                for r in range(0, rows){
                    let next_val = iter.next();
                    match next_val{
                        Some(val) => dst_slice[offset_dst + r] = val,
                        None => break 'outer
                    };
                }
                offset_dst += stride;
            }
        }
        // return
        mat
    }
}

/// Main methods of a matrix
impl<T:Num+Ord+Copy> Mat<T> {
    pub fn num_rows(&self) -> uint {
        self.rows
    }

    pub fn num_cols(&self) -> uint {
        self.cols
    }

    pub fn size (&self)-> (uint, uint){
        (self.rows, self.cols)
    }
    pub fn num_cells(&self)->uint {
        self.rows * self.cols
    }
    /// Returns the number of actual memory elements 
    /// per column stored in the memory
    pub fn stride (&self)->uint {
        self.rows
    }

    pub fn capacity(&self)-> uint {
        self.stride() * self.cols
    }

    /// Returns an unsafe pointer to the matrix's 
    /// buffer.
    #[inline]
    pub fn as_ptr(&self)-> *const T{
        self.ptr as *const T
    }

    /// Returns a mutable unsafe pointer to
    /// the matrix's underlying buffer
    #[inline]
    pub fn as_mut_ptr(&mut self) -> *mut T{
        self.ptr
    }

    #[inline]
    pub fn get(&self, r : uint, c : uint) -> T  {
        unsafe {
            *self.ptr.offset(self.cell_to_offset(r, c) as int)
        }
    }

    /// Converts an index to cell address (row, column)
    #[inline]
    pub fn index_to_cell(&self, index : uint) -> (uint, uint){
        let c = index / self.rows;
        let r = index - c*self.rows;
        (r, c)
    }

    /// Converts a cell address to an index (r, c) to index
    #[inline]
    pub fn cell_to_index(&self, r : uint,  c: uint) -> uint{
        c * self.rows + r
    }

    /// Converts the matrix to vector from standard library
    pub fn to_std_vec(&self) -> Vec<T> {
        let mut vec: Vec<T> = Vec::with_capacity(self.num_cells());
        // We iterate over elements in matrix and push in the vector
        let s = self.as_slice_();
        for c in range(0, self.cols){
            for r in range (0, self.rows){
                let offset = self.cell_to_offset(r, c);
                vec.push(s[offset]);
            }
        } 
        vec
    }
}


impl<T:Num+Ord+Copy> Index<uint,T> for Mat<T> {
    #[inline]
    fn index<'a>(&'a self, index: &uint) -> &'a T {
        // The matrix is column major order
        let (r, c) = self.index_to_cell(*index);
        let offset = c * self.stride() + r;
        unsafe {
            &*self.ptr.offset(offset as int)
        }
    }
}



impl<T:Num+Ord+Copy> Mat<T>{
    /// This function is for internal use only.
    #[inline]
    fn as_mut_slice<'a>(&'a mut self) -> &'a mut [T] {
        unsafe {
            mem::transmute(RawSlice {
                data: self.as_mut_ptr() as *const T,
                len: self.capacity(),
            })
        }
    }
}



impl <T:Num+Ord+Copy> Clone for Mat<T> {

    fn clone(&self )-> Mat<T> {
        let m : Mat<T> = Mat::new(self.rows, self.cols);
        unsafe{
            ptr::copy_memory(m.ptr, self.ptr as *const T, self.capacity());
        }
        m
    }
}

impl <T:Num+Ord+Copy+fmt::Show> fmt::Show for Mat<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = self.as_slice_();
        try!(write!(f, "["));
        // Here we print row by row
        for r in range (0, self.rows) {
           try!(write!(f, "\n  "));
            for c in range (0, self.cols){
                let offset = c * self.stride() + r;
                try!(write!(f, "{} ", s[offset]));
            }
        }
        try!(write!(f, "\n]\n"));
        Ok(())
    }
}

/// Matrix addition support
impl<T:Num+Ord+Copy> ops::Add<Mat<T>, Mat<T>> for Mat<T> {
    fn add(&self, _rhs: &Mat<T>) -> Mat<T> {
        // Validate dimensions are same.
        if self.size() != _rhs.size(){
            fail!(DimensionsMismatch.to_string());
        }
        let result : Mat<T> = Mat::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = _rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in range(0, n){
                let i = i_ as int;
                *pc.offset(i) = *pa.offset(i) + *pb.offset(i);
            }
        }
        result
    }
}

#[unsafe_destructor]
impl<T:Num+Ord+Copy> Drop for Mat<T> {
    fn drop(&mut self) {
        if self.num_cells() != 0 {
            unsafe {
                dealloc(self.ptr, self.capacity())
            }
        }
    }
}


/******************************************************
 *
 *   Private implementation of Mat
 *
 *******************************************************/

impl<T:Num+Ord+Copy> Mat<T> {
    /// Returns a slice into `self`.
    #[inline]
    fn as_slice_<'a>(&'a self) -> &'a [T] {
        unsafe { mem::transmute(RawSlice { data: self.as_ptr(), len: self.capacity() }) }
    }

    /// Maps a cell index to actual offset in the internal buffer
    #[inline]
    fn cell_to_offset(&self, r : uint,  c: uint)-> uint {
        c * self.stride() + r
    } 
}


/***
Private helper functions follow
*/

#[inline]
unsafe fn dealloc<T>(ptr: *mut T, len: uint) {
    if mem::size_of::<T>() != 0 {
        deallocate(ptr as *mut u8,
                   len * mem::size_of::<T>(),
                   mem::min_align_of::<T>())
    }
}


/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/


#[cfg(test)]
mod tests {

    use  super::{Mat, MatI64};

    #[test]
    fn  create_mat0(){
        let m : MatI64 = Mat::new(3, 4);
        assert_eq!(m.num_cells(), 12);
        assert_eq!(m.size(), (3, 4));
        let v : i64 = m.get(0, 0);
        assert_eq!(v, 0i64);
    }


    #[test]
    fn test_from_iter0(){
        let m : MatI64 = Mat::from_iter(4, 4, range(1, 20));
        let b: Vec<i64> = range(1, 17).collect();
        assert!(m.as_slice_() == b.as_slice_());
        let m : MatI64 = Mat::from_iter(4, 4, range(1, 16));
        assert_eq!(m.get(0, 0), 1);
        assert_eq!(m.get(2, 2), 11);
        let mut b: Vec<i64> = range(1, 16).collect();
        b.push(0);
        assert!(m.as_slice_() == b.as_slice_());
    }

    #[test]
    fn test_index0(){
        let m : MatI64 = Mat::from_iter(4, 4, range(1, 20));
        let x = m[4];
        assert_eq!(x, 5);
    }

    #[test]
    fn test_cell_index_mapping(){
        let rows = 20;
        let cols = 10;
        let m : MatI64 = Mat::new(rows, cols);
        assert_eq!(m.index_to_cell(3 + 2*rows), (3, 2));
        assert_eq!(m.cell_to_index(3, 2), 3 + 2*rows);
        assert_eq!(m.cell_to_index(5, 7), 145);
        assert_eq!(m.index_to_cell(145), (5, 7));
    }

    #[test]
    fn test_to_std_vec(){
        let m : MatI64 = Mat::from_iter(4, 3, range(0, 12));
        let v1 = m.to_std_vec();
        let v2 : Vec<i64> = range(0, 12).collect();
        assert_eq!(v1, v2);
    }

    #[test]
    fn test_ones(){
        let m : MatI64 = Mat::ones(4, 2);
        let v = vec![1i64, 1, 1, 1, 1, 1, 1, 1];
        assert_eq!(m.to_std_vec(), v);
    }

    #[test]
    fn test_sum(){
        let m : MatI64 = Mat::ones(4, 2);
        let m2 = m  + m;
        let v = vec![2i64, 2, 2, 2, 2, 2, 2, 2];
        assert_eq!(m2.to_std_vec(), v);
    }

    #[test]
    #[should_fail]
    fn test_sum_fail(){
        let m1 : MatI64 = Mat::ones(4, 2);
        let m2 : MatI64 = Mat::ones(3, 2);
        m1 + m2;
    }


}
