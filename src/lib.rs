#![feature(unsafe_destructor)]
use std::mem;
use std::ptr;
use std::ops;
use std::rt::heap::{allocate, deallocate};
use std::raw::Slice as RawSlice;

// The following is needed for destroying matrix.




#[doc = "
Represents a matrix of numbers.


The numbers are stored in row-major order.
This is the standard in C-like languages.


Currently, the stride is same as the number
of columns. May be changed later on.

"]
pub struct Mat<T> {
    rows : uint,
    cols : uint, 
    ptr : *mut T
}

/// A matrix of 64-bit signed integers
pub type MatI64 = Mat<i64>;
/// A matrix of 64-bit floating point numbers.
pub type MatF64 = Mat<f64>;

/// Errors related to matrix operations
enum MatErr{
    EmptyMatrix,
    DimensionsMismatch,
}

impl<T> Mat<T> {


    pub fn new(rows: uint, cols : uint)-> Mat<T> {
        assert! (mem::size_of::<T>() != 0);
        let size = rows * cols;
        let bytes = size * mem::size_of::<T>();
        let raw = unsafe {
            allocate(bytes, mem::min_align_of::<T>())
        };
        let ptr = raw as *mut T;
        // zero out the memory
        unsafe { ptr::zero_memory(ptr, size)};
        Mat { rows : rows, cols : cols, ptr : ptr}
    }

}

impl<T> Mat<T> {
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
    /// per row stored in the memory
    pub fn stride (&self)->uint {
        self.cols
    }

    pub fn capacity(&self)-> uint {
        self.stride() * self.rows
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
    pub fn get<'a>(&'a self, r : uint, c : uint) -> &'a T  {
        let offset = r * self.stride() + c;
        unsafe {
            &*self.ptr.offset(offset as int)
        }
    }
}


impl<T> Index<uint,T> for Mat<T> {
    #[inline]
    fn index<'a>(&'a self, index: &uint) -> &'a T {
        &self.as_slice()[*index]
    }
}


impl<T> AsSlice<T> for Mat<T> {
    /// Returns a slice into `self`.
    #[inline]
    fn as_slice<'a>(&'a self) -> &'a [T] {
        unsafe { mem::transmute(RawSlice { data: self.as_ptr(), len: self.capacity() }) }
    }
}

impl<T> Mat<T>{
    #[inline]
    pub fn as_mut_slice<'a>(&'a mut self) -> &'a mut [T] {
        unsafe {
            mem::transmute(RawSlice {
                data: self.as_mut_ptr() as *const T,
                len: self.capacity(),
            })
        }
    }
}

impl<T> ops::SliceMut<uint, [T]> for Mat<T> {
    #[inline]
    fn as_mut_slice_<'a>(&'a mut self) -> &'a mut [T] {
        self.as_mut_slice()
    }

    #[inline]
    fn slice_from_or_fail_mut<'a>(&'a mut self, start: &uint) -> &'a mut [T] {
        self.as_mut_slice().slice_from_or_fail_mut(start)
    }

    #[inline]
    fn slice_to_or_fail_mut<'a>(&'a mut self, end: &uint) -> &'a mut [T] {
        self.as_mut_slice().slice_to_or_fail_mut(end)
    }
    #[inline]
    fn slice_or_fail_mut<'a>(&'a mut self, start: &uint, end: &uint) -> &'a mut [T] {
        self.as_mut_slice().slice_or_fail_mut(start, end)
    }
}



#[unsafe_destructor]
impl<T> Drop for Mat<T> {
    fn drop(&mut self) {
        // This is (and should always remain) a no-op if the fields are
        // zeroed (when moving out, because of #[unsafe_no_drop_flag]).
        if self.num_cells() != 0 {
            unsafe {
                for x in self.as_mut_slice().iter() {
                    ptr::read(x);
                }
                dealloc(self.ptr, self.capacity())
            }
        }
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




#[cfg(test)]
mod tests {

    use  super::{Mat, MatI64};

    #[test]
    fn  create_mat0(){
        let m : MatI64 = Mat::new(3, 4);
        assert_eq!(m.num_cells(), 12);
        assert_eq!(m.size(), (3, 4));
        let v : &i64 = m.get(0, 0);
        assert_eq!(*v, 0i64);
    }

}
