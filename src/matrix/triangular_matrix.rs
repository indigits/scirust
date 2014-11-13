#![doc="Implements the triangular matrix data-type
"]


// std imports
use std::ptr;
use std::mem;
use std::num::{Zero, One};
use std::rt::heap::allocate;


// local imports

use matrix::element::{Number};


use matrix::error::*;

use matrix::traits::{MatrixType, Introspection, 
    MatrixBuffer, Extraction, ERO, ECO, Updates};


// complex numbers
use external::complex::Complex32;
use external::complex::Complex64;

use util;

#[doc = "
Represents a triangular square matrix of numbers.


The matrix can be either upper triangular or 
lower triangular (but not both off course).


For an upper triangular matrix, the
numbers are stored in column major order.

For a lower triangular matrix, the numbers
are stored in row major order. 

The difference in ordering is to simplify
cell index to cell offset calculations.


The half of the matrix which is full of zeros,
is not stored in memory. Thus, all calculations
are done smartly.

This design is not too flexible. The memory allocation
is done in the beginning and never changed afterwards.
Thus, the size of the matrix remains same.
"]
pub struct TriangularMatrix<T:Number> {
    /// Number of rows and columns in the matrix
    size : uint,
    /// The pointer to raw data array of the matrix
    ptr : *mut T,
    /// Indicates whether the matrix is upper or lower triangular
    ut_flag : bool
}


/// A matrix of 8-bit signed integers
pub type TriangularMatrixI8 = TriangularMatrix<i8>;
/// A matrix of 16-bit signed integers
pub type TriangularMatrixI16 = TriangularMatrix<i16>;
/// A matrix of 32-bit signed integers
pub type TriangularMatrixI32 = TriangularMatrix<i32>;
/// A matrix of 64-bit signed integers
pub type TriangularMatrixI64 = TriangularMatrix<i64>;
/// A matrix of 8-bit unsigned integers
pub type TriangularMatrixU8 = TriangularMatrix<u8>;
/// A matrix of 16-bit unsigned integers
pub type TriangularMatrixU16 = TriangularMatrix<u16>;
/// A matrix of 32-bit unsigned integers
pub type TriangularMatrixU32 = TriangularMatrix<u32>;
/// A matrix of 64-bit unsigned integers
pub type TriangularMatrixU64 = TriangularMatrix<u64>;
/// A matrix of  unsigned integers
pub type TriangularMatrixUInt = TriangularMatrix<uint>;
/// A matrix of 32-bit floating point numbers.
pub type TriangularMatrixF32 = TriangularMatrix<f32>;
/// A matrix of 64-bit floating point numbers.
pub type TriangularMatrixF64 = TriangularMatrix<f64>;
/// A matrix of 32-bit complex numbers numbers.
pub type TriangularMatrixC32 = TriangularMatrix<Complex32>;
/// A matrix of 64-bit complex numbers numbers.
pub type TriangularMatrixC64 = TriangularMatrix<Complex64>;


/// Static functions for creating  a triangular matrix
impl<T:Number> TriangularMatrix<T> {


#[doc = "Constructs a new matrix of given size (uninitialized).
"]
    pub fn new(size: uint,
        ut_flag : bool )-> TriangularMatrix<T> {
        let capacity = size * (size + 1) / 2;
        if capacity == 0 {
            // Support for empty matrices.
            // The buffer remains a null pointer.
            return TriangularMatrix{
                size : size, 
                ut_flag : ut_flag,
                ptr : ptr::null_mut()
            };
        }
        let bytes = capacity * mem::size_of::<T>();
        let raw = unsafe {
            allocate(bytes, mem::min_align_of::<T>())
        };
        let ptr = raw as *mut T;
        TriangularMatrix {size : size,
                ut_flag : ut_flag, 
                ptr : ptr}
    }    

    /// Constructs a triangular matrix of all zeros
    pub fn zeros(size: uint, 
        ut : bool)-> TriangularMatrix<T> {
        let m : TriangularMatrix<T> = TriangularMatrix::new(size, ut);
        // zero out the memory
        unsafe { ptr::zero_memory(m.ptr, m.capacity())};
        m
    }
}


/// Core methods for all matrix types
impl<T:Number> MatrixType<T> for TriangularMatrix<T> {

    /// Returns the number of rows in the matrix
    fn num_rows(&self) -> uint {
        self.size
    }

    /// Returns the number of columns in the matrix
    fn num_cols(&self) -> uint {
        self.size
    }

    /// Returns the size of matrix in an (r, c) tuple
    fn size (&self)-> (uint, uint){
        (self.size, self.size)
    }

    /// Returns the number of cells in matrix
    fn num_cells(&self)->uint {
        self.size * self.size
    }

    /// Returns if the matrix is an identity matrix
    fn is_identity(&self) -> bool {
        let o : T = One::one();
        let z  : T = Zero::zero();
        let ptr = self.ptr;
        // Check the 0, 0 element
        if unsafe{*ptr} != o {
            return false;
        }
        let mut offset = 1;
        for i in range(1, self.size){
            for _ in range (0, i){
                let v = unsafe {*ptr.offset(offset)};
                if v != z {
                    return false;
                }
                // Move on to next element in the row (l t) 
                // or column (u t).
                offset += 1;
            }
            // Check the diagonal element
            if unsafe {*ptr.offset(offset)} != o {
                return false;
            }
            // Skip the diagonal element
            offset += 1;
        }
        true
    }

    /// Returns if the matrix is a diagonal matrix
    fn is_diagonal(&self) -> bool {
        let z  : T = Zero::zero();
        let ptr = self.ptr;
        // We ignore the a[0,0] position. 
        // We start working with 2nd row (l t) or column (u t)
        let mut offset = 1;
        for i in range(1, self.size){
            for _ in range (0, i){
                let v = unsafe {*ptr.offset(offset)};
                if v != z {
                    return false;
                }
                // Move on to next element in the row (l t) 
                // or column (u t).
                offset += 1;
            }
            // Skip the diagonal element
            offset += 1;
        }
        true
    }

    /// Returns if the matrix is lower triangular 
    fn is_lt(&self) -> bool {
        ! self.ut_flag 
    }

    /// Returns if the matrix is upper triangular 
    fn is_ut(&self) -> bool {
        self.ut_flag
    }

    /// Returns if the matrix is triangular
    fn is_triangular(&self) -> bool {
        true
    }

    fn get(&self, r : uint, c : uint) -> T  {
        // These assertions help in checking matrix boundaries
        debug_assert!(r < self.size);
        debug_assert!(c < self.size);
        unsafe {
            *self.ptr.offset(self.cell_to_offset(r, c) as int)
        }
    }

    fn set(&mut self, r : uint, c : uint, value : T) {
        // These assertions help in checking matrix boundaries
        debug_assert!(r < self.size);
        debug_assert!(c < self.size);
        unsafe {
            *self.ptr.offset(self.cell_to_offset(r, c) as int) = value;
        }
    }

}

/// Introspection support
impl<T:Number> Introspection for TriangularMatrix<T> {
    /// Indicates if the matrix is a triangular matrix
    fn is_triangular_matrix_type(&self) -> bool {
        true
    }
}

/// Buffer access
impl<T:Number> MatrixBuffer<T> for TriangularMatrix<T> {
    /// Fake implementation
    fn stride (&self)->uint {
        -1 as uint
    }

    /// Returns an unsafe pointer to the matrix's 
    /// buffer.
    #[inline]
    fn as_ptr(&self)-> *const T{
        self.ptr as *const T
    }

    /// Returns a mutable unsafe pointer to
    /// the matrix's underlying buffer
    #[inline]
    fn as_mut_ptr(&mut self) -> *mut T{
        self.ptr
    }

    /// Maps a cell index to actual offset in the internal buffer
    #[inline]
    fn cell_to_offset(&self, r : uint,  c: uint)-> int {
        debug_assert!(if self.ut_flag { 
            r <= c && c < self.size 
        } 
        else {
            c <= r && r < self.size
        });
        let offset = if self.ut_flag {
            c * (c + 1) / 2 + r

        }
        else {
            r * (r + 1) /2 + c
        };
        offset as int
    } 

}


/// Main methods of a triangular matrix
impl<T:Number> TriangularMatrix<T> {

    /// Returns the capacity of the matrix 
    /// i.e. the number of elements it can hold
    pub fn capacity(&self)-> uint {
        let n = self.size;
        n * (n + 1) / 2
    }

}

#[unsafe_destructor]
impl<T:Number> Drop for TriangularMatrix<T> {
    fn drop(&mut self) {
        if self.num_cells() != 0 {
            unsafe {
                util::memory::dealloc(self.ptr, self.capacity())
            }
        }
    }
}

