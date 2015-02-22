#![doc="Implements the triangular matrix data-type
"]


// std imports
use std::ptr;
use std::mem;
use std::fmt;
use std::rt::heap::allocate;


// local imports
use number::{Entry, Zero, One};
use number::{Number};
use matrix::matrix::Matrix;
//use matrix::error::SRError;

use matrix::traits::{Shape, NumberMatrix,
    Introspection, 
    MatrixBuffer, Extraction};


// complex numbers
//use number::Complex32;
//use number::Complex64;

use util;
use discrete::mod_n;

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
pub struct TriangularMatrix<T:Entry> {
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
//pub type TriangularMatrixC32 = TriangularMatrix<Complex32>;
/// A matrix of 64-bit complex numbers numbers.
//pub type TriangularMatrixC64 = TriangularMatrix<Complex64>;


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

    /// Constructs a matrix of all ones.
    pub fn ones(size: uint, ut : bool)-> TriangularMatrix<T> {
        let m : TriangularMatrix<T> = TriangularMatrix::new(size, ut);
        // fill with ones
        let ptr = m.ptr;
        let o : T = One::one();
        let n = m.capacity() as int;
        for i in range(0i, n){
            unsafe{
                *ptr.offset(i) = o;
            }
        }
        m
    }
}


/// Core methods for all matrix types
impl<T:Entry> Shape<T> for TriangularMatrix<T> {

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

    fn get(&self, r : uint, c : uint) -> T  {
        // These assertions help in checking matrix boundaries
        debug_assert!(r < self.size);
        debug_assert!(c < self.size);
        if (self.ut_flag && r > c) || (!self.ut_flag && r < c) {
            let v : T = Zero::zero();
            return v;
        }
        unsafe {
            // TODO : Optimize this
            self.ptr.offset(self.cell_to_offset(r, c)).as_ref().unwrap().clone()
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

/// Implementation of methods for matrices of numbers
impl<T:Number> NumberMatrix<T> for TriangularMatrix<T> {
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
    #[inline]
    fn is_lt(&self) -> bool {
        ! self.ut_flag 
    }

    /// Returns if the matrix is upper triangular 
    #[inline]
    fn is_ut(&self) -> bool {
        self.ut_flag
    }

    /// A triangular matrix is never symmetric
    #[inline]
    fn is_symmetric(&self) -> bool{
        false
    }


    /// Returns if the matrix is triangular
    #[inline]
    fn is_triangular(&self) -> bool {
        true
    }

    /// Computes the trace of the matrix
    fn trace(&self) -> T{
        if self.is_empty() {
            return Zero::zero()
        }
        let mut offset = 0i;
        let ptr = self.as_ptr();
        let mut result = unsafe{*ptr};
        for i in range(1, self.smaller_dim()){
            offset += (i + 1) as int;
            result = result + unsafe{*ptr.offset(offset)};
        }
        result
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
impl<T:Entry> MatrixBuffer<T> for TriangularMatrix<T> {

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



/// Implement extraction API for triangular matrix 
impl <T:Number> Extraction<T> for TriangularMatrix<T> {

    /// Returns the r'th row vector
    fn row(&self, r : int) -> Matrix<T> {
        // Lets ensure that the row value is mapped to
        // a value in the range [0, rows - 1]
        let r = mod_n(r, self.num_rows() as int);        
        let mut result : Matrix<T> = Matrix::new(1, self.num_cols());
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        let z : T  = Zero::zero();
        let mut dst_offset = 0i;
        let n = self.size;
        if self.ut_flag {
            // The triangle is stored in column major order.
            // We have zeros in the beginning r zeros
            for _ in range(0, r) {
                unsafe{
                    *pd.offset(dst_offset) = z;
                }
                dst_offset += 1;
            }
            for c in range(r, n){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    *pd.offset(dst_offset) = *ps.offset(src_offset);
                }
                dst_offset += 1;
            }
        }
        else {
            // the lower triangle is stored in row major order.
            // r-th row contains r + 1 entries. Rest are zero.
            for c in range(0, r + 1) {
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    *pd.offset(dst_offset) = *ps.offset(src_offset);
                }
                dst_offset += 1;
            }
            for _ in range(r + 1, n){
                unsafe{
                    *pd.offset(dst_offset) = z;
                }
                dst_offset += 1;
            }
       }
        result
    }

    /// Returns the c'th column vector
    fn col(&self, c : int) -> Matrix<T>{
        // Lets ensure that the col value is mapped to
        // a value in the range [0, cols - 1]
        let c = mod_n(c, self.num_cols() as int);        
        let mut result : Matrix<T> = Matrix::new(self.num_rows(), 1);
        let pd = result.as_mut_ptr();
        let ps = self.as_ptr();
        let z : T  = Zero::zero();
        let mut dst_offset = 0i;
        let n = self.size;
        if self.ut_flag {
            // The upper triangle is stored in column major order.
            // c-th col contains c + 1 entries. Rest are zero.
            for r in range(0, c + 1) {
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    *pd.offset(dst_offset) = *ps.offset(src_offset);
                }
                dst_offset += 1;
            }
            for _ in range(c + 1, n){
                unsafe{
                    *pd.offset(dst_offset) = z;
                }
                dst_offset += 1;
            }
        }
        else {
            // the lower triangle is stored in row major order.
            // We have zeros in the beginning c zeros
            for _ in range(0, c) {
                unsafe{
                    *pd.offset(dst_offset) = z;
                }
                dst_offset += 1;
            }
            for r in range(c, n){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    *pd.offset(dst_offset) = *ps.offset(src_offset);
                }
                dst_offset += 1;
            }
       }
        result
    }

    /// Extract a submatrix from the matrix
    /// rows can easily repeat if the number of requested rows is higher than actual rows
    /// cols can easily repeat if the number of requested cols is higher than actual cols
    fn sub_matrix(&self, start_row : int, 
        start_col : int , 
        num_rows: uint, 
        num_cols : uint) -> Matrix<T>{
        let r = mod_n(start_row, self.num_rows() as int);        
        let c = mod_n(start_col, self.num_cols() as int);
        let mut result : Matrix<T> = Matrix::new(num_rows, num_cols);
        let pd = result.as_mut_ptr();
        let mut dc = 0;
        for c in range(c, c + num_cols).map(|x | x % self.num_cols()) {
            let mut dr = 0;
            for r in range(r , r + num_rows).map(|x|  x % self.num_rows()) {
                let v = self.get(r, c);
                let dst_offset = result.cell_to_offset(dr, dc);
                unsafe{
                    *pd.offset(dst_offset) = v;
                }
                dr += 1;
            }
            dc += 1;
        }
        result
    }

    /// Returns the upper triangular part of the matrix
    fn ut_matrix(&self)->Matrix<T>{
        unimplemented!();
    }

    /// Returns the lower triangular part of the matrix
    fn lt_matrix(&self)->Matrix<T>{
        unimplemented!();
    }
}


/// Formatting of the triangular matrix on screen
impl <T:Number> fmt::Debug for TriangularMatrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // We need to find out the number of characters needed
        // to show each value.
        // maximum value
        let cap = self.capacity();
        let mut strings : Vec<String> = Vec::with_capacity(cap);
        let mut max_len : uint = 0;
        let ptr = self.ptr;
        for i in range(0i, cap as int){
            let v = unsafe {*ptr.offset(i)};
            let s = v.to_string();
            let slen = s.len();
            strings.push(s);
            if slen > max_len {
                max_len = slen;
            }
        }
        try!(write!(f, "["));
        // Here we print row by row
        let n = self.size;
        for r in range (0, n) {
           try!(write!(f, "\n  "));
            for c in range (0, n){
                if self.ut_flag {
                    if r > c {
                        for _ in range(0, max_len + 1){
                            try!(write!(f, " "));
                        }
                        try!(write!(f, "0"));
                        continue;
                    }
                }
                else {
                    if r < c {
                        for _ in range(0, max_len + 1){
                            try!(write!(f, " "));
                        }
                        try!(write!(f, "0"));
                        continue;
                    }
                }
                // This is something from within the matrix.
                let offset = self.cell_to_offset(r, c);
                let ref s = strings[offset as uint];
                let extra = max_len + 2 - s.len();
                for _ in range(0, extra){
                    try!(write!(f, " "));
                }
                try!(write!(f, "{}", s));
            }
        }
        try!(write!(f, "\n]"));
        Ok(())
    }
}

impl <T:Number> fmt::Display for TriangularMatrix<T> {
    /// Display and Debug versions are same
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}



/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/

#[cfg(test)]
mod tests {

    use super::*;
    use matrix::traits::*;
    use matrix::constructors::*;

    #[test]
    fn test_create_0(){
        let m : TriangularMatrixI64 = TriangularMatrix::new(4, true);
        println!("{}", m);
        assert_eq!(m.size(), (4, 4));
    }

    #[test]
    fn test_create_ones(){
        let n = 4;
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, true);
        println!("{}", m);
        for r in range(0, n) {
            for c in range(0, n){
                if r > c {
                    assert_eq!(m.get(r, c), 0);
                }
                else {
                    assert_eq!(m.get(r, c), 1);
                }
            }
        }
        assert!(m.is_ut());
        assert!(!m.is_lt());
        assert!(m.is_triangular());
        assert!(m.is_triangular_matrix_type());
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, false);
        println!("{}", m);
        for r in range(0, n) {
            for c in range(0, n){
                if r < c {
                    assert_eq!(m.get(r, c), 0);
                }
                else {
                    assert_eq!(m.get(r, c), 1);
                }
            }
        }
        assert_eq!(m.num_rows(), n);
        assert_eq!(m.num_cols(), n);
        assert_eq!(m.size(), (n, n));
        assert_eq!(m.num_cells(), n * n);
        assert_eq!(m.capacity(), n * (n + 1)/ 2);
        assert!(!m.is_ut());
        assert!(m.is_lt());
        assert!(m.is_triangular());
        assert!(m.is_triangular_matrix_type());
    }

    #[test]
    fn test_identity(){
        let n = 4;
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, true);
        assert!(!m.is_identity());
        let mut m : TriangularMatrixI64 = TriangularMatrix::zeros(n, true);
        assert!(!m.is_identity());
        for i in range(0, n){
            m.set(i, i, 1);
        }
        assert!(m.is_identity());
        assert!(m.is_diagonal());
    }

    #[test]
    fn test_is_diagonal(){
        let n = 4;
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, true);
        assert!(!m.is_diagonal());
        let mut m : TriangularMatrixI64 = TriangularMatrix::zeros(n, true);
        assert!(m.is_diagonal());
        for i in range(0, n){
            m.set(i, i, i as i64);
        }
        assert!(!m.is_identity());
        assert!(m.is_diagonal());
    }

    #[test]
    fn test_extract_row(){
        let n = 4;
        // upper triangular 
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, true);
        let r = m.row(0);
        assert_eq!(r, matrix_cw_i64(1,4, [1, 1, 1, 1].as_slice()));
        let r = m.row(1);
        assert_eq!(r, matrix_cw_i64(1,4, [0, 1, 1, 1].as_slice()));
        let r = m.row(2);
        assert_eq!(r, matrix_cw_i64(1,4, [0, 0, 1, 1].as_slice()));
        let r = m.row(3);
        assert_eq!(r, matrix_cw_i64(1,4, [0, 0, 0, 1].as_slice()));
        // lower triangular 
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, false);
        let r = m.row(0);
        assert_eq!(r, matrix_cw_i64(1,4, [1, 0, 0, 0].as_slice()));
        let r = m.row(1);
        assert_eq!(r, matrix_cw_i64(1,4, [1, 1, 0, 0].as_slice()));
        let r = m.row(2);
        assert_eq!(r, matrix_cw_i64(1,4, [1, 1, 1, 0].as_slice()));
        let r = m.row(3);
        assert_eq!(r, matrix_cw_i64(1,4, [1, 1, 1, 1].as_slice()));
    }

    #[test]
    fn test_extract_col(){
        let n = 4;
        // lower triangular
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, false);
        let r = m.col(0);
        assert_eq!(r, matrix_cw_i64(4,1, [1, 1, 1, 1].as_slice()));
        let r = m.col(1);
        assert_eq!(r, matrix_cw_i64(4,1, [0, 1, 1, 1].as_slice()));
        let r = m.col(2);
        assert_eq!(r, matrix_cw_i64(4,1, [0, 0, 1, 1].as_slice()));
        let r = m.col(3);
        assert_eq!(r, matrix_cw_i64(4,1, [0, 0, 0, 1].as_slice()));
        // upper triangular 
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, true);
        let r = m.col(0);
        assert_eq!(r, matrix_cw_i64(4,1, [1, 0, 0, 0].as_slice()));
        let r = m.col(1);
        assert_eq!(r, matrix_cw_i64(4,1, [1, 1, 0, 0].as_slice()));
        let r = m.col(2);
        assert_eq!(r, matrix_cw_i64(4,1, [1, 1, 1, 0].as_slice()));
        let r = m.col(3);
        assert_eq!(r, matrix_cw_i64(4,1, [1, 1, 1, 1].as_slice()));
    }

    #[test]
    fn test_extract_sub_matrix(){
        let n = 4;
        // lower triangular
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, false);
        let m2 = m.sub_matrix(0, 0, 2, 2);
        assert_eq!(m2, matrix_rw_i64(2,2, [
            1, 0,
            1, 1].as_slice()));
    }

    #[test]
    fn test_trace(){
        let n = 4;
        // lower triangular
        let m : TriangularMatrixI64 = TriangularMatrix::ones(n, false);
        assert_eq!(m.trace(), 4);
    }

}



