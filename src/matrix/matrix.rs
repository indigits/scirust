#![doc="Provides the basic Matrix data type
"]

// std imports
use std::mem;
use std::ptr;
use std::ops;
use std::cmp;
use std::fmt;
use num::{Float};
use std::iter::Iterator;
use std::rt::heap::{allocate, deallocate};
use std::ops::{Index};
use std::raw::Slice;

// external imports
use num::traits::{Zero, One, Signed};

// local imports

use discrete::{mod_n};
use algebra::structure::{MagmaBase, 
    CommutativeMonoidAddPartial, 
    CommutativeMonoidMulPartial,
    QuasiGroupAddPartial,
    CommutativeGroupAddPartial,
    FieldPartial};
use error::SRError;
use matrix::iter::*;
use matrix::view::MatrixView;
use matrix::traits::{Shape, NumberMatrix,
    Strided,
    StridedNumberMatrix,
    StridedFloatMatrix,
    Introspection, 
    MatrixBuffer,
    Search};

// complex numbers
//use algebra::Complex32;
//use algebra::Complex64;

use util;


#[doc = "
Represents a matrix of numbers.


The numbers are stored in column-major order.
This is the standard in Fortran and MATLAB.

"]
pub struct Matrix<T:MagmaBase> {
    /// Number of rows in the matrix
    rows : usize,
    /// Number of columns in the matrix
    cols : usize, 
    /// Number of allocated rows 
    xrows : usize, 
    /// Number of allocated columns
    xcols : usize,
    /// The pointer to raw data array of the matrix
    ptr : *mut T
}

/// A matrix of 8-bit signed integers
pub type MatrixI8 = Matrix<i8>;
/// A matrix of 16-bit signed integers
pub type MatrixI16 = Matrix<i16>;
/// A matrix of 32-bit signed integers
pub type MatrixI32 = Matrix<i32>;
/// A matrix of 64-bit signed integers
pub type MatrixI64 = Matrix<i64>;
/// A matrix of 8-bit unsigned integers
pub type MatrixU8 = Matrix<u8>;
/// A matrix of 16-bit unsigned integers
pub type MatrixU16 = Matrix<u16>;
/// A matrix of 32-bit unsigned integers
pub type MatrixU32 = Matrix<u32>;
/// A matrix of 64-bit unsigned integers
pub type MatrixU64 = Matrix<u64>;
/// A matrix of 32-bit floating point numbers.
pub type MatrixF32 = Matrix<f32>;
/// A matrix of 64-bit floating point numbers.
pub type MatrixF64 = Matrix<f64>;
/// A matrix of 32-bit complex numbers numbers.
//pub type MatrixC32 = Matrix<Complex32>;
/// A matrix of 64-bit complex numbers numbers.
//pub type MatrixC64 = Matrix<Complex64>;




/// Static functions for creating  a matrix
impl<T:MagmaBase> Matrix<T> {

#[doc = "Constructs a new matrix of given size (uninitialized).

# Remarks 

The contents of the matrix are not initialized. Hence, 
it doesn't make sense to use this function liberally.
Still the function is internally useful since
different constructor functions need to initialize
the matrix differently.
"]
    pub fn new(rows: usize, cols : usize)-> Matrix<T> {
        debug_assert! (mem::size_of::<T>() != 0);
        let xrows = rows.next_power_of_two();
        let xcols = cols.next_power_of_two();
        let capacity = xrows *  xcols;

        // Support for empty  matrices
        if capacity == 0{
            // We do not allocate any memory for the buffer.
            // We leave it as a NULL pointer.
            return Matrix { rows : rows, 
                cols : cols,
                xrows : xrows,
                xcols : xcols, 
                ptr : ptr::null_mut()
            };
        }
        let bytes = capacity * mem::size_of::<T>();
        let raw = unsafe {
            allocate(bytes, mem::min_align_of::<T>())
        };
        let ptr = raw as *mut T;
        Matrix { rows : rows, 
                cols : cols,
                xrows : xrows,
                xcols : xcols, 
                ptr : ptr}
    }

}

/// Static functions for creating  a matrix of numbers
impl<T:CommutativeMonoidAddPartial> Matrix<T> {
    /// Constructs a scalar matrix
    pub fn from_scalar (scalar : T) -> Matrix <T>{
        let m : Matrix<T> = Matrix::new(1, 1);
        unsafe {*m.ptr = scalar;}
        m
    }

    /// Constructs a matrix of all zeros
    pub fn zeros(rows: usize, cols : usize)-> Matrix<T> {
        let m : Matrix<T> = Matrix::new(rows, cols);
        let ptr = m.ptr;
        // zero out the memory
        let z : T = Zero::zero();
        unsafe {
            for c in 0..cols{
                for r in 0..rows{
                    let offset = m.cell_to_offset(r, c);
                    let p  = ptr.offset(offset as isize);
                    *p = z;
                }
            } 
        }
        m
    }

    #[doc = "Constructs a matrix from a slice of data reading
    data in column wise order.
    "]
    pub fn from_slice_cw(rows: usize, cols : usize, values: &[T]) -> Matrix<T>{
        let mut mat : Matrix<T> = Matrix::new(rows, cols);
        // stride of new matrix
        let stride = mat.stride();
        // get a mutable slice from m
        {
            let dst_slice = mat.as_mut_slice();
            // The number of entries we can copy
            let n_values = values.len();
            let z : T = Zero::zero();
            let mut n = 0;
            let mut offset_dst = 0;
            for _ in 0..cols{
                for r in 0..rows{
                    let v = if n < n_values {
                        values[n]
                    }else{
                        z
                    };
                    dst_slice[offset_dst + r] = v;
                    n+=1;
                }
                offset_dst += stride;
            }
        }
        // return
        mat
    }

    #[doc = "Constructs a matrix from a slice of data reading data in row wise order.

# Remarks 

In source code, when we are constructing matrices from slices, a
slice looks easier to read in row wise order. Thus, this
function should be more useful in constructing matrices
by hand.
"]
    pub fn from_slice_rw(rows: usize, cols : usize, values: &[T]) -> Matrix<T>{
        let mat : Matrix<T> = Matrix::new(rows, cols);
        // get a mutable slice from m
        {
            let ptr = mat.ptr;
            // The number of entries we can copy
            let n_values = values.len();
            let z : T = Zero::zero();
            let mut n = 0;
            for r in 0..rows{
                for c in 0..cols{
                    let v = if n < n_values {
                        values[n]
                    }else{
                        z
                    };
                    let dst_offset = mat.cell_to_offset(r,c);
                    unsafe{
                        *ptr.offset(dst_offset) = v;
                    }
                    n+=1;
                }
            }
        }
        // return
        mat
    }

    pub fn from_iter_cw< A : Iterator<Item=T>>(rows: usize, cols : usize, mut iter: A) -> Matrix<T>{
        let mut mat : Matrix<T> = Matrix::new(rows, cols);
        let stride = mat.stride();
        // get a mutable slice from m
        {
            let dst_slice = mat.as_mut_slice();
            let mut offset_dst = 0;
            let z : T = Zero::zero();
            let mut completed_columns  = 0;
            'outer: for _ in 0..cols{
                for r in 0..rows{
                    let next_val = iter.next();
                    match next_val{
                        Some(val) => dst_slice[offset_dst + r] = val,
                        None => {
                            // Finish this column with zeros
                            for _ in r..rows{
                                dst_slice[offset_dst + r] = z;
                            }
                            completed_columns += 1;
                            offset_dst += stride;
                            break 'outer
                        }
                    };
                }
                completed_columns += 1;
                offset_dst += stride;
            }
            if completed_columns < cols {
                // We  need to fill remaining columns with zeros
                for _ in completed_columns..cols{
                    for r in 0..rows{
                        dst_slice[offset_dst + r] = z;
                    }
                    completed_columns += 1;
                    offset_dst += stride;
                }
            }
        }
        // return
        mat
    }


    /// Builds a matrix from an iterator reading numbers in a 
    /// row-wise order
    pub fn from_iter_rw< A : Iterator<Item=T>>(rows: usize, cols : usize, 
        iter: A) -> Matrix<T>{
        let m : Matrix<T> = Matrix::new(rows, cols);
        let ptr = m.ptr;
        let z : T = Zero::zero();
        let mut r = 0;
        let mut c = 0;
        let nc = m.num_cols();
        let nr = m.num_rows();
        for v in iter {
            if c == nc {
                c = 0;
                r = r + 1;
            }
            if r == nr {
                break;
            }
            let dst_offset = m.cell_to_offset(r, c);
            unsafe{
                *ptr.offset(dst_offset) = v;
            }
            c += 1;
        }
        loop {
            if c == nc {
                c = 0;
                r = r + 1;
            }
            if r == nr {
                break;
            }
            let dst_offset = m.cell_to_offset(r, c);
            unsafe{
                *ptr.offset(dst_offset) = z;
            }
            c += 1;
        }
        // return
        m
    }

    /// Construct a diagonal matrix from a vector
    pub fn diag_from_vec(v : &Matrix<T>) -> Matrix<T>{
        if !v.is_vector(){
            panic!(SRError::IsNotAVector.to_string());
        }
        let n = v.num_cells();
        let m : Matrix<T> = Matrix::zeros(n, n);
        let src = v.ptr;
        let dst = m.ptr;
        // Copy the elements of v in the vector
        for r in 0..n{
            let offset = m.cell_to_offset(r, r);
            unsafe{
                *dst.offset(offset) =  *src.offset(r as isize);
            }
        }
        m
    }

}

/// Static functions for creating  a matrix of numbers
impl<T:CommutativeMonoidAddPartial+One> Matrix<T> {

    /// Constructs a matrix of all ones.
    pub fn ones(rows: usize, cols : usize)-> Matrix<T> {
        let m : Matrix<T> = Matrix::new(rows, cols);
        // fill with ones
        let ptr = m.ptr;
        let o : T = One::one();
        unsafe {
            for c in 0..cols{
                for r in 0..rows{
                    let offset = m.cell_to_offset(r, c);
                    let p  = ptr.offset(offset as isize);
                    *p = o;
                }
            } 
        }
        m
    }


    /// Constructs an identity matrix
    pub fn identity(rows: usize, cols : usize) -> Matrix<T> {
        let m : Matrix<T> = Matrix::zeros(rows, cols);
        // fill with ones
        let ptr = m.ptr;
        let one : T = One::one();
        let n = cmp::min(rows, cols);
        for i in 0..n{
            let offset = m.cell_to_offset(i, i);
            unsafe{
                *ptr.offset(offset) = one;
            }
        }
        m
    }

    /// Constructs a unit vector
    /// (1, 0, 0), (0, 1, 0), (0, 0, 1), etc.
    pub fn unit_vector( length : usize, dim : usize) -> Matrix<T> {
        let mut m : Matrix<T> = Matrix::zeros(length, 1);
        m.set(dim, 0, One::one());
        m
    }

}


/// Core methods for all matrix types
impl<T:MagmaBase> Shape<T> for Matrix<T> {

    /// Returns the number of rows in the matrix
    fn num_rows(&self) -> usize {
        self.rows
    }

    /// Returns the number of columns in the matrix
    fn num_cols(&self) -> usize {
        self.cols
    }

    /// Returns the size of matrix in an (r, c) tuple
    fn size (&self)-> (usize, usize){
        (self.rows, self.cols)
    }

    /// Returns the number of cells in matrix
    fn num_cells(&self)->usize {
        self.rows * self.cols
    }


    fn get(&self, r : usize, c : usize) -> T  {
        // These assertions help in checking matrix boundaries
        debug_assert!(r < self.rows);
        debug_assert!(c < self.cols);
        let offset = self.cell_to_offset(r, c);
        unsafe {
            let v = self.ptr.offset(offset);
            // TODO : Optimize this
            v.as_ref().unwrap().clone()
        }
    }

    fn set(&mut self, r : usize, c : usize, value : T) {
        // These assertions help in checking matrix boundaries
        debug_assert!(r < self.rows);
        debug_assert!(c < self.cols);
        unsafe {
            *self.ptr.offset(self.cell_to_offset(r, c)) = value;
        }
    }
}

/// Methods available to number matrices
impl<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial> NumberMatrix<T> for Matrix<T> {

    /// Returns if the matrix is an identity matrix
    fn is_identity(&self) -> bool {
        let o : T = One::one();
        let z  : T = Zero::zero();
        let ptr = self.ptr;
        for c in 0..self.cols{
            for r in 0..self.rows{
                let offset = self.cell_to_offset(r, c);
                let v = unsafe {*ptr.offset(offset)};
                if r == c {
                    if v != o {
                        return false;
                    }
                }else if v != z {
                    return false;
                }

            }
        } 
        true
    }

    /// Returns if the matrix is a diagonal matrix
    fn is_diagonal(&self) -> bool {
        let z  : T = Zero::zero();
        let ptr = self.ptr;
        for c in 0..self.cols{
            for r in 0..self.rows{
                if r != c {
                    let offset = self.cell_to_offset(r, c);
                    let v = unsafe {*ptr.offset(offset)};
                    if v != z {
                        return false;
                    }
                }
            }
        } 
        true
    }

    /// Returns if the matrix is lower triangular 
    fn is_lt(&self) -> bool {
        let z  : T = Zero::zero();
        let ptr = self.ptr;
        for c in 0..self.cols{
            for r in 0..c{
                let offset = self.cell_to_offset(r, c);
                let v = unsafe {*ptr.offset(offset)};
                if v != z {
                    return false;
                }
            }
        } 
        true
    }

    /// Returns if the matrix is upper triangular 
    fn is_ut(&self) -> bool {
        let z  : T = Zero::zero();
        let ptr = self.ptr;
        for c in 0..self.cols{
            for r in c+1..self.rows{
                let offset = self.cell_to_offset(r, c);
                let v = unsafe {*ptr.offset(offset)};
                if v != z {
                    return false;
                }
            }
        } 
        true
    }

    /// Returns if the matrix is symmetric
    fn is_symmetric(&self) -> bool{
        if !self.is_square(){
            // Only square matrices can be symmetric
            return false;
        }
        // size of the square matrix
        let n = self.num_rows();
        for i in 0..n{
            for j in (i + 1)..n{
                if self.get(i, j) != self.get(j, i) {
                    return false;
                }
            }
        }
        true
    }

    /// Returns the trace of the matrix
    fn trace(&self) -> T{
        if self.is_empty() {
            return Zero::zero()
        }
        let mut result = self.get(0, 0);
        let stride = self.stride() as isize;
        let mut offset = stride;
        let ptr = self.as_ptr();
        for i in 1..self.smaller_dim(){
            result = result + unsafe{*ptr.offset(offset + i as isize)};
            offset += stride;
        }
        result
    }
}

/// Introspection support
impl<T> Introspection for Matrix<T> {
    /// This is a standard matrix object
    fn is_standard_matrix_type(&self) -> bool {
        true
    }
}

/// Strided buffer
impl<T:MagmaBase> Strided for Matrix<T> {

    /// Returns the number of actual memory elements 
    /// per column stored in the memory
    fn stride (&self)->usize {
        self.xrows
    }

}

impl<T:FieldPartial> StridedNumberMatrix<T> for Matrix<T> {
}

impl<T:FieldPartial+Float> StridedFloatMatrix<T> for Matrix<T> {
}



/// Buffer access
impl<T:MagmaBase> MatrixBuffer<T> for Matrix<T> {

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
    fn cell_to_offset(&self, r : usize,  c: usize)-> isize {
        (c * self.stride() + r) as isize
    } 

}

/// Main methods of a matrix
impl<T:MagmaBase> Matrix<T> {

    /// Returns the capacity of the matrix 
    /// i.e. the number of elements it can hold
    pub fn capacity(&self)-> usize {
        self.xrows * self.xcols
    }



}


/// Functions to construct new matrices out of a matrix and other conversions
impl<T:CommutativeMonoidAddPartial+One> Matrix<T> {



    /// Returns an iterator over a specific row of matrix
    pub fn row_iter(&self, r : isize) -> RowIterator<T>{
        let r = mod_n(r, self.rows as isize);        
        // Lets find the offset of the begging of the row
        let offset = self.cell_to_offset(r, 0);
        let iter : RowIterator<T> = RowIterator::new(self.cols,
            self.stride(), unsafe {self.ptr.offset(offset)} as *const T);
        iter
    }

    /// Returns an iterator over a specific column of the matrix
    pub fn col_iter(&self, c : isize) -> ColIterator<T>{
        let c = mod_n(c, self.cols as isize);        
        // Lets find the offset of the begging of the column
        let offset = self.cell_to_offset(0, c);
        let iter : ColIterator<T> = ColIterator::new(self.rows,
            unsafe {self.ptr.offset(offset)} as *const T);
        iter
    }

    /// Returns an iterator over all cells  of the matrix
    pub fn cell_iter(&self) -> CellIterator<T>{
        let iter : CellIterator<T> = CellIterator::new(
            self.rows, self.cols, self.stride(),
            self.ptr as *const T);
        iter
    }


    // Repeats this matrix in both horizontal and vertical directions 
    pub fn repeat_matrix(&self, num_rows : usize, num_cols : usize) -> Matrix<T> {
        let rows = self.rows * num_rows;
        let cols = self.cols * num_cols;
        let result : Matrix<T> = Matrix::new(rows, cols);
        let pd = result.ptr;
        let ps = self.ptr;
        for bc in 0..num_cols{
            let bc_start = bc * self.cols;
            for br in 0..num_rows{
                let br_start =  br * self.rows;
                for c in 0..self.cols {
                    for r in 0..self.rows{
                        let src_offset = self.cell_to_offset(r, c);
                        let dst_offset = result.cell_to_offset(br_start + r, bc_start + c);
                        unsafe{
                            *pd.offset(dst_offset) = *ps.offset(src_offset);
                        }
                    }
                }

            }
        }
        result
    }

    /// Extracts the primary diagonal from the matrix as a vector
    pub fn diagonal_vector(&self) -> Matrix<T> {
        let m  = self.smaller_dim();
        let result : Matrix<T> = Matrix::new(m, 1);
        let src = self.ptr;
        let dst = result.ptr;
        for i in 0..m{
            let offset = self.cell_to_offset(i, i);
            unsafe{
                *dst.offset(i as isize) = *src.offset(offset);
            } 
        }
        result        
    }

    /// Extracts the primary diagonal from the matrix as a matrix of same size
    pub fn diagonal_matrix(&self) -> Matrix<T> {
        let m  = cmp::min(self.rows, self.cols);
        let result : Matrix<T> = Matrix::zeros(self.rows, self.cols);
        let src = self.ptr;
        let dst = result.ptr;
        for i in 0..m{
            let offset = self.cell_to_offset(i, i);
            unsafe{
                *dst.offset(offset) = *src.offset(offset);
            } 
        }
        result        
    }

    /// Returns the upper triangular part of the matrix as a new matrix
    pub fn ut(&self)->Matrix<T>{
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let src = self.ptr;
        let dst = result.ptr;
        let z  : T = Zero::zero();
        for c in 0..self.cols{
            for r in 0..(c+1){
                let offset = self.cell_to_offset(r, c);
                unsafe {*dst.offset(offset) = *src.offset(offset);}
            }
            for r in (c+1)..self.rows{
                let offset = self.cell_to_offset(r, c);
                unsafe {*dst.offset(offset) = z;}
            }
        }
        result
    }

    /// Returns the lower triangular part of the matrix as a new matrix
    pub fn lt(&self)->Matrix<T>{
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let src = self.ptr;
        let dst = result.ptr;
        let z  : T = Zero::zero();
        for c in 0..self.cols{
            for r in 0..c{
                let offset = self.cell_to_offset(r, c);
                unsafe {*dst.offset(offset) = z;}
            }
            for r in c..self.rows{
                let offset = self.cell_to_offset(r, c);
                unsafe {*dst.offset(offset) = *src.offset(offset);}
            }
        }
        result
    }

    /// Returns the matrix with permuted rows
    pub fn permuted_rows(&self, permutation : &MatrixU16)->Matrix<T>{
        debug_assert!(permutation.is_col());
        debug_assert_eq!(permutation.num_cells(), self.num_rows());
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let src = self.ptr;
        let dst = result.ptr;
        for c in 0..self.cols{
            let src_start = self.cell_to_offset(0, c);
            let dst_start = result.cell_to_offset(0, c);
            for r in 0..self.num_rows(){
                let src_offset = src_start + permutation[r] as isize;
                let dst_offset = dst_start +  r as isize;
                unsafe {*dst.offset(dst_offset) = *src.offset(src_offset);}
            }
        }
        result
    }

    /// Returns the matrix with permuted columns
    pub fn permuted_cols(&self, permutation : &MatrixU16)->Matrix<T>{
        debug_assert!(permutation.is_col());
        debug_assert_eq!(permutation.num_cells(), self.num_cols());
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let src = self.ptr;
        let dst = result.ptr;
        for c in 0..self.cols{
            debug_assert!((permutation[c] as usize) < self.num_cols());
            let mut src_offset = self.cell_to_offset(0, permutation[c] as usize);
            let mut dst_offset = result.cell_to_offset(0, c);
            for _ in 0..self.num_rows(){
                unsafe {*dst.offset(dst_offset) = *src.offset(src_offset);}
                src_offset += 1;
                dst_offset += 1;
            }
        }
        result
    }

}

impl<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial> Matrix<T> {
    /// Computes power of a matrix
    /// Returns a new matrix
    pub fn pow(&self, exp : usize) -> Matrix<T>{
        if !self.is_square() {
            panic!(SRError::IsNotSquareMatrix.to_string());
        }
        if exp == 0 {
            return Matrix::identity(self.rows, self.cols);
        }
        let mut result = self.clone();
        for _ in 0..(exp -1){
            result = &result * self;
        }
        result
    }




    /// Inner product or dot product of two vectors
    /// Both must be column vectors
    /// Both must have same length
    /// result = a' * b.
    pub fn inner_prod(&self, other : &Matrix<T>) -> T {
        debug_assert!(self.is_col());
        debug_assert!(other.is_col());
        debug_assert!(self.num_cells() == other.num_cells());
        let mut result : T =  Zero::zero();
        let pa = self.ptr;
        let pb = other.ptr;
        for i in 0..self.num_rows(){
            let ii = i as isize;
            let va = unsafe{*pa.offset(ii)};
            let vb = unsafe{*pb.offset(ii)};
            result = result + va * vb;
        }
        result
    }

    /// Outer product of two vectors
    /// Both must be column vectors
    /// Both must have same length
    /// result = a * b'.
    pub fn outer_prod(&self, other : &Matrix<T>) -> Matrix<T> {
        debug_assert!(self.is_col());
        debug_assert!(other.is_col());
        debug_assert!(self.num_cells() == other.num_cells());
        let n = self.num_rows();
        let result : Matrix<T> =  Matrix::new(n, n);
        let pa = self.ptr;
        let pb = other.ptr;
        let pc = result.ptr;
        for r in 0..n{
            for c in 0..n{
                let va = unsafe{*pa.offset(r as isize)};
                let vb = unsafe{*pb.offset(c as isize)};
                let offset = result.cell_to_offset(r, c);
                unsafe{
                    *pc.offset(offset) = va * vb;
                }
            }
        }
        result
    }

}

//+Neg<Output=T>
impl<T:CommutativeGroupAddPartial> Matrix<T> {
    /// Computes the unary minus of a matrix
    pub fn unary_minus(&self)-> Matrix<T> {
        let result : Matrix<T> = Matrix::new(self.cols, self.rows);
        let pa = self.ptr;
        let pc = result.ptr;
        for r in 0..self.rows{
            for c in 0..self.cols{
                let offset = self.cell_to_offset(r, c);
                unsafe {
                    *pc.offset(offset) = -*pa.offset(offset);
                }
            }
        }
        result
    }
}


/// These methods modify the matrix itself
impl<T:MagmaBase> Matrix<T> {

    /// Appends one or more columns at the end of matrix
    pub fn append_columns(&mut self, 
        other : &Matrix<T>
        )-> &mut Matrix<T> {
        let cols = self.cols;
        self.insert_columns(cols, other)
    }

    /// Prepends one or more columns at the beginning of matrix
    pub fn prepend_columns(&mut self,
        other : &Matrix<T> 
        )-> &mut Matrix<T> {
        self.insert_columns(0, other)
    }

    /// Inserts columns at the specified location
    pub fn insert_columns(&mut self,
        index  : usize,
        other : &Matrix<T> 
        )-> &mut Matrix<T> {
        debug_assert_eq!(self.num_rows() , other.num_rows());
        debug_assert!(other.num_cols() > 0);

        // check the capacity.
        let new_cols = self.cols + other.cols;
        if self.xcols < new_cols {
            // We need to reallocate memory.
            let rows = self.rows;
            self.reallocate(rows, new_cols);
        }
        else{
            //println!("Reallocation is not needed.");
        }
        // Now create space for other matrix to be fitted in
        self.create_column_space(index, other.num_cols());
        // Finally copy the column data from the matrix.
        let src_ptr = other.ptr;
        let dst_ptr = self.ptr;
        for i in 0..other.num_cols(){
            let src_offset = other.cell_to_offset(0, i);
            let dst_offset = self.cell_to_offset(0, i + index);
            for j in 0..self.rows{
                let jj = j as isize;
                unsafe {
                    *dst_ptr.offset(dst_offset + jj) = *src_ptr.offset(src_offset + jj);
                }
            }
        }
        // Update the count of columns
        self.cols += other.num_cols();
        self
    }


    /// Appends one or more rows at the bottom of matrix
    pub fn append_rows(&mut self, 
        other : &Matrix<T>
        )-> &mut Matrix<T> {
        let rows = self.rows;
        self.insert_rows(rows, other)
    }

    /// Prepends one or more rows at the top of matrix
    pub fn prepend_rows(&mut self,
        other : &Matrix<T> 
        )-> &mut Matrix<T> {
        self.insert_rows(0, other)
    }

    /// Inserts rows at the specified location
    pub fn insert_rows(&mut self,
        index  : usize,
        other : &Matrix<T> 
        )-> &mut Matrix<T> {
        // Make sure that the dimensions are compatible.
        debug_assert_eq!(self.num_cols() , other.num_cols());
        // Make sure that there is data to be added.
        debug_assert!(other.num_rows() > 0);

        // check the capacity.
        let new_rows = self.rows + other.rows;
        if self.xrows < new_rows {
            // We need to reallocate memory.
            let cols = self.cols;
            self.reallocate(new_rows, cols);
        }
        // Now create space for other matrix to be fitted in
        self.create_row_space(index, other.num_rows());
        // Finally copy the row data from the other matrix.
        let src_ptr = other.ptr;
        let dst_ptr = self.ptr;
        let src_stride = other.stride();
        let dst_stride = self.stride();
        for i in 0..other.num_rows(){
            let src_offset = other.cell_to_offset(i, 0);
            let dst_offset = self.cell_to_offset(i + index, 0);
            for j in 0..self.cols{
                unsafe {
                    *dst_ptr.offset(dst_offset + (j*dst_stride) as isize) = 
                    *src_ptr.offset(src_offset + (j*src_stride) as isize);
                }
            }
        }
        // Update the count of rows
        self.rows += other.num_rows();
        self
    }



    #[doc = "
    Reallocates the underlying buffer.
    Assumes that requested capacity is greater than zero.

    ## Notes

    A simple reallocate cannot be used. When the 
    number of rows changes, the stride also changes.
    Thus, the matrix elements need to be moved around.
    "] 
    fn reallocate(&mut self, rows : usize, cols : usize){
        let new_xrows = rows.next_power_of_two();
        let new_xcols = cols.next_power_of_two();
        let old_capacity = self.xrows * self.xcols;
        let new_capacity = new_xrows *  new_xcols;
        if old_capacity >= new_capacity{
            // Nothing to do.
            return;
        }
        assert!(new_capacity > 0);
        let new_bytes = new_capacity * mem::size_of::<T>();
        println!("Allocating {} bytes", new_bytes);
        let raw = unsafe {
            allocate(new_bytes, mem::min_align_of::<T>())
        };
        let dst_ptr = raw as *mut T;
        let src_ptr = self.ptr;
        // Copy data from source to destination.
        let src_stride = self.xrows;
        let dst_stride = new_xrows;
        for c in 0..self.cols{
            for r in 0..self.rows{
                let src_offset = (c * src_stride + r) as isize;
                let dst_offset = (c * dst_stride + r) as isize;
                unsafe{
                    *dst_ptr.offset(dst_offset) = *src_ptr.offset(src_offset);
                }
            }

        }
        // Drop older data.
        unsafe {
            let old_bytes = old_capacity * mem::size_of::<T>();
            //println!("Allocating {} bytes", old_bytes);
            deallocate(self.ptr as *mut u8,
                   old_bytes,
                   mem::min_align_of::<T>());
        }
        // Keep new data
        self.ptr = dst_ptr;
        self.xrows = new_xrows;
        self.xcols = new_xcols;
    }

    //// Moves column data around and creates space for new columns
    fn create_column_space(&mut self, start: usize, count :usize){
        // The end must not be beyond capacity
        debug_assert!(count + self.cols <= self.xcols);
        debug_assert!(start + count <= self.xcols);
        if start >= self.cols {
            // Nothing to move.
            return;
        }
        let capacity = self.capacity() as isize;
        // count columns starting from start column need to be
        // shifted by count.
        let ptr = self.ptr;
        // Number of columns to shift
        let cols_to_shift =  self.cols - (start + count - 1);
        //println!("start: {:?} count: {:?}, cols to shift {:?}", start, count, cols_to_shift );
        let mut cur_col = self.cols;
        for _ in 0..cols_to_shift{
            cur_col -= 1;
            let dst_col = cur_col + count;
            //println!("src_col: {:?} dst_col: {:?}", cur_col, dst_col);
            let src_offset = self.cell_to_offset(0, cur_col);
            let dst_offset = self.cell_to_offset(0, dst_col);
            debug_assert!(src_offset < capacity);
            debug_assert!(dst_offset < capacity);
            for i in 0..self.rows{
                let ii = i as isize;
                // Some validations
                debug_assert!(src_offset + ii < capacity);
                debug_assert!(dst_offset + ii < capacity);
                unsafe {
                    *ptr.offset(dst_offset + ii) = *ptr.offset(src_offset + ii);
                }
            }
        }
    }

    //// Moves rows around and creates space for new rows
    fn create_row_space(&mut self, start: usize, count :usize){
        // The end must not be beyond capacity
        assert!(start + count <= self.xrows);
        if start >= self.rows {
            // Nothing to move.
            return;
        }
        let capacity = self.capacity() as isize;
        // count rows starting from start row need to be
        // shifted by count.
        let ptr = self.ptr;
        // Number of rows to shift
        let rows_to_shift =  self.rows - (start + count - 1);
        let stride = self.stride();
        let mut cur_row = self.rows;
        for _ in 0..rows_to_shift{
            cur_row -= 1;
            let dst_row = cur_row + count;
            let src_offset = self.cell_to_offset(cur_row, 0);
            let dst_offset = self.cell_to_offset(dst_row, 0);
            debug_assert!(src_offset < capacity);
            debug_assert!(dst_offset < capacity);
            for i in 0..self.cols{
                let ii = (i*stride) as isize;
                // Some validations
                debug_assert!(src_offset + ii < capacity);
                debug_assert!(dst_offset + ii < capacity);
                unsafe {
                    *ptr.offset(dst_offset + ii) = 
                    *ptr.offset(src_offset + ii);
                }
            }
        }
    }

}



/// Implementation of Matrix search operations.
impl<T:MagmaBase+Signed+PartialOrd> Search<T> for Matrix<T> {
}



/// Views of a matrix
impl<T:MagmaBase> Matrix<T> {
    /// Creates a view on the matrix
    pub fn view(&self, start_row : usize, start_col : usize , num_rows: usize, num_cols : usize) -> MatrixView<T> {
        let result : MatrixView<T> = MatrixView::new(self, start_row, start_col, num_rows, num_cols);
        result
    }
}




/// These functions are available only for types which support
/// ordering [at least partial ordering for floating point numbers].
impl<T:CommutativeMonoidAddPartial+PartialOrd> Matrix<T> {


    // Returns the minimum scalar value with location
    pub fn min_scalar(&self) -> (T, usize, usize){
        if self.is_empty(){
            panic!(SRError::EmptyMatrix.to_string());
        }
        let mut v = self.get(0, 0);
        let ps = self.ptr;
        // The location
        let mut rr = 0;
        let mut cc = 0;
        for c in 0..self.cols{
            for r in 0..self.rows{
                let src_offset = self.cell_to_offset(r, c);
                let s = unsafe{*ps.offset(src_offset)};
                if s < v { 
                    v = s;
                    rr = r;
                    cc = c;
                }
            }
        }
        (v, rr, cc)
    }

    // Returns the maximum scalar value with location
    pub fn max_scalar(&self) -> (T, usize, usize){
        if self.is_empty(){
            panic!(SRError::EmptyMatrix.to_string());
        }
        let mut v = self.get(0, 0);
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.ptr;
        for c in 0..self.cols{
            for r in 0..self.rows{
                let src_offset = self.cell_to_offset(r, c);
                let s = unsafe{*ps.offset(src_offset)};
                if s > v { 
                    v = s;
                    rr = r;
                    cc = c;
                }
            }
        }
        (v, rr, cc)
    }
    /// Returns the minimum scalar value
    pub fn min_scalar_value(&self) -> T{
        let (v , _, _) = self.min_scalar();
        v
    }    
    /// Returns the maximum scalar value
    pub fn max_scalar_value(&self) -> T{
        let (v , _, _) = self.max_scalar();
        v
    }    
}

impl<T:MagmaBase+Signed+PartialOrd> Matrix<T> {

    // Returns the absolute minimum scalar value
    pub fn min_abs_scalar(&self) -> (T, usize, usize){
        if self.is_empty(){
            panic!(SRError::EmptyMatrix.to_string());
        }
        let mut v = self.get(0, 0).abs();
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.ptr;
        for c in 0..self.cols{
            for r in 0..self.rows{
                let src_offset = self.cell_to_offset(r, c);
                let s = unsafe{*ps.offset(src_offset)}.abs();
                if s < v { 
                    v = s;
                    rr = r;
                    cc = c;
                };
            }
        }
        (v, rr, cc)
    }

    // Returns the maximum scalar value
    pub fn max_abs_scalar(&self) -> (T, usize, usize){
        if self.is_empty(){
            panic!(SRError::EmptyMatrix.to_string());
        }
        let mut v = self.get(0, 0).abs();
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.ptr;
        for c in 0..self.cols{
            for r in 0..self.rows{
                let src_offset = self.cell_to_offset(r, c);
                let s = unsafe{*ps.offset(src_offset)}.abs();
                if s > v { 
                    v  = s;
                    rr = r;
                    cc = c;
                }
            }
        }
        (v, rr, cc)
    }    

    /// Returns the absolute minimum scalar value
    pub fn min_abs_scalar_value(&self) -> T{
        let (v , _, _) = self.min_abs_scalar();
        v
    }    
    /// Returns the absolute maximum scalar value
    pub fn max_abs_scalar_value(&self) -> T{
        let (v , _, _) = self.max_abs_scalar();
        v
    }    
}

/// These functions are available only for integer matrices
impl<T:MagmaBase + Signed> Matrix<T> {

    /// Returns if an integer matrix is a logical matrix
    //// i.e. all cells are either 0s or 1s.
    pub fn is_logical(&self) -> bool {
        let z : T = Zero::zero();
        let o : T = One::one();
        for c in 0..self.cols{
            for r in 0..self.rows{
                let offset = self.cell_to_offset(r, c);
                unsafe{ 
                    let v = *self.ptr.offset(offset);
                    if v != z && v != o {
                        return false;
                    }
                }
            }
        }
        true
    }
}


/// These functions are available only for floating point matrices
impl<T:FieldPartial+Float> Matrix<T> {
    /// Returns a matrix showing all the cells which are finite
    pub fn is_finite(&self) -> Matrix<u8>{
        let m : Matrix<u8> = Matrix::ones(self.rows, self.cols);
        for c in 0..self.cols{
            for r in 0..self.rows{
                let offset = self.cell_to_offset(r, c);
                unsafe{ 
                    let v = *self.ptr.offset(offset);
                    *m.ptr.offset(offset) = v.is_finite() as u8;
                }
            }
        }
        m
    }

    /// Returns a matrix showing all the cells which are infinite
    pub fn is_infinite(&self) -> Matrix<u8>{
        let m : Matrix<u8> = Matrix::ones(self.rows, self.cols);
        for c in 0..self.cols{
            for r in 0..self.rows{
                let offset = self.cell_to_offset(r, c);
                unsafe{ 
                    let v = *self.ptr.offset(offset);
                    *m.ptr.offset(offset) = v.is_infinite() as u8;
                }
            }
        }
        m
    }
}



impl<T:MagmaBase> Index<usize> for Matrix<T> {
    type Output = T;
    #[inline]
    fn index<'a>(&'a self, index: usize) -> &'a T {
        // The matrix is column major order
        let (r, c) = self.index_to_cell(index);
        let offset = c * self.stride() + r;
        unsafe {
            &*self.ptr.offset(offset as isize)
        }
    }
}



impl<T:MagmaBase> Matrix<T>{
    /// This function is for internal use only.
    #[inline]
    fn as_mut_slice<'a>(&'a mut self) -> &'a mut [T] {
        unsafe {
            mem::transmute(Slice {
                data: self.as_mut_ptr() as *const T,
                len: self.capacity(),
            })
        }
    }
}



/// Implementation of Clone interface
impl <T:MagmaBase> Clone for Matrix<T> {

    /// Creates a clone of the matrix
    fn clone(&self )-> Matrix<T> {
        let m : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = m.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in 0..n{
                let i = i_ as isize;
                *pb.offset(i) = *pa.offset(i);
            }
        }
        m
    }
}

impl <T:MagmaBase> fmt::Debug for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // We need to find out the number of characters needed
        // to show each value.
        // maximum value
        let mut strings : Vec<String> = Vec::with_capacity(self.num_cells());
        let mut n : usize = 0;
        let ptr = self.ptr;
        for r in 0..self.rows{
            for c in 0..self.cols{
                let offset = self.cell_to_offset(r, c);
                let v = unsafe {*ptr.offset(offset)};
                let s = format!("{:?}", v);
                let slen = s.len();
                strings.push(s);
                if slen > n {
                    n = slen;
                }
            }
        }
        //println!("{}, {}, {}", v, s, n);
        try!(write!(f, "["));
        // Here we print row by row
        for r in 0..self.rows {
           try!(write!(f, "\n  "));
            for c in 0..self.cols{
                let offset =  r *self.cols + c;
                let ref s = strings[offset];
                let extra = n + 2 - s.len();
                for _ in 0..extra{
                    try!(write!(f, " "));
                }
                try!(write!(f, "{}", s));
            }
        }
        try!(write!(f, "\n]"));
        Ok(())
    }
}

impl <T:MagmaBase> fmt::Display for Matrix<T> {
    /// Display and Debug versions are same
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}


/// Matrix addition support
impl<'a, 'b, T:CommutativeMonoidAddPartial> ops::Add<&'b Matrix<T>> for &'a Matrix<T> {
    type Output = Matrix<T>;
    fn add(self, rhs: &'b Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in 0..n{
                let i = i_ as isize;
                *pc.offset(i) = *pa.offset(i) + *pb.offset(i);
            }
        }
        result
    }
}


/// Matrix subtraction support
impl<'a, 'b, T:QuasiGroupAddPartial> ops::Sub<&'b Matrix<T>> for &'a Matrix<T>{
    type Output = Matrix<T>;
    fn sub(self, rhs: &'b Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in 0..n{
                let i = i_ as isize;
                *pc.offset(i) = *pa.offset(i) - *pb.offset(i);
            }
        }
        result
    }
}




/// Matrix equality check support
impl<T:MagmaBase> cmp::PartialEq for Matrix<T>{
    fn eq(&self, other: &Matrix<T>) -> bool {
        let pa = self.ptr as *const  T;
        let pb = other.ptr as *const  T;
        for c in 0..self.cols{
            for r in 0..self.rows{
                let offset_a = self.cell_to_offset(r, c);
                let offset_b = other.cell_to_offset(r, c);
                let va = unsafe{*pa.offset(offset_a)};
                let vb = unsafe{*pb.offset(offset_b)};
                if va != vb {
                    return false;
                }
            }
        }
        true
    }
}

// Element wise operations.
impl<T:CommutativeMonoidAddPartial> Matrix<T> {
    /// Adds matrices element by element
    pub fn add_elt(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in 0..n{
                let i = i_ as isize;
                *pc.offset(i) = *pa.offset(i) + *pb.offset(i);
            }
        }
        result
    }
}

impl<T:CommutativeMonoidAddPartial+ops::Sub<Output=T>> Matrix<T> {
    /// Subtracts matrices element by element
    pub fn sub_elt(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in 0..n{
                let i = i_ as isize;
                *pc.offset(i) = *pa.offset(i) - *pb.offset(i);
            }
        }
        result
    }
}
impl<T:CommutativeMonoidMulPartial> Matrix<T> {
    /// Multiplies matrices element by element
    pub fn mul_elt(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in 0..n{
                let i = i_ as isize;
                *pc.offset(i) = *pa.offset(i) * *pb.offset(i);
            }
        }
        result
    }

    /// Computs power of matrix elements
    pub fn pow_elt(&self, n : usize) -> Matrix<T> {
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pc = result.ptr;
        let cap = self.capacity();
        unsafe{
            for i_ in 0..cap{
                let i = i_ as isize;
                let v = *pa.offset(i);
                let mut result : T = One::one();
                for _ in 1..n{
                    result = result * v; 
                }
                *pc.offset(i) = result;
            }
        }
        result
    }
}

impl<T:CommutativeMonoidMulPartial+ops::Div<Output=T>> Matrix<T> {
    /// Divides matrices element by element
    pub fn div_elt(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in 0..n{
                let i = i_ as isize;
                *pc.offset(i) = *pa.offset(i) / *pb.offset(i);
            }
        }
        result
    }

}

impl<T:MagmaBase> Drop for Matrix<T> {
    fn drop(&mut self) {
        if self.num_cells() != 0 {
            unsafe {
                util::memory::dealloc(self.ptr, self.capacity())
            }
        }
    }
}

/******************************************************
 *
 *   Utility functions for debugging of Matrix
 *
 *******************************************************/

impl<T:MagmaBase> Matrix<T> {
    pub fn print_state(&self){
        let capacity = self.capacity();
        let bytes = capacity * mem::size_of::<T>();
        println!("Rows: {}, Cols: {}, XRows : {}, XCols {} , Capacity: {}, Bytes; {}, Buffer: {:p}, End : {:p}", 
            self.rows, self.cols, self.xrows, self.xcols, 
            capacity, bytes, self.ptr, unsafe {
                self.ptr.offset(capacity as isize)
            });
    }
}


/******************************************************
 *
 *   Private implementation of Matrix
 *
 *******************************************************/

impl<T:MagmaBase> Matrix<T> {
    /// Returns a slice into `self`.
    //#[inline]
    pub fn as_slice_<'a>(&'a self) -> &'a [T] {
        unsafe { mem::transmute(Slice { data: self.as_ptr(), len: self.capacity() }) }
    }

}



/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/


#[cfg(test)]
mod test {

    use  super::{Matrix, MatrixI64, MatrixF64};
    use matrix::constructors::*;
    use matrix::traits::*;
    use num::Float;

    #[test]
    fn  test_create_mat0(){
        let m : MatrixI64 = Matrix::new(3, 4);
        assert_eq!(m.num_cells(), 12);
        assert_eq!(m.size(), (3, 4));
        assert_eq!(m.smaller_dim(), 3);
        assert_eq!(m.larger_dim(), 4);
        // NOTE: we cannot use the contents of the
        // matrix as they are uninitialized.
        let m : MatrixI64 = Matrix::zeros(3, 4);
        let v : i64 = m.get(0, 0);
        assert_eq!(v, 0i64);
    }

    #[test]
    fn test_from_slice_rw0(){
        let  m : MatrixI64 = Matrix::from_slice_rw(3, 3, &[
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9
            ]);
        assert_eq!(m.get(0, 0), 1);
        assert_eq!(m.get(0, 1), 2);
        assert_eq!(m.get(0, 2), 3);
        assert_eq!(m.get(1, 0), 4);
        assert_eq!(m.get(1, 1), 5);
        assert_eq!(m.get(1, 2), 6);
        assert_eq!(m.get(2, 0), 7);
        assert_eq!(m.get(2, 1), 8);
    }


    #[test]
    fn test_from_iter_cw0(){
        for _ in 0..100{
            let m : MatrixI64 = Matrix::from_iter_cw(4, 4, 1..20);
            let b: Vec<i64> = (1..17).collect();
            assert!(m.as_slice_() == b.as_slice());
            let m : MatrixI64 = Matrix::from_iter_cw(4, 8, 1..16);
            assert_eq!(m.get(0, 0), 1);
            assert_eq!(m.get(2, 2), 11);
            let mut b: Vec<i64> = (1..16).collect();
            for _ in 0..17{
                b.push(0);
            }
            assert_eq!(m.as_slice_(), b.as_slice());
        }
    }

    #[test]
    fn test_index0(){
        let m : MatrixI64 = Matrix::from_iter_cw(4, 4, (1..20));
        let x = m[4];
        assert_eq!(x, 5);
    }

    #[test]
    fn test_cell_index_mapping(){
        let rows = 20;
        let cols = 10;
        let m : MatrixI64 = Matrix::new(rows, cols);
        assert_eq!(m.index_to_cell(3 + 2*rows), (3, 2));
        assert_eq!(m.cell_to_index(3, 2), 3 + 2*rows);
        assert_eq!(m.cell_to_index(5, 7), 145);
        assert_eq!(m.index_to_cell(145), (5, 7));
    }

    #[test]
    fn test_to_std_vec(){
        let m : MatrixI64 = Matrix::from_iter_cw(4, 3, 0..12);
        let v1 = m.to_std_vec();
        let v2 : Vec<i64> = (0..12).collect();
        assert_eq!(v1, v2);
    }

    #[test]
    fn test_ones(){
        let m : MatrixI64 = Matrix::ones(4, 2);
        let v = vec![1i64, 1, 1, 1, 1, 1, 1, 1];
        assert_eq!(m.to_std_vec(), v);
    }

    #[test]
    fn test_sum(){
        let m : MatrixI64 = Matrix::ones(4, 2);
        let m2 = &m  + &m;
        let v = vec![2i64, 2, 2, 2, 2, 2, 2, 2];
        assert_eq!(m2.to_std_vec(), v);
    }

    #[test]
    #[should_panic]
    fn test_sum_fail(){
        let m1 : MatrixI64 = Matrix::ones(4, 2);
        let m2 : MatrixI64 = Matrix::ones(3, 2);
        &m1 + &m2;
    }

    #[test]
    fn test_sub(){
        let m : MatrixI64 = Matrix::ones(4, 2);
        let m3 = &(&m  + &m) + &m;
        let m2 = &m3 - &m;  
        let v = vec![2i64, 2, 2, 2, 2, 2, 2, 2];
        assert_eq!(m2.to_std_vec(), v);
    }

    #[test]
    fn test_sub_float(){
        let m : MatrixF64 = Matrix::ones(4, 2);
        let m3 = &(&m  + &m) + &m;
        let m2 = &m3 - &m;  
        let v = vec![2f64, 2., 2., 2., 2., 2., 2., 2.];
        assert_eq!(m2.to_std_vec(), v);
    }
    #[test]
    #[should_panic]
    fn test_sub_fail(){
        let m1 : MatrixI64 = Matrix::ones(4, 2);
        let m2 : MatrixI64 = Matrix::ones(3, 2);
        &m1 - &m2;
    }


    #[test]
    fn test_eq(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(2, 2, 0..4);
        let m2 : MatrixI64 = Matrix::from_iter_cw(2, 2, 0..4);
        assert_eq!(m1, m2);
        let v = vec![1.0f64, 2., 3., 4.];
        let m1 : MatrixF64 = Matrix::from_slice_cw(2, 2, v.as_slice());
        let m2 : MatrixF64 = Matrix::from_slice_cw(2, 2, v.as_slice());
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_from_scalar(){
        let  m : MatrixI64  = Matrix::from_scalar(2);
        assert!(m.is_scalar());
        assert_eq!(m.to_std_vec(), vec![2]);
        assert_eq!(m.to_scalar(), 2);
    }
    #[test]
    fn test_rep_mat(){
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 2,  0..4);
        let m2 = m.repeat_matrix(2, 2);
        assert_eq!(m2.num_cells(), 16);
        assert_eq!(m2.num_rows(), 4);
        assert_eq!(m2.num_cols(), 4);
        assert_eq!(m2.to_std_vec(), vec![0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3]);
    }

    #[test]
    fn test_is_vector(){
        let m : MatrixI64 = Matrix::new(3,1);
        assert!(m.is_vector());
        let m : MatrixI64 = Matrix::new(1,4);
        assert!(m.is_vector());
        let m : MatrixI64 = Matrix::new(1,1);
        assert!(!m.is_vector());
        let m : MatrixI64 = Matrix::new(3,3);
        assert!(!m.is_vector());
    }


    #[test]
    fn test_is_empty(){
        let m : MatrixI64 = Matrix::new(3,1);
        assert!(!m.is_empty());
        let m : MatrixI64 = Matrix::new(4, 0);
        assert!(m.is_empty());
        let m : MatrixI64 = Matrix::new(0, 4);
        assert!(m.is_empty());
        let m : MatrixI64 = Matrix::new(0, 0);
        assert!(m.is_empty());
    }

    #[test]
    fn test_is_finite(){
        let v = vec![0f64, Float::nan(), Float::nan(), 
        Float::infinity(), Float::neg_infinity(), 2., 3., 4.];
        let m : MatrixF64 = Matrix::from_slice_cw(2, 4, v.as_slice());
        let f = m.is_finite();
        assert_eq!(f.to_std_vec(), vec![1, 0, 0, 0, 0, 
            1, 1, 1]);
        let f = m.is_infinite();
        assert_eq!(f.to_std_vec(), vec![0, 0, 0, 1, 1, 
            0, 0, 0]);
    }


    #[test]
    fn test_is_logical(){
        let m : MatrixI64 = Matrix::from_iter_cw(4, 4, (0..16).map(|x| x % 2));
        assert!(m.is_logical());
        let m = &m + &m;
        assert!(!m.is_logical());
    }

    #[test]
    fn test_max(){
        // Note the first 4 entries in v form the first column
        // in the matrix.
        let v = vec![4, 5, 11, 12, 
        0, 1, 2, 3, 
        19, 17, 3, 1, 
        7, 8, 5, 6];
        let m : MatrixI64 = Matrix::from_slice_cw(4, 4, v.as_slice());
        let m2 = m.max_row_wise();
        assert!(m2.is_vector());
        assert!(m2.is_col());
        assert_eq!(m2.to_std_vec(), vec![19, 17, 11, 12]);

        let m2 = m.min_row_wise();
        assert!(m2.is_vector());
        assert!(m2.is_col());
        assert_eq!(m2.to_std_vec(), vec![0, 1, 2, 1]);

        let m2 = m.max_col_wise();
        assert!(m2.is_vector());
        assert!(m2.is_row());
        assert_eq!(m2.to_std_vec(), vec![12, 3, 19, 8]);
        


        let m2 = m.min_col_wise();
        assert!(m2.is_vector());
        assert!(m2.is_row());
        assert_eq!(m2.to_std_vec(), vec![4, 0, 1, 5]);

        assert_eq!(m.min_scalar_value(), 0);
        assert_eq!(m.max_scalar_value(), 19);
    }

    #[test]
    fn test_mul_elt(){
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 2,  0..4);
        let m2 = m.mul_elt(&m);
        let m3  : MatrixI64 = Matrix::from_iter_cw(2, 2,  (0..4).map(|x| x*x));
        assert_eq!(m2, m3);
    }


    #[test]
    fn test_div_elt(){
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 2,  (1..20));
        let m2 = &m + &m;
        let m3 = m2.div_elt(&m);
        assert_eq!(m3.to_std_vec(), vec![2,2,2,2]);
    }


    #[test]
    fn test_identity(){
        let m : MatrixI64 = Matrix::identity(3, 2);
        assert_eq!(m.to_std_vec(), vec![1, 0, 0, 0, 1, 0]);
        assert!(m.is_identity());
        let m : MatrixI64 = Matrix::identity(2, 2);
        assert_eq!(m.to_std_vec(), vec![1, 0, 0, 1]);
        assert!(m.is_identity());
        let m : MatrixI64 = Matrix::identity(2, 3);
        assert_eq!(m.to_std_vec(), vec![1, 0, 0, 1, 0, 0]);
        assert!(m.is_identity());
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 2,  0..4);
        assert!(!m.is_identity());
    }

    #[test]
    fn test_is_square(){
        let m : MatrixI64 = Matrix::new(3,4);
        assert!(!m.is_square());
        let m : MatrixI64 = Matrix::new(100,100);
        assert!(m.is_square());
    }

    #[test]
    fn test_pow(){
        let m : MatrixI64  = Matrix::identity(4, 4);
        let m2 = m.copy_mul_scalar(2);
        assert!(m.is_square());
        let m4 = m2.pow(4);
        let m16 = m.copy_mul_scalar(16);
        assert_eq!(m4, m16); 
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 2,  0..4);
        assert!(m.is_square());
        let m3 = m.pow(3);
        assert!(m3.is_square());
        assert_eq!(m3.to_std_vec(), vec![6, 11, 22, 39]);
        assert_eq!(m.pow(1).to_std_vec(), vec![0, 1, 2, 3]);
        assert_eq!(m.pow(2), &m * &m);
        // This is UGLY. I wish I could fix it.
        let x = &(&(&(&(&(&(&(&(&m * &m) * &m) * &m) * &m) * &m) * &m) * &m) * &m) * &m;
        assert_eq!(m.pow(10), x);
    }

    #[test]
    fn test_unary_minus(){
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 2,  0..4);
        let z : MatrixI64 = Matrix::zeros(2,2);
        let m2 = m.unary_minus();
        let m3 = &z - &m;
        assert_eq!(m2, m3);
    }


    #[test]
    fn test_diag_from_vector(){
        // repeat the same test a hundred times
        for _ in 0..100{
            let v  : MatrixI64 = Matrix::from_iter_cw(4, 1, (20..30));
            assert!(v.is_vector());
            let m = Matrix::diag_from_vec(&v);
            assert!(!m.is_empty());
            assert!(!m.is_vector());
            println!("{:?}", m);
            assert!(m.is_diagonal());
            assert_eq!(m.num_cells(), 16);
            let mut m2 : MatrixI64 = Matrix::zeros(4, 4);
            m2.set(0, 0, 20);
            m2.set(1, 1, 21);
            m2.set(2, 2, 22);
            m2.set(3, 3, 23);
            assert_eq!(m, m2);
        }
    }

    #[test]
    fn test_diagonal(){
        let m  : MatrixI64 = Matrix::from_iter_cw(4, 5, (10..30));
        let v = m.diagonal_vector();
        assert!(v.is_vector());
        assert_eq!(v.num_cells(), 4);
        let v2 : MatrixI64 = Matrix::from_slice_cw(4, 1, vec![10, 15, 20, 25].as_slice());
        assert_eq!(v, v2);
    }


    #[test]
    fn test_triangular(){
        let m = matrix_cw_f64(3, 3, &[1., 0., 0., 
            4., 5., 0.,
            6., 2., 3.]);
        println!("m: {:?}", m);
        assert!(m.is_ut());
        assert!(!m.is_lt());
        assert!(m.is_triangular());
        let m  = m.transpose();
        assert!(!m.is_ut());
        assert!(m.is_lt());
        assert!(m.is_triangular());
    }

    #[test]
    fn test_trace(){
        let m = matrix_cw_f64(3, 3, &[1., 0., 0., 
            4., 5., 0.,
            6., 2., 3.]);
        assert_eq!(m.trace(), 9.);
    }

    #[test]
    fn test_extract_triangular(){
        let m = matrix_rw_i64(3,3,&[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9
            ]);
        let mu = matrix_rw_i64(3,3,&[
            1, 2, 3,
            0, 5, 6,
            0, 0, 9
            ]);
        let ml = matrix_rw_i64(3,3,&[
            1, 0, 0,
            4, 5, 0,
            7, 8, 9
            ]);
        assert_eq!(m.ut(), mu);
        assert_eq!(m.lt(), ml);
    }

    #[test]
    fn test_extract_diagonal_matrix(){
        let m = matrix_rw_i64(3,3,&[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9
            ]);
        let md = matrix_rw_i64(3,3,&[
            1, 0, 0,
            0, 5, 0,
            0, 0, 9
            ]);
        assert_eq!(m.diagonal_matrix(), md);
        assert_eq!(m.diagonal_vector(), vector_i64(&[1, 5, 9]));
        // More columns than rows
        let m = matrix_rw_i64(3,4,&[
            1, 2, 3, 11,
            4, 5, 6, 12,
            7, 8, 9, 19
            ]);
        let md = matrix_rw_i64(3,4,&[
            1, 0, 0, 0,
            0, 5, 0, 0,
            0, 0, 9, 0
            ]);
        assert_eq!(m.diagonal_matrix(), md);
        assert_eq!(m.diagonal_vector(), vector_i64(&[1, 5, 9]));
        // More rows than columns
        let m = matrix_rw_i64(4,3,&[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9,
            10, 11, 12
            ]);
        let md = matrix_rw_i64(4,3,&[
            1, 0, 0,
            0, 5, 0,
            0, 0, 9,
            0, 0, 0
            ]);
        assert_eq!(m.diagonal_matrix(), md);
        assert_eq!(m.diagonal_vector(), vector_i64(&[1, 5, 9]));
    }

    #[test]
    fn test_unit_vector(){
        let v : MatrixI64 = Matrix::unit_vector(10, 3);
        let mut m : MatrixI64 = Matrix::zeros(10, 1);
        m.set(3, 0, 1);
        assert_eq!(v, m);
    }

    #[test]
    fn test_row_col_iter(){
        let m  : MatrixI64 = Matrix::from_iter_cw(4, 5, (10..30));
        let r = m.row_iter(0);
        let v : Vec<i64> = r.collect();
        assert_eq!(v, vec![10, 14, 18, 22, 26]);
        let r = m.col_iter(2);
        let v : Vec<i64> = r.collect();
        assert_eq!(v, vec![18, 19, 20, 21]);
        let m  : MatrixI64 = Matrix::from_iter_cw(3, 2, (10..30));
        let r = m.cell_iter();
        let v : Vec<i64> = r.collect();
        assert_eq!(v, vec![10, 11, 12, 13, 14, 15]);
    }

    #[test]
    fn test_add_columns(){
        let mut m1 :  MatrixI64 = Matrix::from_iter_cw(2, 3, (11..100));
        let m2 :  MatrixI64 = Matrix::from_iter_cw(2, 4, (11..100));
        let m3 : MatrixI64  = Matrix::from_iter_cw(2, 1, (17..100));
        let m4 :  MatrixI64 = Matrix::from_iter_cw(2, 5, (9..100));
        let m5 : MatrixI64  = Matrix::from_iter_cw(2, 1, (9..100));
        // matrix size 2x3.
        println!("{}", m1);
        // adding 2x1 to get 2x4.
        m1.append_columns(&m3);
        println!("{}", m1);
        assert_eq!(m1, m2);
        // adding 2x1 to get 2x5.
        m1.prepend_columns(&m5);
        println!("{}", m1);
        assert_eq!(m1, m4);
        let m6 : MatrixI64  = Matrix::from_iter_cw(2, 1, (5..100));
        m1.prepend_columns(&m6);
        println!("{}", m1);
        let m7 : MatrixI64  = Matrix::from_iter_cw(2, 1, (7..100));
        m1.insert_columns(1, &m7);
        println!("{}", m1);
    }


    #[test]
    fn test_add_rows(){
        let mut m1 :  MatrixI64 = Matrix::from_iter_cw(2, 3, (11..100)).transpose();
        let m2 :  MatrixI64 = Matrix::from_iter_cw(2, 4, (11..100)).transpose();
        let m3 : MatrixI64  = Matrix::from_iter_cw(2, 1, (17..100)).transpose();
        let m4 :  MatrixI64 = Matrix::from_iter_cw(2, 5, (9..100)).transpose();
        let m5 : MatrixI64  = Matrix::from_iter_cw(2, 1, (9..100)).transpose();
        println!("m1: {}", m1);
        m1.print_state();
        println!("m2: {}", m2);
        println!("m3: {}", m3);
        println!("m4: {}", m4);
        m4.print_state();
        m1.append_rows(&m3);
        println!("m1 = [m1 ; m3] {}", m1);
        m1.print_state();
        println!("m4: {}", m4);
        m4.print_state();
        //assert_eq!(m1, m2);
        println!("\nm4 {}", m4);
        m4.print_state();
        println!("m5 {}", m5);
        m1.prepend_rows(&m5);
        println!("m1 = [m5 ; m1] {}", m1);
        m1.print_state();
        println!("\nm4 {}", m4);
        m4.print_state();
        assert_eq!(m1, m4);
        let m6 : MatrixI64  = Matrix::from_iter_cw(2, 1, (5..100)).transpose();
        m1.prepend_rows(&m6);
        println!("{}", m1);
        let m7 : MatrixI64  = Matrix::from_iter_cw(2, 1, (7..100)).transpose();
        m1.insert_rows(1, &m7);
        println!("{}", m1);
    }

    #[test]
    fn test_inner_product(){
        let m1 : MatrixI64 = Matrix::from_slice_cw(3, 1, vec![2, 1, 1].as_slice());
        let m2 : MatrixI64 = Matrix::from_slice_cw(3, 1, vec![1, 1, 2].as_slice());
        let result = m1.inner_prod(&m2);
        assert_eq!(result, 5);
    }


    #[test]
    fn test_outer_product(){
        let m1 : MatrixI64 = Matrix::from_slice_cw(3, 1, vec![2, 1, 1].as_slice());
        let m2 : MatrixI64 = Matrix::from_slice_cw(3, 1, vec![1, 1, 2].as_slice());
        let result = m1.outer_prod(&m2);
        let m3 : MatrixI64 = Matrix::from_slice_cw(3, 3, vec![2, 1, 1, 2, 1, 1, 4, 2, 2].as_slice());
        assert_eq!(result, m3);
    }

    #[test]
    fn test_max_abs_scalar_in_row(){
        let m = matrix_rw_i64(3, 3, &[
            2, 93, 9, 
            2, 71, -79, 
            -83, 62, 6]);
        assert_eq!(m.max_abs_scalar_in_row(0, 0, 3), (93, 1));
        assert_eq!(m.max_abs_scalar_in_row(1, 0, 3), (79, 2));
        assert_eq!(m.max_abs_scalar_in_row(2, 0, 3), (83, 0));
        assert_eq!(m.max_abs_scalar_in_row(1, 0, 2), (71, 1));
        assert_eq!(m.max_abs_scalar_in_row(2, 1, 2), (62, 1));
    }

    #[test]
    fn test_max_abs_scalar_in_col(){
        let m = matrix_cw_i64(3, 3, &[
            2, 93, 9, 
            2, 71, -79, 
            -83, 62, 6]);
        assert_eq!(m.max_abs_scalar_in_col(0, 0, 3), (93, 1));
        assert_eq!(m.max_abs_scalar_in_col(1, 0, 3), (79, 2));
        assert_eq!(m.max_abs_scalar_in_col(2, 0, 3), (83, 0));
        assert_eq!(m.max_abs_scalar_in_col(1, 0, 2), (71, 1));
        assert_eq!(m.max_abs_scalar_in_col(2, 1, 2), (62, 1));
    }

    #[test]
    fn test_permute_rows(){
        let m = matrix_rw_i64(4, 3, &[
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9,
            10, 11, 12
            ]);
        let permutation = vector_u16(&[0, 3, 1, 2]);
        let m = m.permuted_rows(&permutation);
        let m2 = matrix_rw_i64(4, 3, &[
            1, 2, 3, 
            10, 11, 12,
            4, 5, 6,
            7, 8, 9,
            ]);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_permute_cols(){
        let m = matrix_rw_i64(4, 3, &[
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9,
            10, 11, 12
            ]);
        let permutation = vector_u16(&[2, 0, 1]);
        let m = m.permuted_cols(&permutation);
        let m2 = matrix_rw_i64(4, 3, &[
            3, 1, 2,
            6, 4, 5,
            9, 7, 8,
            12, 10, 11,
            ]);
        assert_eq!(m, m2);
    }
}

/******************************************************
 *
 *   Bench marks follow.
 *
 *******************************************************/


#[cfg(test)]
mod bench {
    //extern crate test;
    //use self::test::Bencher;

}
