#![doc="Provides the basic Matrix data type
"]

// std imports
use std::mem;
use std::ops;
use std::cmp;
use std::fmt;
use std::vec::Vec;
use num::{Float};
use std::iter::Iterator;
use std::ops::{Index};

// external imports
use num::traits::{Zero, One, Signed};
use num::complex::{Complex32, Complex64};

// local imports
use mod_n; 
use cell_to_loc;

use sralgebra::{MagmaBase, 
    CommutativeMonoidAddPartial, 
    CommutativeMonoidMulPartial,
    QuasiGroupAddPartial,
    CommutativeGroupAddPartial,
    FieldPartial};
use error::SRError;

use iter::*;
use view::MatrixView;
use traits::{Shape, NumberMatrix,
    Strided,
    StridedNumberMatrix,
    StridedFloatMatrix,
    Introspection, 
    MatrixBuffer,
    Search};




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
    /// The pointer to raw data array of the matrix
    //ptr : *mut T,
    /// Underlying vector
    vec : Vec<T>
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
pub type MatrixC32 = Matrix<Complex32>;
/// A matrix of 64-bit complex numbers numbers.
pub type MatrixC64 = Matrix<Complex64>;




/// Static functions for creating  a matrix
impl<T:MagmaBase> Matrix<T> {

#[doc = "Constructs a new matrix of given size (uninitialized).

# Remarks 

The contents of the matrix are not initialized. Hence, 
it doesn't make sense to use this function liberally.
Still the function is internally useful since
different constructor functions need to initialize
the matrix differently.

The underlying vector has 0 length and full capacity here.
The underlying vector must be grown to full capacity before
the matrix can be practically used.
"]
    fn new(rows: usize, cols : usize)-> Matrix<T> {
        debug_assert! (mem::size_of::<T>() != 0);
        let capacity = rows *  cols;

        // // Support for empty  matrices
        // if capacity == 0{
        //     // We do not allocate any memory for the buffer.
        //     // We leave it as a NULL pointer.
        //     return Matrix { rows : rows, 
        //         cols : cols,
        //         ptr : ptr::null_mut()
        //     };
        // }
        let vec = Vec::with_capacity(capacity);
        debug_assert_eq!(vec.len(), 0);
        //let ptr = vec.as_mut_ptr();
        Matrix { rows : rows, 
                cols : cols, 
                //ptr : ptr, 
                vec : vec}
    }

    #[doc = "Constructs a new matrix of given size
    with uninitialized data. 

    This is fine since matrices are supported with types which don't have a destructor.
    "]
    pub fn new_uninitialized(rows: usize, cols : usize)-> Matrix<T> {
        let mut m : Matrix<T> = Matrix::new(rows, cols);
        let capacity  = rows * cols;
        unsafe{m.vec.set_len(capacity)};
        debug_assert_eq!(m.vec.len(), m.capacity());
        m
    }

    #[doc = "Constructs a new matrix of given size initialized
    with a specific value at all its entries.
    "]
    pub fn new_with(rows: usize, cols : usize, value: T)-> Matrix<T> {
        let mut m : Matrix<T> = Matrix::new(rows, cols);
        let capacity  = rows * cols;
        m.vec.extend((0..capacity).map(|_| value));
        debug_assert_eq!(m.vec.len(), m.capacity());
        m
    }
}


/// Static functions for creating  a matrix of numbers
impl<T:CommutativeMonoidAddPartial> Matrix<T> {
    /// Constructs a scalar matrix
    pub fn from_scalar (scalar : T) -> Matrix <T>{
        let mut m : Matrix<T> = Matrix::new(1, 1);
        m.vec.push(scalar);
        debug_assert_eq!(m.vec.len(), m.capacity());
        m
    }

    /// Constructs a matrix of all zeros
    pub fn zeros(rows: usize, cols : usize)-> Matrix<T> {
        let mut m : Matrix<T> = Matrix::new(rows, cols);
        // zero out the memory
        let z : T = Zero::zero();
        let n  =m.stride() * cols;
        m.vec.extend((0..n).map(|_| z));
        debug_assert_eq!(m.vec.len(), m.capacity());
        m
    }

    #[doc = "Constructs a matrix from a slice of data reading
    data in column wise order.
    "]
    pub fn from_slice_cw(rows: usize, cols : usize, values: &[T]) -> Matrix<T>{
        let mut mat : Matrix<T> = Matrix::new(rows, cols);
        // stride of new matrix
        let stride = mat.stride();
        {
            // the empty vector in which data will be filled.
            let ref mut vec = mat.vec;
            // The number of entries we can copy
            let n_values = values.len();
            // zero value for unused locations
            let z : T = Zero::zero();
            // position in source slice
            let mut n = 0;
            for _ in 0..cols{
                // fill the column
                for _ in 0..rows{
                    let v = if n < n_values {
                        values[n]
                    }else{
                        z
                    };
                    vec.push(v);
                    n+=1;
                }
                // fill the rest of column 
                for _ in rows..stride{
                    vec.push(z);
                }
            }
        }
        debug_assert_eq!(mat.vec.len(), mat.capacity());
        // return the filled matrix
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
        let mut mat : Matrix<T> = Matrix::new(rows, cols);
        // stride of new matrix
        let stride = mat.stride();
        {
            // the empty vector in which data will be filled.
            let ref mut vec = mat.vec;
            // get a mutable slice from m
            // The number of entries we can copy
            let n_values = values.len();
            // zero value for unused locations
            let z : T = Zero::zero();
            for c in 0..cols{
                for r in 0..rows{
                    // position in source slice
                    let n = r * cols + c;
                    // pick value from source
                    let v = if n < n_values {
                        values[n]
                    }else{
                        z
                    };
                    // put the value in vector
                    vec.push(v);
                }
                // fill the rest of column 
                for _ in rows..stride{
                    vec.push(z);
                }
            }
        }
        debug_assert_eq!(mat.vec.len(), mat.capacity());
        // return
        mat
    }

    pub fn from_iter_cw< A : Iterator<Item=T>>(rows: usize, cols : usize, mut iter: A) -> Matrix<T>{
        let mut mat : Matrix<T> = Matrix::new(rows, cols);
        let stride = mat.stride();
        {
            let ref mut vec = mat.vec;
            // zero value to be filled in unused locations
            let z : T = Zero::zero();
            let mut completed_columns  = 0;
            'outer: for _ in 0..cols{
                for r in 0..rows{
                    let next_val = iter.next();
                    match next_val{
                        Some(val) => vec.push(val),
                        None => {
                            // Finish this column with zeros
                            for _ in r..rows{
                                vec.push(z);
                            }
                            completed_columns += 1;
                            // fill the rest of column 
                            for _ in rows..stride{
                                vec.push(z);
                            }
                            break 'outer
                        }
                    };
                }
                completed_columns += 1;
            }
            if completed_columns < cols {
                // We  need to fill remaining columns with zeros
                for _ in completed_columns..cols{
                    for _ in 0..stride{
                        vec.push(z);
                    }
                    completed_columns += 1;
                }
            }
        }
        debug_assert_eq!(mat.vec.len(), mat.capacity());
        // return
        mat
    }


    /// Builds a matrix from an iterator reading numbers in a 
    /// row-wise order
    pub fn from_iter_rw< A : Iterator<Item=T>>(rows: usize, cols : usize, 
        iter: A) -> Matrix<T>{
        let mut m : Matrix<T> = Matrix::zeros(rows, cols);
        let nc = m.num_cols();
        let nr = m.num_rows();
        let stride = m.stride();
        {
            let ref mut vec  = m.vec;
            let mut r = 0;
            let mut c = 0;
            for v in iter {
                if c == nc {
                    c = 0;
                    r = r + 1;
                }
                if r == nr {
                    break;
                }
                let dst_offset = cell_to_loc(stride, r, c);
                vec[dst_offset] = v;
                c += 1;
            }
        }
        debug_assert_eq!(m.vec.len(), m.capacity());
        // return
        m
    }

    /// Construct a diagonal matrix from a vector
    pub fn diag_from_vec(v : &Matrix<T>) -> Matrix<T>{
        if !v.is_vector(){
            panic!(SRError::IsNotAVector.to_string());
        }
        let n = v.num_cells();
        let mut m : Matrix<T> = Matrix::zeros(n, n);
        let stride = m.stride();
        {
            let ref src = v.vec;
            let ref mut dst = m.vec;
            // Copy the elements of v in the vector
            for r in 0..n{
                let offset = cell_to_loc(stride, r, r);
                dst[offset] = src[r];
            }
        }
        debug_assert_eq!(m.vec.len(), m.capacity());
        m
    }

}

/// Static functions for creating  a matrix of numbers
impl<T:CommutativeMonoidAddPartial+One> Matrix<T> {

    /// Constructs a matrix of all ones.
    pub fn ones(rows: usize, cols : usize)-> Matrix<T> {
        let mut m : Matrix<T> = Matrix::new(rows, cols);
        let stride = m.stride();
        // fill with ones
        let o : T = One::one();
        let z : T = Zero::zero();
        {
            let ref mut vec = m.vec;
            for _ in 0..cols{
                for _ in 0..rows{
                    vec.push(o);
                }
                for _ in rows..stride{
                    vec.push(z);
                }
            } 
        }
        debug_assert_eq!(m.vec.len(), m.capacity());
        m
    }


    /// Constructs an identity matrix
    pub fn identity(rows: usize, cols : usize) -> Matrix<T> {
        let mut m : Matrix<T> = Matrix::zeros(rows, cols);
        let stride = m.stride();
        {
            // fill with ones
            let ref mut vec = m.vec;
            let one : T = One::one();
            let n = cmp::min(rows, cols);
            for i in 0..n{
                let offset = cell_to_loc(stride, i, i);
                vec[offset] = one;
            }
        }
        debug_assert_eq!(m.vec.len(), m.capacity());
        m
    }

    /// Constructs a unit vector
    /// (1, 0, 0), (0, 1, 0), (0, 0, 1), etc.
    pub fn unit_vector( length : usize, dim : usize) -> Matrix<T> {
        let mut m : Matrix<T> = Matrix::zeros(length, 1);
        m.set(dim, 0, One::one());
        debug_assert_eq!(m.vec.len(), m.capacity());
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


    unsafe fn get_unchecked(&self, r : usize, c : usize) -> T  {
        // These assertions help in checking matrix boundaries
        debug_assert!(r < self.rows);
        debug_assert!(c < self.cols);
        let offset = self.cell_to_location(r, c);
        self.vec[offset]
    }

    fn set(&mut self, r : usize, c : usize, value : T) {
        // These assertions help in checking matrix boundaries
        debug_assert!(r < self.rows);
        debug_assert!(c < self.cols);
        let location = self.cell_to_location(r, c);
        self.vec[location] = value;
    }
}

/// Methods available to number matrices
impl<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial> NumberMatrix<T> for Matrix<T> {

    /// Returns if the matrix is an identity matrix
    fn is_identity(&self) -> bool {
        let o : T = One::one();
        let z  : T = Zero::zero();
        {
            let ref vec = self.vec;
            for c in 0..self.cols{
                for r in 0..self.rows{
                    let offset = self.cell_to_location(r, c);
                    let v = vec[offset];
                    if r == c {
                        if v != o {
                            return false;
                        }
                    }else if v != z {
                        return false;
                    }

                }
            }
        } 
        true
    }

    /// Returns if the matrix is a diagonal matrix
    fn is_diagonal(&self) -> bool {
        let z  : T = Zero::zero();
        let ref vec = self.vec;
        for c in 0..self.cols{
            for r in 0..self.rows{
                if r != c {
                    let offset = self.cell_to_location(r, c);
                    let v = vec[offset];
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
        let ref vec = self.vec;
        for c in 0..self.cols{
            for r in 0..c{
                let offset = self.cell_to_location(r, c);
                let v = vec[offset];
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
        let ref vec = self.vec;
        for c in 0..self.cols{
            for r in c+1..self.rows{
                let offset = self.cell_to_location(r, c);
                let v = vec[offset];
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
                unsafe {
                    if self.get_unchecked(i, j) != self.get_unchecked(j, i) {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Returns the trace of the matrix
    fn trace(&self) -> T{
        let mut result: T = Zero::zero();
        for e in self.diagonal_iter(){
            result = result + e;
        }
        result
    }
}


/// Introspection support
impl<T : MagmaBase> Introspection for Matrix<T> {
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
        self.rows
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
        self.vec.as_ptr()
    }

    /// Returns a mutable unsafe pointer to
    /// the matrix's underlying buffer
    #[inline]
    fn as_mut_ptr(&mut self) -> *mut T{
        self.vec.as_mut_ptr()
    }

    /// Maps a cell index to actual offset in the internal buffer
    #[inline]
    fn cell_to_offset(&self, r : usize,  c: usize)-> isize {
        (c * self.stride() + r) as isize
    } 

    /// Maps a cell index to actual location in the internal vector
    #[inline]
    fn cell_to_location(&self, r : usize,  c: usize)-> usize {
        (c * self.stride() + r)
    } 

}

/// Main methods of a matrix
impl<T:MagmaBase> Matrix<T> {

    /// Returns the capacity of the matrix 
    /// i.e. the number of elements it can hold
    pub fn capacity(&self)-> usize {
        self.rows * self.cols
    }

    /// Reshapes the matrix
    pub fn reshape(&mut self, rows: usize, cols : usize) -> bool {
        let new_capacity = rows * cols;
        if new_capacity != self.capacity(){
            false
        }
        else{
            // Reshape the matrix
            self.rows = rows;
            self.cols = cols;
            true
        }
    }

}



use std::fmt::Debug;

/// Functions to access matrix elements safely and without bounds checking.
impl<T: Debug + Clone + Copy + PartialEq> Matrix<T> {

    /// Returns an iterator over a specific row of matrix
    pub fn row_iter(&self, r : isize) -> RowIterator<T>{
        let r = mod_n(r, self.rows as isize);        
        // Lets find the offset of the begging of the row
        let offset = self.cell_to_offset(r, 0);
        let iter : RowIterator<T> = RowIterator::new(self.cols,
            self.stride(), unsafe {self.vec.as_ptr().offset(offset)} as *const T);
        iter
    }

    /// Returns an iterator over a specific column of the matrix
    pub fn col_iter(&self, c : isize) -> ColIterator<T>{
        let c = mod_n(c, self.cols as isize);        
        // Lets find the offset of the begging of the column
        let offset = self.cell_to_offset(0, c);
        let iter : ColIterator<T> = ColIterator::new(self.rows,
            unsafe {self.vec.as_ptr().offset(offset)} as *const T);
        iter
    }

    /// Returns an iterator over all cells  of the matrix
    pub fn cell_iter(&self) -> CellIterator<T>{
        let iter : CellIterator<T> = CellIterator::new(
            self.rows, self.cols, self.stride(),
            self.vec.as_ptr() as *const T);
        iter
    }

    /// Provide the main diagonal elements
    pub fn diagonal_iter(&self) -> DiagIterator<T>{
        DiagIterator::new(self.smaller_dim(),self.stride(), self.vec.as_ptr())
    }
}

/// Functions to construct new matrices out of a matrix and other conversions
impl<T:CommutativeMonoidAddPartial+One> Matrix<T> {

    // Repeats this matrix in both horizontal and vertical directions 
    pub fn repeat_matrix(&self, num_rows : usize, num_cols : usize) -> Matrix<T> {
        let rows = self.rows * num_rows;
        let cols = self.cols * num_cols;
        let mut result : Matrix<T> = Matrix::zeros(rows, cols);
        let src_stride = self.stride();
        let dst_stride = result.stride();
        {
            let ref mut pd = result.vec;
            let ref ps = self.vec;
            for bc in 0..num_cols{
                let bc_start = bc * self.cols;
                for br in 0..num_rows{
                    let br_start =  br * self.rows;
                    for c in 0..self.cols {
                        for r in 0..self.rows{
                            let src_offset = cell_to_loc(src_stride, r, c);
                            let dst_offset = cell_to_loc(dst_stride, br_start + r, bc_start + c);
                            pd[dst_offset] = ps[src_offset];
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
        let mut result : Matrix<T> = Matrix::new(m, 1);
        {
            let ref mut dst = result.vec;
            for (i, e) in (0..m).zip(self.diagonal_iter()){
                dst[i] = e;
            }
        }
        result        
    }

    /// Extracts the primary diagonal from the matrix as a matrix of same size
    pub fn diagonal_matrix(&self) -> Matrix<T> {
        let m  = cmp::min(self.rows, self.cols);
        let mut result : Matrix<T> = Matrix::zeros(self.rows, self.cols);
        {
            let ref src = self.vec;
            let ref mut dst = result.vec;
            for i in 0..m{
                let offset = self.cell_to_location(i, i);
                dst[offset] = src[offset];
            }
        }
        result        
    }

    /// Returns the upper triangular part of the matrix as a new matrix
    pub fn ut(&self)->Matrix<T>{
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref src = self.vec;
            let ref mut dst = result.vec;
            let z  : T = Zero::zero();
            for c in 0..self.cols{
                for r in 0..(c+1){
                    let offset = self.cell_to_location(r, c);
                    dst[offset] = src[offset];
                }
                for r in (c+1)..self.rows{
                    let offset = self.cell_to_location(r, c);
                    dst[offset] = z;
                }
            }
        }
        result
    }

    /// Returns the lower triangular part of the matrix as a new matrix
    pub fn lt(&self)->Matrix<T>{
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref src = self.vec;
            let ref mut dst = result.vec;
            let z  : T = Zero::zero();
            for c in 0..self.cols{
                for r in 0..c{
                    let offset = self.cell_to_location(r, c);
                    dst[offset] = src[offset];
                }
                for r in c..self.rows{
                    let offset = self.cell_to_location(r, c);
                    dst[offset] = z;
                }
            }
        }
        result
    }

    /// Returns the matrix with permuted rows
    pub fn permuted_rows(&self, permutation : &MatrixU16)->Matrix<T>{
        debug_assert!(permutation.is_col());
        debug_assert_eq!(permutation.num_cells(), self.num_rows());
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let res_stride = result.stride();
        {
            let ref src = self.vec;
            let ref mut dst = result.vec;
            for c in 0..self.cols{
                let src_start = self.cell_to_location(0, c);
                let dst_start = cell_to_loc(res_stride, 0, c);
                for r in 0..self.num_rows(){
                    let src_offset = src_start + (permutation[r] as usize);
                    let dst_offset = dst_start +  r;
                    dst[dst_offset] = src[src_offset];
                }
            }
        }
        result
    }

    /// Returns the matrix with permuted columns
    pub fn permuted_cols(&self, permutation : &MatrixU16)->Matrix<T>{
        debug_assert!(permutation.is_col());
        debug_assert_eq!(permutation.num_cells(), self.num_cols());
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let res_stride = result.stride();
        {
            let ref src = self.vec;
            let ref mut dst = result.vec;
            for c in 0..self.cols{
                debug_assert!((permutation[c] as usize) < self.num_cols());
                let mut src_offset = self.cell_to_location(0, permutation[c] as usize);
                let mut dst_offset = cell_to_loc(res_stride, 0, c);
                for _ in 0..self.num_rows(){
                    dst[dst_offset] = src[src_offset];
                    src_offset += 1;
                    dst_offset += 1;
                }
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
        {
            let ref pa = self.vec;
            let ref pb = other.vec;
            for i in 0..self.num_rows(){
                let va = pa[i];
                let vb = pb[i];
                result = result + va * vb;
            }
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
        let mut result : Matrix<T> =  Matrix::new(n, n);
        let stride = result.stride();
        {
            let ref pa = self.vec;
            let ref pb = other.vec;
            let ref mut pc = result.vec;
            let z = Zero::zero();
            for c in 0..n{
                for r in 0..n{
                    let va = pa[r];
                    let vb = pb[c];
                    pc.push(va * vb);
                }
                for _ in n..stride{
                    pc.push(z);
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
        let mut result : Matrix<T> = Matrix::new(self.cols, self.rows);
        {
            let ref pa = self.vec;
            let ref mut pc = result.vec;
            for v in pa.iter(){
                let x = -(*v);
                pc.push(x);
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

        // Now create space for other matrix to be fitted in
        {
            // Original matrix X = [A B]
            // Final matrix Y = [A C B]
            // a columns in A , c columns in C and b columns in B.
            // total original columns = a + b
            // total final columns = a + b + c
            // if b = 0, then all new columns are pushed in the end.
            // if b = 1, then 1st column of C goes to last column of X, 
            //    remaining c-1 columns go to the new location and
            //    the last column of X goes to last column of Y.
            // order of operation: 
            // append  tail c-1 columns of C, then the last column of X, then first column of C
            // generalized: append c-b columns of C, then b columns of X then b columns of C
            // identity columns to move.
            // move those columns 
            // This algorithm is complicated. TBD later.
        }
        // check the capacity.
        let old_cols = self.cols;
        let count = other.cols;
        let stride = self.stride();
        let old_capacity = self.capacity();
        let rows = self.rows;
        self.reallocate(rows, old_cols + count);
        let new_capacity = self.capacity();
        // just get the first element of second matrix
        let z = other.vec[0];
        let extra = new_capacity - old_capacity;
        // copy it around
        self.vec.extend((0..extra).map(|_| z));
        debug_assert_eq!(self.vec.len(), self.capacity());
        // now move things around
        self.create_column_space(index, other.num_cols());
        debug_assert_eq!(self.vec.len(), self.capacity());
        // Finally copy the column data from the matrix.
        {
            let ref src = other.vec;
            let ref mut dst = self.vec;
            for i in 0..other.num_cols(){
                let src_offset = other.cell_to_location(0, i);
                let dst_offset = cell_to_loc(stride, 0, i + index);
                for j in 0..self.rows{
                    dst[dst_offset + j] = src[src_offset + j];
                }
            }
        }
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
        // We need to reallocate memory.
        let cols = self.cols;
        self.reallocate(new_rows, cols);
        // Now create space for other matrix to be fitted in
        self.create_row_space(index, other.num_rows());
        // Finally copy the row data from the other matrix.
        let src_stride = other.stride();
        let dst_stride = self.stride();
        {
            let ref src = other.vec;
            let ref mut dst = self.vec;
            for i in 0..other.num_rows(){
                let src_offset = other.cell_to_location(i, 0);
                let dst_offset = cell_to_loc(dst_stride, i + index, 0);
                for j in 0..self.cols{
                    dst[dst_offset + (j*dst_stride) ] = src[src_offset + (j*src_stride)];
                }
            }
        }
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
        let old_capacity = self.rows * self.cols;
        let new_capacity = rows *  cols;
        if old_capacity >= new_capacity{
            // Nothing to do.
            return;
        }
        assert!(new_capacity > 0);
        let len = self.vec.len();
        debug_assert_eq!(len, self.capacity());
        let extra = new_capacity - len;
        self.vec.reserve_exact(extra);
        self.rows = rows;
        self.cols = cols;
        debug_assert_eq!(self.vec.capacity(), new_capacity);
    }

    //// Moves column data around and creates space for new columns
    fn create_column_space(&mut self, start: usize, count :usize){
        // We need to reallocate memory.
        // The end must not be beyond capacity
        let new_cols = self.cols;
        let old_cols = new_cols - count;
        if start >= new_cols {
            // Nothing to move.
            return;
        }
        let capacity = self.capacity();
        // count columns starting from start column need to be
        // shifted by count.
        // Number of columns to shift
        let cols_to_shift =  old_cols - (start + count - 1);
        let stride = self.stride();
        {
            let ref mut vec = self.vec;
            //println!("start: {:?} count: {:?}, cols to shift {:?}", start, count, cols_to_shift );
            let mut cur_col = old_cols;
            for _ in 0..cols_to_shift{
                cur_col -= 1;
                let dst_col = cur_col + count;
                //println!("src_col: {:?} dst_col: {:?}", cur_col, dst_col);
                let src_offset = cell_to_loc(stride, 0, cur_col);
                let dst_offset = cell_to_loc(stride, 0, dst_col);
                debug_assert!(src_offset < capacity);
                debug_assert!(dst_offset < capacity);
                for i in 0..stride{
                    // Some validations
                    debug_assert!(dst_offset + i < capacity);
                    let src_value = vec[src_offset + i];
                    vec[dst_offset + i] = src_value;
                }
            }
        }
        debug_assert_eq!(self.vec.capacity(), self.capacity());
        debug_assert_eq!(self.vec.len(), self.capacity());
    }

    //// Moves rows around and creates space for new rows
    fn create_row_space(&mut self, start: usize, count :usize){
        // The end must not be beyond capacity
        let new_rows = self.rows;
        let old_rows = new_rows - count;
        if start >= old_rows {
            // Nothing to move.
            return;
        }
        let capacity = self.capacity();
        // Number of rows to shift
        let rows_to_shift =  old_rows - (start + count - 1);
        let stride = self.stride();
        // count rows starting from start row need to be
        // shifted by count.
        {
            let ref mut vec = self.vec;
            let mut cur_row = old_rows;
            for _ in 0..rows_to_shift{
                cur_row -= 1;
                let dst_row = cur_row + count;
                let src_offset = cell_to_loc(stride, cur_row, 0);
                let dst_offset = cell_to_loc(stride, dst_row, 0);
                debug_assert!(src_offset < capacity);
                debug_assert!(dst_offset < capacity);
                for i in 0..self.cols{
                    let ii = i*stride;
                    // Some validations
                    debug_assert!(src_offset + ii < capacity);
                    debug_assert!(dst_offset + ii < capacity);
                    vec[dst_offset + ii] = vec[src_offset + ii];
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
        let mut v = unsafe {self.get_unchecked(0, 0)};
        let ps = self.vec.as_ptr();
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
        let mut v = unsafe {self.get_unchecked(0, 0)};
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.vec.as_ptr();
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
        let mut v = unsafe {self.get_unchecked(0, 0)}.abs();
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.vec.as_ptr();
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
        let mut v = unsafe {self.get_unchecked(0, 0)}.abs();
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.vec.as_ptr();
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
        let ptr = self.vec.as_ptr();
        for c in 0..self.cols{
            for r in 0..self.rows{
                let offset = self.cell_to_offset(r, c);
                unsafe{ 
                    let v = *ptr.offset(offset);
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
        let mut m : Matrix<u8> = Matrix::ones(self.rows, self.cols);
        {
            let ptr = self.vec.as_ptr();
            let mptr = m.vec.as_mut_ptr();
            for c in 0..self.cols{
                for r in 0..self.rows{
                    let offset = self.cell_to_offset(r, c);
                    unsafe{ 
                        let v = *ptr.offset(offset);
                        *mptr.offset(offset) = v.is_finite() as u8;
                    }
                }
            }
        }
        m
    }

    /// Returns a matrix showing all the cells which are infinite
    pub fn is_infinite(&self) -> Matrix<u8>{
        let mut m : Matrix<u8> = Matrix::ones(self.rows, self.cols);
        {
            let ptr = self.vec.as_ptr();
            let mptr = m.vec.as_mut_ptr();
            for c in 0..self.cols{
                for r in 0..self.rows{
                    let offset = self.cell_to_offset(r, c);
                    unsafe{ 
                        let v = *ptr.offset(offset);
                        *mptr.offset(offset) = v.is_infinite() as u8;
                    }
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
        &self.vec[offset]
    }
}



/// Implementation of Clone interface
impl <T:MagmaBase> Clone for Matrix<T> {

    /// Creates a clone of the matrix
    fn clone(&self )-> Matrix<T> {
        let mut m : Matrix<T> = Matrix::new(self.rows, self.cols);
        m.vec.extend(self.vec.as_slice());
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
        let ref vec = self.vec;
        for r in 0..self.rows{
            for c in 0..self.cols{
                let offset = self.cell_to_location(r, c);
                let s = format!("{:?}", vec[offset]);
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
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref a_vec = self.vec;
            let ref b_vec = rhs.vec;
            let ref mut result_vec = result.vec;
            let n = self.capacity();
            for i in 0..n{
                result_vec.push(a_vec[i] + b_vec[i]);
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
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref a_vec = self.vec;
            let ref b_vec = rhs.vec;
            let ref mut result_vec = result.vec;
            let n = self.capacity();
            for i in 0..n{
                result_vec.push(a_vec[i] - b_vec[i]);
            }
        }
        result
    }
}




/// Matrix equality check support
impl<T:MagmaBase> cmp::PartialEq for Matrix<T>{
    fn eq(&self, other: &Matrix<T>) -> bool {
        let pa = self.vec.as_ptr();
        let pb = other.vec.as_ptr();
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
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref a_vec = self.vec;
            let ref b_vec = rhs.vec;
            let ref mut result_vec = result.vec;
            let n = self.capacity();
            for i in 0..n{
                result_vec.push(a_vec[i] + b_vec[i]);
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
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref a_vec = self.vec;
            let ref b_vec = rhs.vec;
            let ref mut result_vec = result.vec;
            let n = self.capacity();
            for i in 0..n{
                result_vec.push(a_vec[i] - b_vec[i]);
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
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref a_vec = self.vec;
            let ref b_vec = rhs.vec;
            let ref mut result_vec = result.vec;
            let n = self.capacity();
            for i in 0..n{
                result_vec.push(a_vec[i] * b_vec[i]);
            }
        }
        result
    }

    /// Computs power of matrix elements
    pub fn pow_elt(&self, n : usize) -> Matrix<T> {
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref pa = self.vec;
            let ref mut pc = result.vec;
            let cap = self.capacity();
            for i in 0..cap{
                let v = pa[i];
                let mut result : T = One::one();
                for _ in 1..n{
                    result = result * v; 
                }
                pc[i] = result;
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
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        {
            let ref a_vec = self.vec;
            let ref b_vec = rhs.vec;
            let ref mut result_vec = result.vec;
            let n = self.capacity();
            for i in 0..n{
                result_vec.push(a_vec[i] / b_vec[i]);
            }
        }
        result
    }

}

impl<T:MagmaBase> Drop for Matrix<T> {
    fn drop(&mut self) {
        self.vec.truncate(0);
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
        let ptr = self.vec.as_ptr();
        println!("Rows: {}, Cols: {}, Capacity: {}, Bytes; {}, Buffer: {:p}, End : {:p}", 
            self.rows, self.cols,
            capacity, bytes, ptr, unsafe {
                ptr.offset(capacity as isize)
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
        self.vec.as_slice()
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
    use constructors::*;
    use traits::*;
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
        let v : i64 = m.get(0, 0).unwrap();
        assert_eq!(v, 0i64);
    }

    #[test]
    fn test_from_slice_rw0(){
        let  m : MatrixI64 = Matrix::from_slice_rw(3, 3, &[
            1, 2, 3, 
            4, 5, 6,
            7, 8, 9
            ]);
        assert_eq!(m.get(0, 0).unwrap(), 1);
        assert_eq!(m.get(0, 1).unwrap(), 2);
        assert_eq!(m.get(0, 2).unwrap(), 3);
        assert_eq!(m.get(1, 0).unwrap(), 4);
        assert_eq!(m.get(1, 1).unwrap(), 5);
        assert_eq!(m.get(1, 2).unwrap(), 6);
        assert_eq!(m.get(2, 0).unwrap(), 7);
        assert_eq!(m.get(2, 1).unwrap(), 8);
    }


    #[test]
    fn test_from_iter_cw0(){
        for _ in 0..100{
            let m : MatrixI64 = Matrix::from_iter_cw(4, 4, 1..20);
            let b: Vec<i64> = (1..17).collect();
            assert!(m.as_slice_() == b.as_slice());
            let m : MatrixI64 = Matrix::from_iter_cw(4, 8, 1..16);
            assert_eq!(m.get(0, 0).unwrap(), 1);
            assert_eq!(m.get(2, 2).unwrap(), 11);
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
    fn test_reshape(){
        let v = vec![1, 2, 3, 4, 
        5, 6, 7, 8,
        9, 10, 11, 12];
        let mut m1 : MatrixI64 = Matrix::from_slice_cw(8, 1, v.as_slice());
        m1.reshape(2, 4);
        let m2 : MatrixI64 = Matrix::from_slice_cw(2, 4, v.as_slice());
        assert_eq!(m1, m2);
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
        println!("m: {} v: {}", m, v);
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
        println!("{}", m.trace());
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
