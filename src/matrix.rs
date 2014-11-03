// library imports
use std::mem;
use std::ptr;
use std::ops;
use std::cmp;
use std::fmt;
use std::num;
use std::num::{One, Zero};
use std::iter::Iterator;
use std::rt::heap::{allocate, deallocate};
use std::raw::Slice as RawSlice;


// srmat imports

use discrete::{mod_n};
use matelt::{MatrixElt};
use materr::*;
use matiter::*;
use matview::MatrixView;


// The following is needed for destroying matrix.




#[doc = "
Represents a matrix of numbers.


The numbers are stored in column-major order.
This is the standard in Fortran and MATLAB.


Currently, the stride is same as the number
of rows. May be changed later on.

For a row vector, the stride  = 1.

"]
pub struct Matrix<T:MatrixElt> {
    /// Number of rows in the matrix
    rows : uint,
    /// Number of columns in the matrix
    cols : uint, 
    /// Number of allocated rows 
    xrows : uint, 
    /// Number of allocated columns
    xcols : uint,
    /// The pointer to raw data array of the matrix
    ptr : *mut T
}

/// A matrix of 64-bit signed integers
pub type MatrixI64 = Matrix<i64>;
/// A matrix of 64-bit floating point numbers.
pub type MatrixF64 = Matrix<f64>;




/// Static functions for creating  a matrix
impl<T:MatrixElt> Matrix<T> {

    /// Constructs a scalar matrix
    pub fn from_scalar (scalar : T) -> Matrix <T>{
        let m : Matrix<T> = Matrix::new(1, 1);
        unsafe {*m.ptr = scalar;}
        m
    }

    pub fn new(rows: uint, cols : uint)-> Matrix<T> {
        assert! (mem::size_of::<T>() != 0);
        let xrows = num::next_power_of_two(rows);
        let xcols = num::next_power_of_two(cols);
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

    /// Constructs a matrix of all zeros
    pub fn zeros(rows: uint, cols : uint)-> Matrix<T> {
        let m : Matrix<T> = Matrix::new(rows, cols);
        // zero out the memory
        unsafe { ptr::zero_memory(m.ptr, m.capacity())};
        m
    }


    /// Constructs a matrix of all ones.
    pub fn ones(rows: uint, cols : uint)-> Matrix<T> {
        let m : Matrix<T> = Matrix::new(rows, cols);
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


    /// Constructs an identity matrix
    pub fn identity(rows: uint, cols : uint) -> Matrix<T> {
        let m : Matrix<T> = Matrix::zeros(rows, cols);
        // fill with ones
        let ptr = m.ptr;
        let one : T = One::one();
        let n = cmp::min(rows, cols);
        for i in range(0, n){
            let offset = m.cell_to_offset(i, i);
            unsafe{
                *ptr.offset(offset) = one;
            }
        }
        m
    }


    pub fn from_slice(rows: uint, cols : uint, values: &[T]) -> Matrix<T>{
        let mut mat : Matrix<T> = Matrix::new(rows, cols);
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

    pub fn from_iter< A : Iterator<T>>(rows: uint, cols : uint, mut iter: A) -> Matrix<T>{
        let mut mat : Matrix<T> = Matrix::new(rows, cols);
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

    /// Construct a diagonal matrix from a vector
    pub fn diag_from_vec(v : &Matrix<T>) -> Matrix<T>{
        if !v.is_vector(){
            fail!(NotAVector.to_string());
        }
        let n = v.num_cells();
        let m : Matrix<T> = Matrix::new(n, n);
        let src = v.ptr;
        let dst = m.ptr;
        // Copy the elements of v in the vector
        for r in range(0, n){
            let offset = m.cell_to_offset(r, r);
            unsafe{
                *dst.offset(offset) =  *src.offset(r as int);
            }
        }
        m
    }

    /// Constructs a unit vector
    /// (1, 0, 0), (0, 1, 0), (0, 0, 1), etc.
    pub fn unit_vector( length : uint, dim : uint) -> Matrix<T> {
        let mut m : Matrix<T> = Matrix::zeros(length, 1);
        m.set(dim, 0, One::one());
        m
    }
}

/// Main methods of a matrix
impl<T:MatrixElt> Matrix<T> {
    //// Returns the number of rows in the matrix
    pub fn num_rows(&self) -> uint {
        self.rows
    }

    /// Returns the number of columns in the matrix
    pub fn num_cols(&self) -> uint {
        self.cols
    }

    /// Returns the size of matrix in an (r, c) tuple
    pub fn size (&self)-> (uint, uint){
        (self.rows, self.cols)
    }
    /// Returns the number of cells in matrix
    pub fn num_cells(&self)->uint {
        self.rows * self.cols
    }
    /// Indicates if the matrix is a row vector
    pub fn is_row(&self) -> bool {
        self.rows == 1
    }

    /// Indicates if the matrix is a column vector
    pub fn is_col(&self) -> bool {
        self.cols == 1
    }

    /// Indicates if the matrix is a scalar actually
    pub fn is_scalar(&self) -> bool {
        self.num_cells() == 1
    }

    /// Indicates if the matrix is a vector
    pub fn is_vector(&self) -> bool {
        (self.rows == 1) ^ (self.cols == 1)
    } 

    /// Indicates if the matrix is empty
    pub fn is_empty(&self) -> bool {
        self.rows * self.cols == 0
    }

    /// Indicates if the matrix is square
    pub fn is_square(&self) -> bool {
        self.rows == self.cols
    }

    /// Returns the number of actual memory elements 
    /// per column stored in the memory
    pub fn stride (&self)->uint {
        self.xrows
    }

    /// Returns the capacity of the matrix 
    /// i.e. the number of elements it can hold
    pub fn capacity(&self)-> uint {
        self.xrows * self.xcols
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

    #[inline]
    pub fn set(&mut self, r : uint, c : uint, value : T) {
        unsafe {
            *self.ptr.offset(self.cell_to_offset(r, c) as int) = value;
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

    /// Returns if the matrix is an identity matrix
    pub fn is_identity(&self) -> bool {
        let o : T = One::one();
        let z  : T = Zero::zero();
        let ptr = self.ptr;
        for c in range(0, self.cols){
            for r in range (0, self.rows){
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
    pub fn is_diagonal(&self) -> bool {
        let z  : T = Zero::zero();
        let ptr = self.ptr;
        for c in range(0, self.cols){
            for r in range (0, self.rows){
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
}


/// Functions to construct new matrices out of a matrix and other conversions
impl<T:MatrixElt> Matrix<T> {
    /// Converts the matrix to vector from standard library
    pub fn to_std_vec(&self) -> Vec<T> {
        let mut vec: Vec<T> = Vec::with_capacity(self.num_cells());
        // We iterate over elements in matrix and push in the vector
        let ptr = self.ptr;
        for c in range(0, self.cols){
            for r in range (0, self.rows){
                let offset = self.cell_to_offset(r, c);
                vec.push(unsafe{*ptr.offset(offset)});
            }
        } 
        vec
    }

    /// Converts the matrix to a scalar 
    pub fn to_scalar(&self) -> T {
        if !self.is_scalar() {
            fail! (DimensionsMismatch.to_string());
        }
        self.get(0, 0)
    }

    /// Returns the r'th row vector
    pub fn row(&self, r : int) -> Matrix<T> {
        // Lets ensure that the row value is mapped to
        // a value in the range [0, rows - 1]
        let r = mod_n(r, self.rows as int);        
        let result : Matrix<T> = Matrix::new(1, self.cols);
        let pd = result.ptr;
        let ps = self.ptr;
        for c in range(0, self.cols){
            let src_offset = self.cell_to_offset(r, c);
            let dst_offset = result.cell_to_offset(0, c);
            unsafe{
                *pd.offset(dst_offset) = *ps.offset(src_offset);
            }
        }
        result
    }

    /// Returns the c'th column vector
    pub fn col(&self, c : int) -> Matrix<T> {
        // Lets ensure that the col value is mapped to
        // a value in the range [0, cols - 1]
        let c = mod_n(c, self.cols as int);        
        let result : Matrix<T> = Matrix::new(self.rows, 1);
        let pd = result.ptr;
        let ps = self.ptr;
        for r in range(0, self.rows){
            let src_offset = self.cell_to_offset(r, c);
            let dst_offset = result.cell_to_offset(r, 0);
            unsafe{
                *pd.offset(dst_offset) = *ps.offset(src_offset);
            }
        }
        result
    }

    /// Returns an iterator over a specific row of matrix
    pub fn row_iter(&self, r : int) -> RowIterator<T>{
        let r = mod_n(r, self.rows as int);        
        // Lets find the offset of the begging of the row
        let offset = self.cell_to_offset(r, 0);
        let iter : RowIterator<T> = RowIterator::new(self.cols,
            self.stride(), unsafe {self.ptr.offset(offset)} as *const T);
        iter
    }

    /// Returns an iterator over a specific column of the matrix
    pub fn col_iter(&self, c : int) -> ColIterator<T>{
        let c = mod_n(c, self.cols as int);        
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


    /// Extract a submatrix from the matrix
    /// rows can easily repeat if the number of requested rows is higher than actual rows
    /// cols can easily repeat if the number of requested cols is higher than actual cols
    pub fn sub_matrix(&self, start_row : int, start_col : int , num_rows: uint, num_cols : uint) -> Matrix<T> {
        let r = mod_n(start_row, self.rows as int);        
        let c = mod_n(start_col, self.cols as int);
        let result : Matrix<T> = Matrix::new(num_rows, num_cols);
        let pd = result.ptr;
        let ps = self.ptr;
        let mut dc = 0;
        for c in range(c, c + num_cols).map(|x | x % self.cols) {
            let mut dr = 0;
            for r in range(r , r + num_rows).map(|x|  x % self.rows) {
                let src_offset = self.cell_to_offset(r, c);
                let dst_offset = result.cell_to_offset(dr, dc);
                unsafe{
                    *pd.offset(dst_offset) = *ps.offset(src_offset);
                }
                dr += 1;
            }
            dc += 1;
        }
        result
    }

    // Repeats this matrix in both horizontal and vertical directions 
    pub fn repeat_matrix(&self, num_rows : uint, num_cols : uint) -> Matrix<T> {
        let rows = self.rows * num_rows;
        let cols = self.cols * num_cols;
        let result : Matrix<T> = Matrix::new(rows, cols);
        let pd = result.ptr;
        let ps = self.ptr;
        for bc in range(0, num_cols){
            let bc_start = bc * self.cols;
            for br in range(0, num_rows){
                let br_start =  br * self.rows;
                for c in range(0, self.cols) {
                    for r in range(0, self.rows){
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
    pub fn diagonal(&self) -> Matrix<T> {
        let m  = cmp::min(self.rows, self.cols);
        let result : Matrix<T> = Matrix::new(m, 1);
        let src = self.ptr;
        let dst = result.ptr;
        for i in range(0, m){
            let offset = self.cell_to_offset(i, i);
            unsafe{
                *dst.offset(i as int) = *src.offset(offset);
            } 
        }
        result        
    }

    /// Add the matrix by a scalar
    /// Returns a new matrix
    pub fn add_scalar(&self, rhs: T) -> Matrix<T> {
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pc = result.ptr;
        for r in range(0, self.rows){
            for c in range(0, self.cols){
                let offset = self.cell_to_offset(r, c);
                unsafe {
                    *pc.offset(offset) = *pa.offset(offset) + rhs;
                }
            }
        }
        result
    }


    /// Multiply the matrix by a scalar
    /// Returns a new matrix
    pub fn mul_scalar(&self, rhs: T) -> Matrix<T> {
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pc = result.ptr;
        for r in range(0, self.rows){
            for c in range(0, self.cols){
                let offset = self.cell_to_offset(r, c);
                unsafe {
                    *pc.offset(offset) = *pa.offset(offset) * rhs;
                }
            }
        }
        result
    }

    /// Divide the matrix by a scalar
    /// Returns a new matrix
    pub fn div_scalar(&self, rhs: T) -> Matrix<T> {
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pc = result.ptr;
        for r in range(0, self.rows){
            for c in range(0, self.cols){
                let offset = self.cell_to_offset(r, c);
                unsafe {
                    *pc.offset(offset) = *pa.offset(offset) / rhs;
                }
            }
        }
        result
    }

    /// Computes power of a matrix
    /// Returns a new matrix
    pub fn pow(&self, exp : uint) -> Matrix<T>{
        if !self.is_square() {
            fail!(NonSquareMatrix.to_string());
        }
        if exp == 0 {
            return Matrix::identity(self.rows, self.cols);
        }
        let mut result = self.clone();
        for _ in range(0, exp -1){
            result = result * *self;
        }
        result
    }

    /// Computes the transpose of a matrix.
    /// This doesn't involve complex conjugation.
    /// Returns a new matrix
    pub fn transpose(&self) -> Matrix <T>{
        let result : Matrix<T> = Matrix::new(self.cols, self.rows);
        let pa = self.ptr;
        let pc = result.ptr;
        for r in range(0, self.rows){
            for c in range(0, self.cols){
                let src_offset = self.cell_to_offset(r, c);
                let dst_offset = result.cell_to_offset(c, r);
                unsafe {
                    *pc.offset(dst_offset) = *pa.offset(src_offset);
                }
            }
        }
        result
    }

    /// Computes the unary minus of a matrix
    pub fn unary_minus(&self)-> Matrix<T> {
        let result : Matrix<T> = Matrix::new(self.cols, self.rows);
        let pa = self.ptr;
        let pc = result.ptr;
        for r in range(0, self.rows){
            for c in range(0, self.cols){
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
impl<T:MatrixElt> Matrix<T> {

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
        index  : uint,
        other : &Matrix<T> 
        )-> &mut Matrix<T> {
        assert_eq!(self.num_rows() , other.num_rows());
        assert!(other.num_cols() > 0);

        // check the capacity.
        let new_cols = self.cols + other.cols;
        if self.xcols < new_cols {
            // We need to reallocate memory.
            let rows = self.rows;
            self.reallocate(rows, new_cols);
        }
        // Now create space for other matrix to be fitted in
        self.create_column_space(index, other.num_cols());
        // Finally copy the column data from the matrix.
        let src_ptr = other.ptr;
        let dst_ptr = self.ptr;
        for i in range(0, other.num_cols()){
            let src_offset = other.cell_to_offset(0, i);
            let dst_offset = self.cell_to_offset(0, i + index);
            for j in range(0, self.rows){
                let jj = j as int;
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
        index  : uint,
        other : &Matrix<T> 
        )-> &mut Matrix<T> {
        // Make sure that the dimensions are compatible.
        assert_eq!(self.num_cols() , other.num_cols());
        // Make sure that there is data to be added.
        assert!(other.num_rows() > 0);

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
        for i in range(0, other.num_rows()){
            let src_offset = other.cell_to_offset(i, 0);
            let dst_offset = self.cell_to_offset(i + index, 0);
            for j in range(0, self.cols){
                unsafe {
                    *dst_ptr.offset(dst_offset + (j*dst_stride) as int) = 
                    *src_ptr.offset(src_offset + (j*src_stride) as int);
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
    fn reallocate(&mut self, rows : uint, cols : uint){
        let new_xrows = num::next_power_of_two(rows);
        let new_xcols = num::next_power_of_two(cols);
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
        for c in range(0, self.cols){
            for r in range(0, self.rows){
                let src_offset = (c * src_stride + r) as int;
                let dst_offset = (c * dst_stride + r) as int;
                unsafe{
                    *dst_ptr.offset(dst_offset) = *src_ptr.offset(src_offset);
                }
            }

        }
        // Drop older data.
        unsafe {
            let old_bytes = old_capacity * mem::size_of::<T>();
            println!("Allocating {} bytes", old_bytes);
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
    fn create_column_space(&mut self, start: uint, count :uint){
        // The end must not be beyond capacity
        assert!(start + count <= self.xcols);
        if start >= self.cols {
            // Nothing to move.
            return;
        }
        let capacity = self.capacity() as int;
        // count columns starting from start column need to be
        // shifted by count.
        let mut cur_col = self.cols - 1;
        let ptr = self.ptr;
        // Number of columns to shift
        let cols_to_shift =  self.cols - (start + count - 1);
        for _ in range(0, cols_to_shift){
            let dst_col = cur_col + count;
            let src_offset = self.cell_to_offset(0, cur_col);
            let dst_offset = self.cell_to_offset(0, dst_col);
            assert!(src_offset < capacity);
            assert!(dst_offset < capacity);
            for i in range(0, self.rows){
                let ii = i as int;
                // Some validations
                assert!(src_offset + ii < capacity);
                assert!(dst_offset + ii < capacity);
                unsafe {
                    *ptr.offset(dst_offset + ii) = *ptr.offset(src_offset + ii);
                }
            }
            cur_col -= 1;
        }
    }

    //// Moves rows around and creates space for new rows
    fn create_row_space(&mut self, start: uint, count :uint){
        // The end must not be beyond capacity
        assert!(start + count <= self.xrows);
        if start >= self.rows {
            // Nothing to move.
            return;
        }
        let capacity = self.capacity() as int;
        // count rows starting from start row need to be
        // shifted by count.
        let mut cur_row = self.rows - 1;
        let ptr = self.ptr;
        // Number of rows to shift
        let rows_to_shift =  self.rows - (start + count - 1);
        let stride = self.stride();
        for _ in range(0, rows_to_shift){
            let dst_row = cur_row + count;
            let src_offset = self.cell_to_offset(cur_row, 0);
            let dst_offset = self.cell_to_offset(dst_row, 0);
            assert!(src_offset < capacity);
            assert!(dst_offset < capacity);
            for i in range(0, self.cols){
                let ii = (i*stride) as int;
                // Some validations
                assert!(src_offset + ii < capacity);
                assert!(dst_offset + ii < capacity);
                unsafe {
                    *ptr.offset(dst_offset + ii) = 
                    *ptr.offset(src_offset + ii);
                }
            }
            cur_row -= 1;
        }
    }
}

/// Views of a matrix
impl<T:MatrixElt> Matrix<T> {
    /// Creates a view on the matrix
    pub fn view(&self, start_row : int, start_col : int , num_rows: uint, num_cols : uint) -> MatrixView<T> {
        let result : MatrixView<T> = MatrixView::new(self, start_row, start_col, num_rows, num_cols);
        result
    }
}

/// These functions are available only for types which support
/// ordering [at least partial ordering for floating point numbers].
impl<T:MatrixElt+PartialOrd> Matrix<T> {
    /// Returns a column vector consisting of maximum over each row
    pub fn max_row_wise(&self) -> Matrix<T>{
        // Pick the first column
        let result = self.col(0);
        let pd = result.ptr;
        let ps = self.ptr;
        for r  in range(0, self.rows){
            let dst_offset = result.cell_to_offset(r, 0);
            for c in range (1, self.cols){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    let s = *ps.offset(src_offset);
                    let d = *pd.offset(dst_offset);
                    *pd.offset(dst_offset) = if s > d { s } else {d};

                }
            }
        }
        result
    }


   /// Returns a row vector consisting of maximum over each column
    pub fn max_col_wise(&self) -> Matrix<T>{
        // Pick the first row
        let result = self.row(0);
        let pd = result.ptr;
        let ps = self.ptr;
        for c in range (0, self.cols){
            let dst_offset = result.cell_to_offset(0, c);
            for r  in range(1, self.rows){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    let s = *ps.offset(src_offset);
                    let d = *pd.offset(dst_offset);
                    *pd.offset(dst_offset) = if s > d { s } else {d};

                }
            }
        }
        result
    }


    /// Returns a column vector consisting of minimum over each row
    pub fn min_row_wise(&self) -> Matrix<T>{
        // Pick the first column
        let result = self.col(0);
        let pd = result.ptr;
        let ps = self.ptr;
        for r  in range(0, self.rows){
            let dst_offset = result.cell_to_offset(r, 0);
            for c in range (1, self.cols){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    let s = *ps.offset(src_offset);
                    let d = *pd.offset(dst_offset);
                    *pd.offset(dst_offset) = if s < d { s } else {d};

                }
            }
        }
        result
    }


   /// Returns a row vector consisting of minimum over each column
    pub fn min_col_wise(&self) -> Matrix<T>{
        // Pick the first row
        let result = self.row(0);
        let pd = result.ptr;
        let ps = self.ptr;
        for c in range (0, self.cols){
            let dst_offset = result.cell_to_offset(0, c);
            for r  in range(1, self.rows){
                let src_offset = self.cell_to_offset(r, c);
                unsafe{
                    let s = *ps.offset(src_offset);
                    let d = *pd.offset(dst_offset);
                    *pd.offset(dst_offset) = if s < d { s } else {d};
                }
            }
        }
        result
    }

    // Returns the minimum scalar value
    pub fn min_scalar(&self) -> T{
        if self.is_empty(){
            fail!(EmptyMatrix.to_string());
        }
        let mut v = self.get(0, 0);
        let ps = self.ptr;
        for c in range(0, self.cols){
            for r in range(0, self.rows){
                let src_offset = self.cell_to_offset(r, c);
                let s = unsafe{*ps.offset(src_offset)};
                v = if s < v { s} else { v };
            }
        }
        v
    }

    // Returns the maximum scalar value
    pub fn max_scalar(&self) -> T{
        if self.is_empty(){
            fail!(EmptyMatrix.to_string());
        }
        let mut v = self.get(0, 0);
        let ps = self.ptr;
        for c in range(0, self.cols){
            for r in range(0, self.rows){
                let src_offset = self.cell_to_offset(r, c);
                let s = unsafe{*ps.offset(src_offset)};
                v = if s > v { s} else { v };
            }
        }
        v
    }    
}

/// These functions are available only for integer matrices
impl<T:MatrixElt+Int> Matrix<T> {

    /// Returns if an integer matrix is a logical matrix
    //// i.e. all cells are either 0s or 1s.
    pub fn is_logical(&self) -> bool {
        let z : T = Zero::zero();
        let o : T = One::one();
        for c in range(0, self.cols){
            for r in range(0, self.rows){
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
impl<T:MatrixElt+Float> Matrix<T> {
    /// Returns a matrix showing all the cells which are finite
    pub fn is_finite(&self) -> Matrix<u8>{
        let m : Matrix<u8> = Matrix::ones(self.rows, self.cols);
        for c in range(0, self.cols){
            for r in range(0, self.rows){
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
        for c in range(0, self.cols){
            for r in range(0, self.rows){
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






impl<T:MatrixElt> Index<uint,T> for Matrix<T> {
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



impl<T:MatrixElt> Matrix<T>{
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



/// Implementation of Clone interface
impl <T:MatrixElt> Clone for Matrix<T> {

    /// Creates a clone of the matrix
    fn clone(&self )-> Matrix<T> {
        let m : Matrix<T> = Matrix::new(self.rows, self.cols);
        unsafe{
            //TODO: this may not work if xrows and xcols are different.
            ptr::copy_memory(m.ptr, self.ptr as *const T, self.capacity());
        }
        m
    }
}

impl <T:MatrixElt> fmt::Show for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let ptr = self.ptr;
        try!(write!(f, "["));
        // Here we print row by row
        for r in range (0, self.rows) {
           try!(write!(f, "\n  "));
            for c in range (0, self.cols){
                let offset = self.cell_to_offset(r, c);
                let v = unsafe {*ptr.offset(offset)};
                try!(write!(f, "{} ", v));
            }
        }
        try!(write!(f, "\n]"));
        Ok(())
    }
}

/// Matrix addition support
impl<T:MatrixElt> ops::Add<Matrix<T>, Matrix<T>> for Matrix<T> {
    fn add(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            fail!(DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
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


/// Matrix subtraction support
impl<T:MatrixElt> ops::Sub<Matrix<T>, Matrix<T>> for Matrix<T>{
    fn sub(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            fail!(DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in range(0, n){
                let i = i_ as int;
                *pc.offset(i) = *pa.offset(i) - *pb.offset(i);
            }
        }
        result
    }
}


#[doc = "
This is not yet implemented.

We wish to support three different use cases ``m = m1 * m2`` 
where m1 and m2 are both matrices ``m = m1 * s`` 
where s is a scalar and ``m = s * m1``
where s is a scalar.

Currently, Rust doesn't
support multiple trait implementations for
the same type. It does get compiled, but one
faces the error: multiple applicable methods in scope.

The workaround is to define a separate trait in the
middle. 
    
"]
pub trait MatrixMul<Result> {
}


/// Matrix multiplication support
impl<T:MatrixElt> ops::Mul<Matrix<T>, Matrix<T>> for Matrix<T>{
    fn mul(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions match for multiplication
        if self.cols != rhs.rows{
            fail!(DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, rhs.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let zero : T = Zero::zero();
        unsafe {
            for r in range(0, self.rows){
                for c in range(0, rhs.cols){
                    let mut sum = zero;
                    for j in range(0, self.cols){
                        let lhs_offset = self.cell_to_offset(r, j);
                        let rhs_offset = rhs.cell_to_offset(j, c);
                        let term = *pa.offset(lhs_offset) * *pb.offset(rhs_offset);
                        sum = sum + term;
                    }
                    let dst_offset = result.cell_to_offset(r, c);
                    *pc.offset(dst_offset)  = sum;
                }
            }
        }
        result
    }
}


/// Matrix equality check support
impl<T:MatrixElt> cmp::PartialEq for Matrix<T>{
    fn eq(&self, other: &Matrix<T>) -> bool {
        let pa = self.ptr as *const  T;
        let pb = other.ptr as *const  T;
        for c in range(0, self.cols){
            for r in range(0, self.rows){
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
impl<T:MatrixElt> Matrix<T> {
    /// Multiplies matrices element by element
    pub fn mul_elt(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            fail!(DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in range(0, n){
                let i = i_ as int;
                *pc.offset(i) = *pa.offset(i) * *pb.offset(i);
            }
        }
        result
    }

    /// Divides matrices element by element
    pub fn div_elt(&self, rhs: &Matrix<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            fail!(DimensionsMismatch.to_string());
        }
        let result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.ptr;
        let pb = rhs.ptr;
        let pc = result.ptr;
        let n = self.capacity();
        unsafe{
            for i_ in range(0, n){
                let i = i_ as int;
                *pc.offset(i) = *pa.offset(i) / *pb.offset(i);
            }
        }
        result
    }
}

#[unsafe_destructor]
impl<T:MatrixElt> Drop for Matrix<T> {
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
 *   Utility functions for debugging of Matrix
 *
 *******************************************************/

impl<T:MatrixElt> Matrix<T> {
    pub fn print_state(&self){
        let capacity = self.capacity();
        let bytes = capacity * mem::size_of::<T>();
        println!("Rows: {}, Cols: {}, XRows : {}, XCols {} , Capacity: {}, Bytes; {}, Buffer: {:p}, End : {:p}", 
            self.rows, self.cols, self.xrows, self.xcols, 
            capacity, bytes, self.ptr, unsafe {
                self.ptr.offset(capacity as int)
            });
    }
}


/******************************************************
 *
 *   Private implementation of Matrix
 *
 *******************************************************/

impl<T:MatrixElt> Matrix<T> {
    /// Returns a slice into `self`.
    //#[inline]
    pub fn as_slice_<'a>(&'a self) -> &'a [T] {
        unsafe { mem::transmute(RawSlice { data: self.as_ptr(), len: self.capacity() }) }
    }

    /// Maps a cell index to actual offset in the internal buffer
    #[inline]
    fn cell_to_offset(&self, r : uint,  c: uint)-> int {
        (c * self.stride() + r) as int
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

    use  super::{Matrix, MatrixI64, MatrixF64};

    #[test]
    fn  create_mat0(){
        let m : MatrixI64 = Matrix::new(3, 4);
        assert_eq!(m.num_cells(), 12);
        assert_eq!(m.size(), (3, 4));
        let v : i64 = m.get(0, 0);
        assert_eq!(v, 0i64);
    }


    #[test]
    fn test_from_iter0(){
        let m : MatrixI64 = Matrix::from_iter(4, 4, range(1, 20));
        let b: Vec<i64> = range(1, 17).collect();
        assert!(m.as_slice_() == b.as_slice_());
        let m : MatrixI64 = Matrix::from_iter(4, 4, range(1, 16));
        assert_eq!(m.get(0, 0), 1);
        assert_eq!(m.get(2, 2), 11);
        let mut b: Vec<i64> = range(1, 16).collect();
        b.push(0);
        assert_eq!(m.as_slice_(), b.as_slice_());
    }

    #[test]
    fn test_index0(){
        let m : MatrixI64 = Matrix::from_iter(4, 4, range(1, 20));
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
        let m : MatrixI64 = Matrix::from_iter(4, 3, range(0, 12));
        let v1 = m.to_std_vec();
        let v2 : Vec<i64> = range(0, 12).collect();
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
        let m2 = m  + m;
        let v = vec![2i64, 2, 2, 2, 2, 2, 2, 2];
        assert_eq!(m2.to_std_vec(), v);
    }

    #[test]
    #[should_fail]
    fn test_sum_fail(){
        let m1 : MatrixI64 = Matrix::ones(4, 2);
        let m2 : MatrixI64 = Matrix::ones(3, 2);
        m1 + m2;
    }

    #[test]
    fn test_sub(){
        let m : MatrixI64 = Matrix::ones(4, 2);
        let m3 = m  + m + m;
        let m2 = m3 - m;  
        let v = vec![2i64, 2, 2, 2, 2, 2, 2, 2];
        assert_eq!(m2.to_std_vec(), v);
    }

    #[test]
    fn test_sub_float(){
        let m : MatrixF64 = Matrix::ones(4, 2);
        let m3 = m  + m + m;
        let m2 = m3 - m;  
        let v = vec![2f64, 2., 2., 2., 2., 2., 2., 2.];
        assert_eq!(m2.to_std_vec(), v);
    }
    #[test]
    #[should_fail]
    fn test_sub_fail(){
        let m1 : MatrixI64 = Matrix::ones(4, 2);
        let m2 : MatrixI64 = Matrix::ones(3, 2);
        m1 - m2;
    }


    #[test]
    fn test_mult(){
        let m1 : MatrixI64 = Matrix::from_iter(2, 2, range(0, 4));
        let m2 : MatrixI64 = Matrix::from_iter(2, 2, range(0, 4));
        let m3 = m1 * m2;
        let v = vec![2i64, 3, 6, 11];
        assert_eq!(m3.to_std_vec(), v);
    }

    #[test]
    fn test_eq(){
        let m1 : MatrixI64 = Matrix::from_iter(2, 2, range(0, 4));
        let m2 : MatrixI64 = Matrix::from_iter(2, 2, range(0, 4));
        assert_eq!(m1, m2);
        let v = vec![1.0f64, 2., 3., 4.];
        let m1 : MatrixF64 = Matrix::from_slice(2, 2, v.as_slice());
        let m2 : MatrixF64 = Matrix::from_slice(2, 2, v.as_slice());
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_extract_row(){
        let m1 : MatrixI64 = Matrix::from_iter(4, 4, range(0, 16));
        let m2  = m1.row(0);
        assert_eq!(m2.to_std_vec(), vec![0, 4, 8, 12]);
        assert_eq!(m2.num_rows() , 1);
        assert!(m2.is_row());
        assert_eq!(m2.num_cols() , m1.num_cols());
        assert_eq!(m2.num_cells() , m1.num_cols());
        assert_eq!(m1.row(1).to_std_vec(), vec![1, 5, 9, 13]);
        assert_eq!(m1.row(2).to_std_vec(), vec![2, 6, 10, 14]);
        assert_eq!(m1.row(3).to_std_vec(), vec![3, 7, 11, 15]);
        assert_eq!(m1.row(-1).to_std_vec(), vec![3, 7, 11, 15]);
        assert_eq!(m1.row(-2).to_std_vec(), vec![2, 6, 10, 14]);
        assert_eq!(m1.row(-6).to_std_vec(), vec![2, 6, 10, 14]);
    }


    #[test]
    fn test_extract_col(){
        let m1 : MatrixI64 = Matrix::from_iter(4, 4, range(0, 16));
        let m2  = m1.col(0);
        assert_eq!(m2.to_std_vec(), vec![0, 1, 2, 3]);
        assert!(!m2.is_row());
        assert!(m2.is_col());
        assert!(!m2.is_scalar());
        assert_eq!(m2.num_rows() , m1.num_rows());
        assert_eq!(m2.num_cells() , m1.num_rows());
        assert_eq!(m1.col(1).to_std_vec(), vec![4, 5, 6, 7]);
        assert_eq!(m1.col(2).to_std_vec(), vec![8, 9, 10, 11]);
        assert_eq!(m1.col(3).to_std_vec(), vec![12, 13, 14, 15]);
        assert_eq!(m1.col(-1).to_std_vec(), vec![12, 13, 14, 15]);
        assert_eq!(m1.col(-2).to_std_vec(), vec![8, 9, 10, 11]);
        assert_eq!(m1.col(-6).to_std_vec(), vec![8, 9, 10, 11]);
    }

    #[test]
    fn test_from_scalar(){
        let  m : MatrixI64  = Matrix::from_scalar(2);
        assert!(m.is_scalar());
        assert_eq!(m.to_std_vec(), vec![2]);
        assert_eq!(m.to_scalar(), 2);
    }
    #[test]
    fn test_sub_matrix(){
        let m  : MatrixI64 = Matrix::from_iter(4, 4, range(0, 16));
        let m1 = m.sub_matrix(0, 0, 2, 2);
        assert_eq!(m1.num_cells(), 4);
        assert_eq!(m1.num_rows(), 2);
        assert_eq!(m1.num_cols(), 2);
        assert_eq!(m1.to_std_vec(), vec![0, 1, 4, 5]);
        assert_eq!(m.sub_matrix(1, 0, 2, 2).to_std_vec(), vec![1, 2, 5, 6]);
        assert_eq!(m.sub_matrix(4, 4, 2, 2).to_std_vec(), vec![0, 1, 4, 5]);
        assert_eq!(m.sub_matrix(-1, -1, 2, 2).to_std_vec(), vec![15, 12, 3, 0]);
        assert_eq!(m.sub_matrix(-4, -4, 2, 2).to_std_vec(), vec![0, 1, 4, 5]);
    }

    #[test]
    fn test_rep_mat(){
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4));
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
        let m : MatrixF64 = Matrix::from_slice(2, 4, v.as_slice());
        let f = m.is_finite();
        assert_eq!(f.to_std_vec(), vec![1, 0, 0, 0, 0, 
            1, 1, 1]);
        let f = m.is_infinite();
        assert_eq!(f.to_std_vec(), vec![0, 0, 0, 1, 1, 
            0, 0, 0]);
    }


    #[test]
    fn test_is_logical(){
        let m : MatrixI64 = Matrix::from_iter(4, 4, range(0, 16).map(|x| x % 2));
        assert!(m.is_logical());
        let m = m + m;
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
        let m : MatrixI64 = Matrix::from_slice(4, 4, v.as_slice());
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

        assert_eq!(m.min_scalar(), 0);
        assert_eq!(m.max_scalar(), 19);
    }

    #[test]
    fn test_mul_elt(){
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4));
        let m2 = m.mul_elt(&m);
        let m3  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4).map(|x| x*x));
        assert_eq!(m2, m3);
    }


    #[test]
    fn test_div_elt(){
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(1, 20));
        let m2 = m + m;
        let m3 = m2.div_elt(&m);
        assert_eq!(m3.to_std_vec(), vec![2,2,2,2]);
    }

    #[test]
    fn test_add_scalar (){
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4));
        let m2 = m.add_scalar(2);
        assert_eq!(m2.to_std_vec(), vec![2, 3, 4, 5]);
    }

    #[test]
    fn test_mul_scalar (){
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4));
        let m2 = m.mul_scalar(2);
        assert_eq!(m2.to_std_vec(), vec![0, 2, 4, 6]);
    }

    #[test]
    fn test_div_scalar (){
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4).map(|x| x * 3));
        let m2 = m.div_scalar(3);
        assert_eq!(m2.to_std_vec(), vec![0, 1, 2, 3]);
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
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4));
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
        let m2 = m.mul_scalar(2);
        assert!(m.is_square());
        let m4 = m2.pow(4);
        let m16 = m.mul_scalar(16);
        assert_eq!(m4, m16); 
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4));
        assert!(m.is_square());
        let m3 = m.pow(3);
        assert!(m3.is_square());
        assert_eq!(m3.to_std_vec(), vec![6, 11, 22, 39]);
        assert_eq!(m.pow(1).to_std_vec(), vec![0, 1, 2, 3]);
        assert_eq!(m.pow(2), m * m);
        assert_eq!(m.pow(10), m * m * m * m * m * m * m * m * m * m);
    }

    #[test]
    fn test_transpose(){
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4));
        assert_eq!(m.to_std_vec(), vec![0, 1, 2, 3]);
        assert_eq!(m.transpose().to_std_vec(), vec![0, 2, 1, 3]);
        assert_eq!(m.transpose().transpose().to_std_vec(), vec![0, 1, 2, 3]);
        let m  : MatrixI64 = Matrix::from_iter(2, 3,  range(0, 10));
        assert_eq!(m.transpose().to_std_vec(), vec![
            0, 2, 4, 1, 3, 5]);
        let m4 :  MatrixI64 = Matrix::from_iter(2, 5, range(9, 100));
        let m5 = m4.transpose();
        println!("m4: {}", m4);
        println!("m5: {}", m5);
        assert_eq!(m5.transpose(), m4);
    }

    #[test]
    fn test_unary_minus(){
        let m  : MatrixI64 = Matrix::from_iter(2, 2,  range(0, 4));
        let z : MatrixI64 = Matrix::zeros(2,2);
        let m2 = m.unary_minus();
        let m3 = z - m;
        assert_eq!(m2, m3);
    }


    #[test]
    fn test_diag_from_vector(){
        let v  : MatrixI64 = Matrix::from_iter(4, 1, range(20, 30));
        assert!(v.is_vector());
        let m = Matrix::diag_from_vec(&v);
        assert!(!m.is_empty());
        assert!(!m.is_vector());
        assert!(m.is_diagonal());
        assert_eq!(m.num_cells(), 16);
        let mut m2 : MatrixI64 = Matrix::zeros(4, 4);
        m2.set(0, 0, 20);
        m2.set(1, 1, 21);
        m2.set(2, 2, 22);
        m2.set(3, 3, 23);
        assert_eq!(m, m2);
    }

    #[test]
    fn test_diagonal(){
        let m  : MatrixI64 = Matrix::from_iter(4, 5, range(10, 30));
        let v = m.diagonal();
        assert!(v.is_vector());
        assert_eq!(v.num_cells(), 4);
        let v2 : MatrixI64 = Matrix::from_slice(4, 1, vec![10, 15, 20, 25].as_slice());
        assert_eq!(v, v2);
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
        let m  : MatrixI64 = Matrix::from_iter(4, 5, range(10, 30));
        let mut r = m.row_iter(0);
        let v : Vec<i64> = r.collect();
        assert_eq!(v, vec![10, 14, 18, 22, 26]);
        let mut r = m.col_iter(2);
        let v : Vec<i64> = r.collect();
        assert_eq!(v, vec![18, 19, 20, 21]);
        let m  : MatrixI64 = Matrix::from_iter(3, 2, range(10, 30));
        let mut r = m.cell_iter();
        let v : Vec<i64> = r.collect();
        assert_eq!(v, vec![10, 11, 12, 13, 14, 15]);
    }

    #[test]
    fn test_add_columns(){
        let mut m1 :  MatrixI64 = Matrix::from_iter(2, 3, range(11, 100));
        let m2 :  MatrixI64 = Matrix::from_iter(2, 4, range(11, 100));
        let m3 : MatrixI64  = Matrix::from_iter(2, 1, range(17, 100));
        let m4 :  MatrixI64 = Matrix::from_iter(2, 5, range(9, 100));
        let m5 : MatrixI64  = Matrix::from_iter(2, 1, range(9, 100));
        println!("{}", m1);
        m1.append_columns(&m3);
        println!("{}", m1);
        assert_eq!(m1, m2);
        m1.prepend_columns(&m5);
        println!("{}", m1);
        assert_eq!(m1, m4);
        let m6 : MatrixI64  = Matrix::from_iter(2, 1, range(5, 100));
        m1.prepend_columns(&m6);
        println!("{}", m1);
        let m7 : MatrixI64  = Matrix::from_iter(2, 1, range(7, 100));
        m1.insert_columns(1, &m7);
        println!("{}", m1);
    }


    #[test]
    fn test_add_rows(){
        let mut m1 :  MatrixI64 = Matrix::from_iter(2, 3, range(11, 100)).transpose();
        let m2 :  MatrixI64 = Matrix::from_iter(2, 4, range(11, 100)).transpose();
        let m3 : MatrixI64  = Matrix::from_iter(2, 1, range(17, 100)).transpose();
        let m4 :  MatrixI64 = Matrix::from_iter(2, 5, range(9, 100)).transpose();
        let m5 : MatrixI64  = Matrix::from_iter(2, 1, range(9, 100)).transpose();
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
        let m6 : MatrixI64  = Matrix::from_iter(2, 1, range(5, 100)).transpose();
        m1.prepend_rows(&m6);
        println!("{}", m1);
        let m7 : MatrixI64  = Matrix::from_iter(2, 1, range(7, 100)).transpose();
        m1.insert_rows(1, &m7);
        println!("{}", m1);
    }

}
