// Standard library imports
use std::mem;
use std::ops;
use std::num;
use std::fmt;
use std::ptr;

// srmat imports
use matrix::element::{MatrixElt};
use matrix::matrix::{Matrix};
use matrix::error::*;
//use discrete::*;

#[doc = "
Defines a view on a matrix. 

A view on a matrix is a subset of 
chosen rows and columns.


"]
pub struct MatrixView<'a, T:'a+MatrixElt>{
    // Reference to the associated matrix
    m : &'a Matrix<T>,
    // start row
    start_row : uint,
    // number or rows
    rows  : uint, 
    // skip between rows (by default 1)
    row_skip : uint,
    // start column
    start_col : uint,
    // Number of columns 
    cols : uint,
    // Skip between columns (by default 1)
    col_skip : uint
}


/// Static functions for creating  a view
impl<'a, T:MatrixElt> MatrixView<'a, T> {
    pub fn new(m : &Matrix<T>, start_row : uint, start_col : uint , num_rows: uint, num_cols : uint) -> MatrixView<T> {
        let row_skip = 1u;
        let col_skip = 1u;
        debug_assert!(start_row + row_skip*num_rows <= m.num_rows());
        debug_assert!(start_col + col_skip*num_cols <= m.num_cols());
        let result : MatrixView<T> = MatrixView{
            m : m,
            start_row : start_row,
            start_col : start_col, 
            rows: num_rows,
            cols : num_cols,
            row_skip : row_skip,
            col_skip : col_skip
        };
        result
    }

    /// Returns the start row
    pub fn start_row(&self) -> uint{
        self.start_row
    } 

    /// Returns the number of rows in the view
    pub fn num_rows(&self) -> uint {
        self.rows
    }

    /// Returns the number of columns in the view
    pub fn num_cols(&self) -> uint {
        self.cols
    }


    /// Returns the size of view in an (r, c) tuple
    pub fn size (&self)-> (uint, uint){
        (self.rows, self.cols)
    }
    /// Returns the number of cells in view
    pub fn num_cells(&self)->uint {
        self.rows * self.cols
    }
    /// Indicates if the view is a row vector
    pub fn is_row(&self) -> bool {
        self.rows == 1
    }

    /// Indicates if the view is a column vector
    pub fn is_col(&self) -> bool {
        self.cols == 1
    }

    /// Indicates if the view is a scalar actually
    pub fn is_scalar(&self) -> bool {
        self.num_cells() == 1
    }

    /// Indicates if the view is a vector
    pub fn is_vector(&self) -> bool {
        (self.rows == 1) ^ (self.cols == 1)
    } 

    /// Indicates if the view is empty
    pub fn is_empty(&self) -> bool {
        self.rows * self.cols == 0
    }

    /// Indicates if the view is square
    pub fn is_square(&self) -> bool {
        self.rows == self.cols
    }

    /// Gets an element in the view
    #[inline]
    pub fn get(&self, r : uint, c : uint) -> T  {
        debug_assert!(r < self.rows);
        debug_assert!(c < self.cols);
        let ptr = self.m.as_ptr();
        let offset = self.cell_to_offset(r, c);
        debug_assert!((offset as uint) < self.m.capacity());
        unsafe {
            *ptr.offset(offset)
        }
    }

    /// Sets an element in the view
    #[inline]
    pub fn set(&mut self, r : uint, c : uint, value : T) {
        debug_assert!(r < self.rows);
        debug_assert!(c < self.cols);
        let ptr = self.m.as_ptr();
        // I know more than the compiler
        // I am allowing modification of the underlying buffer
        let ptr : *mut T = unsafe { mem::transmute(ptr) };
        let offset = self.cell_to_offset(r, c);
        debug_assert!((offset as uint) < self.m.capacity());
        unsafe {
            *ptr.offset(offset) = value;
        }
    }

    /// Converts an index to cell address (row, column)
    #[inline]
    pub fn index_to_cell(&self, index : uint) -> (uint, uint){
        debug_assert!(index < self.num_cells());
        let c = index / self.rows;
        let r = index - c*self.rows;
        (r, c)
    }

    /// Converts a cell address to an index (r, c) to index
    #[inline]
    pub fn cell_to_index(&self, r : uint,  c: uint) -> uint{
        debug_assert!(r < self.rows);
        debug_assert!(c < self.cols);
        c * self.rows + r
    }

    /******************************************************
     *
     *   Private implementation of MatrixView
     *
     *******************************************************/

    /// Maps a cell index to actual offset in the internal buffer
    #[inline]
    fn cell_to_offset(&self, r : uint,  c: uint)-> int {
        let r = self.start_row + r * self.row_skip;
        let c = self.start_col + c * self.col_skip;
        (c * self.m.stride() + r) as int
    } 

}



/// Functions to construct new views out of a view and other conversions
impl<'a, T:MatrixElt> MatrixView<'a, T> {
    /// Converts the view to vector from standard library
    pub fn to_std_vec(&self) -> Vec<T> {
        let mut vec: Vec<T> = Vec::with_capacity(self.num_cells());
        // We iterate over elements in matrix and push in the vector
        let ptr = self.m.as_ptr();
        for c in range(0, self.cols){
            for r in range (0, self.rows){
                let offset = self.cell_to_offset(r, c);
                vec.push(unsafe{*ptr.offset(offset)});
            }
        } 
        vec
    }

    /// Returns the view as a new matrix.
    /// Creates a copy of the data.
    pub fn to_matrix(&self,) -> Matrix<T> {
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pd = result.as_mut_ptr();
        let ps = self.m.as_ptr();
        for c in range(0, self.cols) {
            for r in range(0, self.rows) {
                let src_offset = self.cell_to_offset(r, c);
                let dst_offset = result.cell_to_offset(r, c);
                unsafe{
                    *pd.offset(dst_offset) = *ps.offset(src_offset);
                }
            }
        }
        result
    }

}

impl<'a, T:MatrixElt> MatrixView<'a, T> {
}


impl<'a, T:MatrixElt+PartialOrd> MatrixView<'a, T> {
    // Returns the minimum scalar value
    pub fn min_scalar(&self) -> (T, uint, uint){
        if self.is_empty(){
            fail!(EmptyMatrix.to_string());
        }
        let mut v = self.get(0, 0);
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.m.as_ptr();
        for c in range(0, self.cols){
            for r in range(0, self.rows){
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

    // Returns the maximum scalar value
    pub fn max_scalar(&self) -> (T, uint, uint){
        if self.is_empty(){
            fail!(EmptyMatrix.to_string());
        }
        let mut v = self.get(0, 0);
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.m.as_ptr();
        for c in range(0, self.cols){
            for r in range(0, self.rows){
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

impl<'a, T:MatrixElt+PartialOrd+Signed> MatrixView<'a, T> {

    // Returns the absolute minimum scalar value
    pub fn abs_min_scalar(&self) -> (T, uint, uint){
        if self.is_empty(){
            fail!(EmptyMatrix.to_string());
        }
        let mut v = num::abs(self.get(0, 0));
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.m.as_ptr();
        for c in range(0, self.cols){
            for r in range(0, self.rows){
                let src_offset = self.cell_to_offset(r, c);
                let s = num::abs(unsafe{*ps.offset(src_offset)});
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
    pub fn abs_max_scalar(&self) -> (T, uint, uint){
        if self.is_empty(){
            fail!(EmptyMatrix.to_string());
        }
        let mut v = num::abs(self.get(0, 0));
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.m.as_ptr();
        for c in range(0, self.cols){
            for r in range(0, self.rows){
                let src_offset = self.cell_to_offset(r, c);
                let s = num::abs(unsafe{*ps.offset(src_offset)});
                if s > v { 
                    v  = s;
                    rr = r;
                    cc = c;
                }
            }
        }
        (v, rr, cc)
    }    


    /******************************************************
     *
     *   Elementary row / column operations.
     *
     *******************************************************/


    /// Row switching.
    pub fn ero_switch(&mut self, 
        i :  uint,
        j : uint
        )-> &mut MatrixView<'a, T> {
        debug_assert! (i  < self.rows);
        debug_assert! (j  < self.rows);
        let ptr = self.m.as_ptr();
        // I am allowing modification of the underlying buffer
        let ptr : *mut T = unsafe { mem::transmute(ptr) };
        for c in range(0, self.cols){
            let offset_a = self.cell_to_offset(i, c);
            let offset_b = self.cell_to_offset(j, c);
            unsafe {
                ptr::swap(ptr.offset(offset_a), ptr.offset(offset_b));
            }
        }
        self
    }

    /// Row scaling by a factor.
    pub fn ero_scale(&mut self, 
        r :  uint, 
        scale : T
        )-> &mut MatrixView<'a, T> {
        debug_assert! (r  < self.rows);
        let ptr = self.m.as_ptr();
        // I am allowing modification of the underlying buffer
        let ptr : *mut T = unsafe { mem::transmute(ptr) };
        for c in range(0, self.cols){
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
    /// The j-th row can be outside the view also.
    /// This is the row relative to the start of the view.
    pub fn ero_scale_add(&mut self, 
        i :  uint, 
        j :  int, 
        scale : T
        )-> &mut MatrixView<'a, T> {
        debug_assert! (i  < self.rows);
        let m = self.m;
        // Compute j-th row in m (by doing offset)
        let j = j + (self.start_row as int);
        debug_assert! (j  >= 0);
        let j = j as uint;
        debug_assert!(j < m.num_rows());
        let ptr = m.as_ptr();
        // I am allowing modification of the underlying buffer
        let ptr : *mut T = unsafe { mem::transmute(ptr) };
        for c in range(0, self.cols){
            // i-th row from the view
            let offset_a = self.cell_to_offset(i, c);
            // j-th row from the matrix
            let offset_b = m.cell_to_offset(j, c * self.col_skip + self.start_col);
            unsafe {
                let va = *ptr.offset(offset_a);
                let vb = *ptr.offset(offset_b);
                *ptr.offset(offset_a) = va + scale * vb;
            }
        }
        self
    }

}

/// View + View =  Matrix addition
impl<'a, 'b, T:MatrixElt> ops::Add<MatrixView<'b, T>, Matrix<T>> for MatrixView<'a, T> {
    fn add(&self, rhs: &MatrixView<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            fail!(DimensionsMismatch.to_string());
        }
        let mut result : Matrix<T> = Matrix::new(self.rows, self.cols);
        let pa = self.m.as_ptr();
        let pb = rhs.m.as_ptr();
        let pc = result.as_mut_ptr();
        for c in range (0, self.cols){
            for r in range(0, self.rows){
                let dst_offset = result.cell_to_offset(r, c);
                let a_offset  = self.cell_to_offset(r, c);
                let b_offset = rhs.cell_to_offset(r, c);
                unsafe{
                    *pc.offset(dst_offset) = *pa.offset(a_offset) + *pb.offset(b_offset);
                }
            }
        }
        result
    }
}


impl <'a, T:MatrixElt> fmt::Show for MatrixView<'a, T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let ptr = self.m.as_ptr();
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


/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/


#[cfg(test)]
mod test{
    use matrix::{Matrix, MatrixI64};

    #[test]
    fn test_basic(){
        let m1 :  MatrixI64 = Matrix::from_iter(10, 8, range(1, 100));
        let mut v1 = m1.view(2, 3, 4, 4);
        assert_eq!(v1.size(), (4u, 4u));
        assert_eq!(v1.num_rows(), 4);
        assert_eq!(v1.num_cols(), 4);
        assert_eq!(v1.is_scalar(), false);
        assert_eq!(v1.is_vector(), false);
        assert_eq!(v1.is_row(), false);
        assert_eq!(v1.is_col(), false);
        assert_eq!(v1.get(1,1), 44);
        v1.set(1,1, 300);
        assert_eq!(v1.get(1,1), 300);
    }
    #[test]
    fn test_view_multiple(){
        let m1 :  MatrixI64 = Matrix::from_iter(10, 8, range(1, 100));
        // We can create multiple views easily.
        let v1 = m1.view(1, 1, 4, 4);
        let v2 = m1.view(2, 2, 4, 4);
        let v3 = m1.view(2, 3, 4, 4);
        assert_eq!(v1.get(0, 0), 12);
        assert_eq!(v2.get(0, 0), 23);
        assert_eq!(v3.get(0, 0), 33);
    }

    #[test]
    fn test_view_matrix_mutability(){
        let mut m1 :  MatrixI64 = Matrix::from_iter(10, 8, range(1, 100));
        {
            let v1 = m1.view(2, 2, 4, 4);
            // The following line doesn't compile since we have borrowed an
            // immutable reference of m1 in v1.
            //m1.set(0, 0, 2);
            assert_eq!(v1.get(0, 0), 23);
        }
        m1.set(0, 0, 2);
    }

    #[test]
    fn test_cell_index_mapping(){
        let rows = 20;
        let cols = 10;
        let m : MatrixI64 = Matrix::new(rows, cols);
        let vrows = 4;
        let vcols = 3;
        let v = m.view(1, 1, vrows, vcols);
        assert_eq!(v.index_to_cell(3 + 2*vrows), (3, 2));
        assert_eq!(v.cell_to_index(3, 2), 3 + 2*vrows);
        assert_eq!(v.cell_to_index(3, 2), 11);
        assert_eq!(v.index_to_cell(11), (3, 2));
    }

    #[test]
    fn test_to_std_vec(){
        let m :  MatrixI64 = Matrix::from_iter(10, 8, range(1, 100));
        let vrows = 4;
        let vcols = 3;
        let mv = m.view(1, 1, vrows, vcols);
        let vv : Vec<i64> = vec![12, 13, 14, 15, 
        22, 23, 24, 25,
        32, 33, 34, 35];
        assert_eq!(mv.to_std_vec(), vv);
    }

    #[test]
    fn test_addition(){
        let m :  MatrixI64 = Matrix::from_iter(10, 10, range(1, 200));
        let v1   = m.view(2, 2, 2, 2); // 23 , 24 , 33, 34
        let v2 = m.view(1, 1, 2, 2);  // 12, 13, 22, 33
        let m2 = v1 + v2; // 
        let m3 : MatrixI64 = Matrix::from_slice(2, 2, vec![35, 37, 55, 57].as_slice());
        assert_eq!(m2, m3);
    }


    #[test]
    fn test_min_max(){
        let m :  MatrixI64 = Matrix::from_iter(20, 20, range(-100, 400));
        let v1   = m.view(2, 2, 6, 6);
        println!("v1 : {}", v1);
        assert_eq!(v1.max_scalar(), (47, 5, 5));
        assert_eq!(v1.min_scalar(), (-58, 0, 0));
        let (abs_max, rr, cc) = v1.abs_max_scalar();
        assert_eq!(abs_max, 58);
        assert_eq!((rr, cc), (0, 0));
        let (abs_min, rr, cc) = v1.abs_min_scalar();
        assert_eq!(abs_min, 2);
        assert_eq!((rr, cc), (0, 3));
    }

}

