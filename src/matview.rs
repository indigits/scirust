// Standard library imports
use std::mem;
use std::ops;

// srmat imports
use matelt::{MatrixElt};
use matrix::{Matrix};
use materr::*;
use discrete::*;

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
    pub fn new(m : &Matrix<T>, start_row : int, start_col : int , num_rows: uint, num_cols : uint) -> MatrixView<T> {
        let r = mod_n(start_row, m.num_rows() as int);        
        let c = mod_n(start_col, m.num_cols() as int);
        let row_skip = 1u;
        let col_skip = 1u;
        assert!(r + row_skip*num_rows <= m.num_rows());
        assert!(c + col_skip*num_cols <= m.num_cols());
        let result : MatrixView<T> = MatrixView{
            m : m,
            start_row : r,
            start_col : c, 
            rows: num_rows,
            cols : num_cols,
            row_skip : row_skip,
            col_skip : col_skip
        };
        result
    }
    //// Returns the number of rows in the view
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
        assert!(r < self.rows);
        assert!(c < self.cols);
        let ptr = self.m.as_ptr();
        let offset = self.cell_to_offset(r, c);
        assert!((offset as uint) < self.m.capacity());
        unsafe {
            *ptr.offset(offset)
        }
    }

    /// Sets an element in the view
    #[inline]
    pub fn set(&mut self, r : uint, c : uint, value : T) {
        assert!(r < self.rows);
        assert!(c < self.cols);
        let ptr = self.m.as_ptr();
        // I know more than the compiler
        // I am allowing modification of the underlying buffer
        let ptr : *mut T = unsafe { mem::transmute(ptr) };
        let offset = self.cell_to_offset(r, c);
        assert!((offset as uint) < self.m.capacity());
        unsafe {
            *ptr.offset(offset) = value;
        }
    }

    /// Converts an index to cell address (row, column)
    #[inline]
    pub fn index_to_cell(&self, index : uint) -> (uint, uint){
        assert!(index < self.num_cells());
        let c = index / self.rows;
        let r = index - c*self.rows;
        (r, c)
    }

    /// Converts a cell address to an index (r, c) to index
    #[inline]
    pub fn cell_to_index(&self, r : uint,  c: uint) -> uint{
        assert!(r < self.rows);
        assert!(c < self.cols);
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
}

impl<'a, T:MatrixElt> MatrixView<'a, T> {
}


impl<'a, T:MatrixElt> MatrixView<'a, T> {
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
}

