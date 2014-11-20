// Standard library imports
use std::mem;
use std::ops;
use std::fmt;
use std::num::Float;
//use std::ptr;

// local imports
use entry::{Entry, One, Zero};
use number::{Number};
use error::SRError;
use matrix::matrix::{Matrix};
use matrix::traits::{Shape, NumberMatrix, 
    Strided,
    StridedNumberMatrix,
    StridedFloatMatrix,
    Introspection, 
    MatrixBuffer, Extraction};

//use discrete::*;

#[doc = "
Defines a view on a matrix. 

A view on a matrix is a subset of 
chosen rows and columns.


"]
pub struct MatrixView<'a, T:'a+Entry>{
    // Reference to the associated matrix
    m : &'a Matrix<T>,
    // start row
    start_row : uint,
    // number or rows
    rows  : uint, 
    // start column
    start_col : uint,
    // Number of columns 
    cols : uint,
}


/// Static functions for creating  a view
impl<'a, T:Entry> MatrixView<'a, T> {
    pub fn new(m : &Matrix<T>, start_row : uint, start_col : uint , num_rows: uint, num_cols : uint) -> MatrixView<T> {
        debug_assert!(start_row + num_rows <= m.num_rows());
        debug_assert!(start_col + num_cols <= m.num_cols());
        let result : MatrixView<T> = MatrixView{
            m : m,
            start_row : start_row,
            start_col : start_col, 
            rows: num_rows,
            cols : num_cols
        };
        result
    }

}

///Basic methods for a view
impl<'a, T:Entry> MatrixView<'a, T> {
    /// Returns the start row
    #[inline]
    pub fn start_row(&self) -> uint{
        self.start_row
    } 
    /// Returns the start column
    #[inline]
    pub fn start_col(&self) -> uint{
        self.start_col
    }
    /// Returns the underlying matrix reference 
    #[inline]
    pub fn matrix(&self)-> &'a Matrix<T>{
        self.m
    }
}


///Basic methods for a view of a matrix of  numbers
impl<'a, T:Number> MatrixView<'a, T> {

    /// Copies data from other view
    // TODO: write tests
    pub fn copy_from(&mut self, rhs: &MatrixView<T>){
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
        }        
        let pd : *mut T = unsafe { mem::transmute(self.m.as_ptr()) };
        let ps = rhs.m.as_ptr();
        for c in range (0, self.cols){
            for r in range(0, self.rows){
                let src_offset  = rhs.cell_to_offset(r, c);
                let dst_offset = self.cell_to_offset(r, c);
                unsafe{
                    *pd.offset(dst_offset) = *ps.offset(src_offset);
                }
            }
        }
    }

    /// Copies data from other view with scale factor
    // TODO: write tests
    pub fn copy_scaled_from(&mut self, rhs: &MatrixView<T>, scale: T){
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
        }        
        let pd : *mut T = unsafe { mem::transmute(self.m.as_ptr()) };
        let ps = rhs.m.as_ptr();
        for c in range (0, self.cols){
            for r in range(0, self.rows){
                let src_offset  = rhs.cell_to_offset(r, c);
                let dst_offset = self.cell_to_offset(r, c);
                unsafe{
                    *pd.offset(dst_offset) = scale * *ps.offset(src_offset);
                }
            }
        }
    }


    /******************************************************
     *
     *   Private implementation of MatrixView
     *
     *******************************************************/

}

/// Strided buffer 
impl <'a, T:Entry> Strided for MatrixView<'a, T> {
    /// Returns the number of actual memory elements 
    /// per column stored in the memory
    fn stride (&self)->uint {
        self.matrix().stride()
    }
}

/// Implement Buffer API for matrix view
impl <'a, T:Entry> MatrixBuffer<T> for MatrixView<'a, T> {

    /// Returns an unsafe pointer to the matrix's 
    /// buffer.
    #[inline]
    fn as_ptr(&self)-> *const T{
        self.matrix().as_ptr()
    }

    /// Returns a mutable unsafe pointer to
    /// the matrix's underlying buffer
    #[inline]
    fn as_mut_ptr(&mut self) -> *mut T{
        let p  = self.matrix().as_ptr();
        let ptr : *mut T = unsafe { mem::transmute(p) };
        ptr
    }

    /// Maps a cell index to actual offset in the internal buffer
    #[inline]
    fn cell_to_offset(&self, r : uint,  c: uint)-> int {
        let r = self.start_row + r;
        let c = self.start_col + c;
        (c * self.m.stride() + r) as int
    } 

    /// Returns the offset of the first cell in the buffer
    #[inline]
    fn start_offset(&self) -> int {
        (self.start_col() * self.stride() + self.start_row()) as int
    }

}


/// Implement extraction API for matrix view 
impl <'a, T:Number> Extraction<T> for MatrixView<'a, T> {

}




/// Implementation of common matrix methods
impl <'a, T:Entry> Shape<T> for MatrixView<'a, T> {
    /// Returns the number of rows in the view
    fn num_rows(&self) -> uint {
        self.rows
    }

    /// Returns the number of columns in the view
    fn num_cols(&self) -> uint {
        self.cols
    }


    /// Returns the size of view in an (r, c) tuple
    fn size (&self)-> (uint, uint){
        (self.rows, self.cols)
    }

    /// Returns the number of cells in view
    fn num_cells(&self)->uint {
        self.rows * self.cols
    }


    /// Gets an element in the view
    #[inline]
    fn get(&self, r : uint, c : uint) -> T  {
        debug_assert!(r < self.rows);
        debug_assert!(c < self.cols);
        let ptr = self.m.as_ptr();
        let offset = self.cell_to_offset(r, c);
        debug_assert!((offset as uint) < self.m.capacity());
        unsafe {
            // TODO : Optimize this
            ptr.offset(offset).as_ref().unwrap().clone()
        }
    }

    /// Sets an element in the view
    #[inline]
    fn set(&mut self, r : uint, c : uint, value : T) {
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

}

/// Implementation of methods related to matrices of numbers
impl <'a, T:Number> NumberMatrix<T> for MatrixView<'a, T> {
    /// Returns if the matrix is an identity matrix
    fn is_identity(&self) -> bool {
        let o : T = One::one();
        let z  : T = Zero::zero();
        let ptr = self.m.as_ptr();
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
    fn is_diagonal(&self) -> bool {
        let z  : T = Zero::zero();
        let ptr = self.m.as_ptr();
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

    /// Returns if the matrix is lower triangular 
    fn is_lt(&self) -> bool {
        let z  : T = Zero::zero();
        let ptr = self.m.as_ptr();
        for c in range(0, self.cols){
            for r in range (0, c){
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
        let ptr = self.m.as_ptr();
        for c in range(0, self.cols){
            for r in range (c+1, self.rows){
                let offset = self.cell_to_offset(r, c);
                let v = unsafe {*ptr.offset(offset)};
                println!("r: {}, c: {}, v: {}", r, c, v);
                if v != z {
                    return false;
                }
            }
        } 
        true
    }

    fn trace(&self) -> T{
        if self.is_empty() {
            return Zero::zero()
        }
        let stride = self.stride() as int;
        let mut offset = self.start_offset();
        let ptr = self.as_ptr();
        let mut result = unsafe {*ptr.offset(offset)};
        for i in range(1, self.smaller_dim()){
            offset += stride;
            result = result + unsafe{*ptr.offset(offset + i as int)};
        }
        result
    }

}

impl<'a, T:Number> StridedNumberMatrix<T> for MatrixView<'a, T> {
}


impl<'a, T:Number+Float> StridedFloatMatrix<T> for MatrixView<'a, T> {
}


/// Functions to construct new views out of a view and other conversions
impl<'a, T:Number> MatrixView<'a, T> {

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

/// Introspection support
impl<'a, T:Number> Introspection for MatrixView<'a, T> {
    /// This is a view inside a matrix
    fn is_matrix_view_type(&self) -> bool {
        true
    }
}


impl<'a, T:Number+PartialOrd> MatrixView<'a, T> {
    // Returns the minimum scalar value
    pub fn min_scalar(&self) -> (T, uint, uint){
        if self.is_empty(){
            panic!(SRError::EmptyMatrix.to_string());
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
            panic!(SRError::EmptyMatrix.to_string());
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


/// View + View =  Matrix addition
impl<'a, 'b, T:Number> ops::Add<MatrixView<'b, T>, Matrix<T>> for MatrixView<'a, T> {
    fn add(&self, rhs: &MatrixView<T>) -> Matrix<T> {
        // Validate dimensions are same.
        if self.size() != rhs.size(){
            panic!(SRError::DimensionsMismatch.to_string());
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


impl <'a, T:Number> fmt::Show for MatrixView<'a, T> {
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
    use matrix::traits::*;
    use matrix::eo::eo_traits::*;
    use matrix::matrix::*;
    use matrix::constructors::*;

    #[test]
    fn test_basic(){
        let m1 :  MatrixI64 = Matrix::from_iter_cw(10, 8, range(1, 100));
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
        assert!(m1.is_standard_matrix_type());
        assert!(!m1.is_matrix_view_type());
        assert!(!v1.is_standard_matrix_type());
        assert!(v1.is_matrix_view_type());
    }
    #[test]
    fn test_view_multiple(){
        let m1 :  MatrixI64 = Matrix::from_iter_cw(10, 8, range(1, 100));
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
        let mut m1 :  MatrixI64 = Matrix::from_iter_cw(10, 8, range(1, 100));
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
        let m :  MatrixI64 = Matrix::from_iter_cw(10, 8, range(1, 100));
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
        let m :  MatrixI64 = Matrix::from_iter_cw(10, 10, range(1, 200));
        let v1   = m.view(2, 2, 2, 2); // 23 , 24 , 33, 34
        let v2 = m.view(1, 1, 2, 2);  // 12, 13, 22, 33
        let m2 = v1 + v2; // 
        let m3 : MatrixI64 = Matrix::from_slice_cw(2, 2, vec![35, 37, 55, 57].as_slice());
        assert_eq!(m2, m3);
    }


    #[test]
    fn test_min_max(){
        let m :  MatrixI64 = Matrix::from_iter_cw(20, 20, range(-100, 400));
        let v1   = m.view(2, 2, 6, 6);
        println!("v1 : {}", v1);
        assert_eq!(v1.max_scalar(), (47, 5, 5));
        assert_eq!(v1.min_scalar(), (-58, 0, 0));
        let (abs_max, rr, cc) = v1.max_abs_scalar();
        assert_eq!(abs_max, 58);
        assert_eq!((rr, cc), (0, 0));
        let (abs_min, rr, cc) = v1.min_abs_scalar();
        assert_eq!(abs_min, 2);
        assert_eq!((rr, cc), (0, 3));
    }


    #[test]
    fn test_extract_row(){
        let m :  MatrixI64 = Matrix::from_iter_cw(20, 20, range(-100, 400));
        let v1   = m.view(2, 2, 6, 6);
        println!("v1 : {}", v1);
        let r1 = v1.row(0);
        let v2 = m.view(2,2, 1, 6);
        assert_eq!(v2.to_matrix(), r1);
    }

    #[test]
    fn test_extract_col(){
        let m :  MatrixI64 = from_range_rw_i64(20, 20, -100, 400);
        let v1   = m.view(2, 2, 6, 6);
        println!("v1 : {}", v1);
        let c1 = v1.col(0);
        let v2 = m.view(2,2, 6, 1);
        assert_eq!(v2.to_matrix(), c1);
    }

    #[test]
    fn test_extract_sub_matrix(){
        let m :  MatrixI64 = from_range_rw_i64(20, 20, -100, 400);
        let v1   = m.view(2, 2, 6, 6);
        let m2  = v1.sub_matrix(2, 2, 4, 4);
        let m3 = m.sub_matrix(4,4, 4, 4);
        println!("m2 : {}", m2);
        assert_eq!(m2, m3);
    }


    #[test]
    fn test_col_switch(){
        let m1 = from_range_rw_i32(10, 20, 0, 500);
        println!("m1: {}", m1);
        let m2 = m1.transpose();
        let mut v1 = m1.view(2, 2, 4, 6);
        let mut v2 = m2.view(2, 2, 6, 4);
        v1.eco_switch(1, 2);
        v2.ero_switch(1, 2);
        println!("m1: {}", m1);
        println!("m2: {}", m2);
        assert_eq!(m1, m2.transpose());
    }


    #[test]
    fn test_col_scale(){
        let m1 = from_range_rw_i32(10, 20, 0, 500);
        println!("m1: {}", m1);
        let m2 = m1.transpose();
        let mut v1 = m1.view(2, 2, 4, 6);
        let mut v2 = m2.view(2, 2, 6, 4);
        v1.eco_scale(1, 2);
        v2.ero_scale(1, 2);
        println!("m1: {}", m1);
        println!("m2: {}", m2);
        assert_eq!(m1, m2.transpose());
    }

    #[test]
    fn test_col_scale_add_0(){
        let m1 = from_range_rw_i32(3, 3, 0, 100);
        println!("m1: {}", m1);
        let m2 = m1.transpose();
        println!("m2: {}", m2);
        let mut v1 = m1.view(1, 1, 1, 2);
        let mut v2 = m2.view(1, 1, 2, 1);
        println!("v1: {}", v1);
        println!("v2: {}", v2);
        v1.eco_scale_add(1, 0, 3);
        v2.ero_scale_add(1, 0, 3);
        println!("v1: {}", v1);
        println!("v2: {}", v2);
        println!("m1: {}", m1);
        println!("m2: {}", m2);
        assert_eq!(m1, m2.transpose());
    }


    #[test]
    fn test_col_scale_add_1(){
        let m1 = from_range_rw_i32(10, 20, 0, 500);
        println!("m1: {}", m1);
        let m2 = m1.transpose();
        let mut v1 = m1.view(2, 2, 4, 6);
        let mut v2 = m2.view(2, 2, 6, 4);
        v1.eco_scale_add(1, 2, 3);
        v2.ero_scale_add(1, 2, 3);
        println!("m1: {}", m1);
        println!("m2: {}", m2);
        assert_eq!(m1, m2.transpose());
    }

    #[test]
    fn test_trace(){
        let m = matrix_cw_f64(3, 3, &[1., 0., 0., 
            4., 5., 0.,
            6., 2., 3.]);
        let v = m.view(0, 0, 2, 2);
        assert_eq!(v.trace(), 6.);
        let v = m.view(1, 1, 2, 2);
        println!("{}", v);
        assert_eq!(v.trace(), 8.);
    }
}

