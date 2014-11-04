// std imports


// local imports
use matrix::element::{MatrixElt};


/******************************************************
 *
 *   Iterator support of Matrix
 *
 *******************************************************/

/// An iterator over the elements of a matrix in a row
pub struct RowIterator<T:MatrixElt>{
    cols : uint,
    stride: uint,
    pos  : uint, 
    ptr : *const T
}



impl <T:MatrixElt> RowIterator<T> {
    /// Creates a new iterator object
    pub fn new(cols: uint, stride : uint, ptr : *const T)-> RowIterator<T>{
        RowIterator{cols : cols, stride : stride, pos : 0, ptr : ptr}
    }
}

impl <T:MatrixElt> Iterator<T> for RowIterator<T> {
    fn next(&mut self) -> Option<T> {
        if self.cols == self.pos{
            // No more data
            return None;
        }
        let offset = self.pos * self.stride;
        self.pos += 1;
        Some(unsafe{*self.ptr.offset(offset as int)})
    }
}

 /// An iterator over the elements of a matrix in a column
pub struct ColIterator<T:MatrixElt>{
    rows : uint,
    pos  : uint, 
    ptr : *const T
}


impl <T:MatrixElt> ColIterator<T> {
    /// Creates a new iterator object
    pub fn new(rows: uint, ptr : *const T)-> ColIterator<T>{
        ColIterator{rows : rows, pos : 0, ptr : ptr}
    }
}



impl <T:MatrixElt> Iterator<T> for ColIterator<T> {
    fn next(&mut self) -> Option<T> {
        if self.rows == self.pos{
            // No more data
            return None;
        }
        let offset = self.pos;
        self.pos += 1;
        Some(unsafe{*self.ptr.offset(offset as int)})
    }
}

/// A column major iterator over the elements of a matrix
pub struct CellIterator<T:MatrixElt>{
    rows : uint,
    cols : uint,
    stride: uint, 
    r  : uint, 
    c : uint,
    ptr : *const T
}

impl <T:MatrixElt> CellIterator<T> {
    /// Creates a new iterator object
    pub fn new(rows: uint, cols: uint, stride : uint, ptr : *const T)-> CellIterator<T>{
        CellIterator{rows : rows, 
            cols : cols, 
            stride : stride, 
            r : 0, c : 0, 
            ptr : ptr}
    }
}


impl <T:MatrixElt> Iterator<T> for CellIterator<T> {
    fn next(&mut self) -> Option<T> {
        if self.cols == self.c{
            // No more data
            return None;
        }
        let offset = self.c * self.stride + self.r;
        self.r += 1;
        if self.r == self.rows {
            self.r = 0;
            self.c += 1;
        }
        Some(unsafe{*self.ptr.offset(offset as int)})
    }
}


