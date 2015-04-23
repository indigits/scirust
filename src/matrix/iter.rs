// std imports

// external imports


// local imports
use algebra::structure::{MagmaBase, FieldPartial};


/******************************************************
 *
 *   Iterator support of Matrix
 *
 *******************************************************/

/// An iterator over the elements of a matrix in a row
pub struct RowIterator<T:MagmaBase>{
    cols : usize,
    stride: usize,
    pos  : usize, 
    ptr : *const T
}



impl <T:MagmaBase> RowIterator<T> {
    /// Creates a new iterator object
    pub fn new(cols: usize, stride : usize, ptr : *const T)-> RowIterator<T>{
        RowIterator{cols : cols, stride : stride, pos : 0, ptr : ptr}
    }
}

impl <T:MagmaBase> Iterator for RowIterator<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        if self.cols == self.pos{
            // No more data
            return None;
        }
        let offset = self.pos * self.stride;
        self.pos += 1;
        Some(unsafe{*self.ptr.offset(offset as isize)})
    }
}

 /// An iterator over the elements of a matrix in a column
pub struct ColIterator<T:MagmaBase>{
    rows : usize,
    pos  : usize, 
    ptr : *const T
}


impl <T:MagmaBase> ColIterator<T> {
    /// Creates a new iterator object
    pub fn new(rows: usize, ptr : *const T)-> ColIterator<T>{
        ColIterator{rows : rows, pos : 0, ptr : ptr}
    }
}



impl <T:MagmaBase> Iterator for ColIterator<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        if self.rows == self.pos{
            // No more data
            return None;
        }
        let offset = self.pos;
        self.pos += 1;
        Some(unsafe{*self.ptr.offset(offset as isize)})
    }
}

/// A column major iterator over the elements of a matrix
pub struct CellIterator<T:MagmaBase>{
    rows : usize,
    cols : usize,
    stride: usize, 
    r  : usize, 
    c : usize,
    ptr : *const T
}

impl <T:MagmaBase> CellIterator<T> {
    /// Creates a new iterator object
    pub fn new(rows: usize, cols: usize, stride : usize, ptr : *const T)-> CellIterator<T>{
        CellIterator{rows : rows, 
            cols : cols, 
            stride : stride, 
            r : 0, c : 0, 
            ptr : ptr}
    }
}


impl <T:MagmaBase> Iterator for CellIterator<T> {
    type Item = T;
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
        Some(unsafe{*self.ptr.offset(offset as isize)})
    }
}


