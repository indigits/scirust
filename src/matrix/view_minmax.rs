// std imports

// external imports
use num::traits::Signed;


// local imports
use error::SRError;
use algebra::structure::CommutativeMonoidAddPartial;
use matrix::view::MatrixView;
use matrix::traits::*;

/// Implementation of min-max with absolute values API for matrix view
impl <'a, T:CommutativeMonoidAddPartial+Signed+PartialOrd> MinMaxAbs<T> for MatrixView<'a, T> {

    // Returns the absolute minimum scalar value
    fn min_abs_scalar(&self) -> (T, usize, usize){
        if self.is_empty(){
            panic!(SRError::EmptyMatrix.to_string());
        }
        let mut v = unsafe {self.get_unchecked(0, 0)}.abs();
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.matrix().as_ptr();
        for c in 0..self.num_cols(){
            for r in 0..self.num_rows(){
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
    fn max_abs_scalar(&self) -> (T, usize, usize){
        if self.is_empty(){
            panic!(SRError::EmptyMatrix.to_string());
        }
        let mut v = unsafe {self.get_unchecked(0, 0)}.abs();
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.matrix().as_ptr();
        for c in 0..self.num_cols(){
            for r in 0..self.num_rows(){
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
}

