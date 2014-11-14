// std imports
use std::num;

// local imports
use error::*;
use matrix::view::MatrixView;
use matrix::traits::*;
use number::Number;


/// Implementation of min-max with absolute values API for matrix view
impl <'a, T:Number+PartialOrd+Signed> MinMaxAbs<T> for MatrixView<'a, T> {

    // Returns the absolute minimum scalar value
    fn min_abs_scalar(&self) -> (T, uint, uint){
        if self.is_empty(){
            fail!(EmptyMatrix.to_string());
        }
        let mut v = num::abs(self.get(0, 0));
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.matrix().as_ptr();
        for c in range(0, self.num_cols()){
            for r in range(0, self.num_rows()){
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
    fn max_abs_scalar(&self) -> (T, uint, uint){
        if self.is_empty(){
            fail!(EmptyMatrix.to_string());
        }
        let mut v = num::abs(self.get(0, 0));
        // The location
        let mut rr = 0;
        let mut cc = 0;
        let ps = self.matrix().as_ptr();
        for c in range(0, self.num_cols()){
            for r in range(0, self.num_rows()){
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
}

