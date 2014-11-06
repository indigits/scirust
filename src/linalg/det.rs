
// Std imports
use std::num::{Zero, One};

// local imports

use matrix::*;

/*
/// A special view in matrix designed for computing determinants
struct DetView<'a, T:'a+MatrixElt>{
    // The rows of the matrix which have been excluded
    pub excluded_rows :Vec<bool>,
    pub excluded_cols : Vec<bool>,
    /// Number of rows or columns in the view
    view_size : uint,
    m : &'a Matrix<T>
}

impl <'a, T:MatrixElt> DetView<'a, T> {
    /// Creates a new view
    pub fn new(m : &Matrix<T>)->DetView<T>{
        let n = m.num_rows();
        let result = DetView{m : m,  
            excluded_rows  : Vec::with_capacity(n),
            excluded_cols : Vec::with_capacity(n),
            view_size : 0
        };
        result

    }
}
*/

#[doc="Returns the determinant of a matrix view.
"]
pub fn  det<T:MatrixElt+Signed>(m : &Matrix<T>)->Result<T,MatrixError>{
    let n = m.num_cols();
    if !m.is_square(){
        return Err(NonSquareMatrix);
    }
    if m.is_scalar(){
        return Ok(m.get(0, 0));
    }
    //println!("m: {}", m);
    let mut result : T = Zero::zero();
    let mut m2 : Matrix<T> = Matrix::new(n-1, n-1);
    let ps = m.as_ptr();
    let pd = m2.as_mut_ptr();
    let mut sign : T = One::one();
    for i in range(0, n){
        for c in range(0, i){
            for r in range(1, n){
                //println!("first: i : {}, c : {}, r: {}", i , c, r);
                let src_offset = m.cell_to_offset(r, c);
                let dst_offset = m2.cell_to_offset(r  - 1, c);
                debug_assert!(src_offset < m.capacity() as int);
                debug_assert!(dst_offset < m2.capacity() as int);
                unsafe {
                    let v = *ps.offset(src_offset);
                    //println!("v = {}", v);
                    *pd.offset(dst_offset) = v;
                }
            }
        }
        for c in range(i+1, n){
            for r in range(1, n){
                //println!("second: i : {}, c : {}, r: {}", i , c, r);
                let src_offset = m.cell_to_offset(r, c);
                let dst_offset = m2.cell_to_offset(r  - 1, c - 1);
                debug_assert!(src_offset < m.capacity() as int);
                debug_assert!(dst_offset < m2.capacity() as int);
                unsafe {
                    let v = *ps.offset(src_offset);
                    //println!("v = {}", v);
                    *pd.offset(dst_offset) = v;
                }
            }
        }
        let ai = m.get(0, i);
        let ai_minor_det =  det(&m2).unwrap();
        //println!("sign: {}, ai: {}, Ai : {}", sign, ai, ai_minor_det);
        result = result + sign * ai * ai_minor_det;
        sign = -sign;
    }
    Ok(result)
}


#[cfg(test)]
mod test{
    use super::*;
    use matrix::*;
    use std::num;

    #[test]
    fn test_det_0(){
        let m = matrix_i64(2,2, [1, 2, 3, 4]);
        let d = det(&m).unwrap();
        assert_eq!(d, -2);
    }


    #[test]
    fn test_det_1(){
        let m = matrix_i64(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]);
        let d = det(&m).unwrap();
        assert_eq!(d, 0);
    }

    #[test]
    fn test_det_hadamard(){
        let m = hadamard(4).unwrap();
        let d = det(&m).unwrap();
        assert_eq!(d, 16.0);

        let m = hadamard(8).unwrap();
        let d = det(&m).unwrap();
        assert_eq!(d, 4096.0);
    }

    #[test]
    fn test_det_hilbert(){
        assert_eq!(hilbert(1).det().unwrap(), 1.0);
        assert!(num::abs(hilbert(2).det().unwrap() - 0.083333333333333) < 1e-10);
        assert!(num::abs(hilbert(4).det().unwrap() - 1.653439153439264e-07) < 1e-10);
        assert!(num::abs(hilbert(8).det().unwrap() - 2.737050310006999e-33) < 1e-10);
    }

}
