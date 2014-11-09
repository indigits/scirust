
// Std imports
use std::num;

// local imports

use matrix::*;

#[doc="Returns the determinant of a matrix.

# Remarks 

The determinant is defined only for square matrices.

The determinant of an empty matrix is 1.  
See http://en.wikipedia.org/wiki/Matrix_(mathematics)#Empty_matrices.

"]
pub fn  det<T:Number+Signed>(m : &Matrix<T>)->Result<T,MatrixError>{
    if !m.is_square(){
        return Err(NonSquareMatrix);
    }
    if m.is_empty(){
        return Ok(num::One::one());
    }
    Ok(det_(m))
}


/// Private implementation of determinant
/// Assumes that matrix is indeed square.
fn det_<T:Number+Signed>(m : &Matrix<T>)->T{
    let a0 = m.get(0, 0);
    debug!("m: {}", m);
    if m.is_scalar(){
        return a0;
    }
    let n = m.num_cols();
    let mut m2 = m.sub_matrix(1,1, n-1, n-1);
    let ps = m.as_ptr();
    let pd = m2.as_mut_ptr();
    let mut sign : T  = num::One::one();
    let mut result = sign * a0 * det_(&m2);
    sign = -sign;
    for c in range(0, n-1){
        for r in range(1, n){
            debug!("r : {}, c : {}", r , c);
            let src_offset = m.cell_to_offset(r, c);
            let dst_offset = m2.cell_to_offset(r - 1, c);
            debug_assert!(src_offset < m.capacity() as int);
            debug_assert!(dst_offset < m2.capacity() as int);
            unsafe {
                let v = *ps.offset(src_offset);
                //debug!("v = {}", v);
                *pd.offset(dst_offset) = v;
            }
        }
        let ai = m.get(0, c+1);
        let ai_minor_det =  det_(&m2);
        debug!("sign: {}, ai: {}, Ai : {}", sign, ai, ai_minor_det);
        result = result + sign * ai * ai_minor_det;
        sign = -sign;
    }
    result
}

#[cfg(test)]
mod test{
    use super::*;
    use matrix::*;
    use std::num;

    #[test]
    fn test_det_0(){
        let m = matrix_cw_i64(2,2, [1, 2, 3, 4]);
        let d = det(&m).unwrap();
        assert_eq!(d, -2);
    }


    #[test]
    fn test_det_1(){
        let m = matrix_cw_i64(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]);
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

    #[test]
    fn test_empty_mat_det(){
        let m : MatrixI64 = Matrix::new(0, 0);
        let d = m.det().unwrap();
        assert_eq!(d, 1);
    }

    #[test]
    fn test_examples_from_testdata(){
        assert_eq!(testdata::square_0().det().unwrap(), -13.);
        assert_eq!(testdata::square_1().det().unwrap(), 6.);
    }

}
