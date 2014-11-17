
// Std imports
use std::num::{Float};

// local imports
use entry::One;
use matrix::*;

#[doc="Returns the determinant of a matrix.

# Remarks 

This is a naive implementation of determinant based on the
definition of determinant.  See det_float for better implementation
for floating point matrices.

Usually determinants based on
elimination or factorization are much faster.

The determinant is defined only for square matrices.

The determinant of an empty matrix is 1.  
See http://en.wikipedia.org/wiki/Matrix_(mathematics)#Empty_matrices.

"]
pub fn  det<T:Signed>(m : &Matrix<T>)->Result<T,SRError>{
    if !m.is_square(){
        return Err(IsNotSquareMatrix);
    }
    if m.is_empty(){
        return Ok(One::one());
    }
    Ok(det_naive(m))
}

#[doc="Returns the determinant of a matrix of floating point numbers.

# Remarks
"]
pub fn  det_float<T:Signed+Float>(m : &Matrix<T>)->Result<T,SRError>{
    if !m.is_square(){
        return Err(IsNotSquareMatrix);
    }
    if m.is_empty(){
        return Ok(One::one());
    }
    Ok(det_ge(&mut m.clone()))
}

/// Private implementation of determinant
/// Assumes that matrix is indeed square.
pub fn det_naive<T:Signed>(m : &Matrix<T>)->T{
    let a0 = m.get(0, 0);
    //debug!("m: {}", m);
    if m.is_scalar(){
        return a0;
    }
    let n = m.num_cols();
    let mut m2 = m.sub_matrix(1,1, n-1, n-1);
    let ps = m.as_ptr();
    let pd = m2.as_mut_ptr();
    let mut sign : T  = One::one();
    let mut result = sign * a0 * det_naive(&m2);
    sign = -sign;
    for c in range(0, n-1){
        for r in range(1, n){
            debug!("r : {}, c : {}", r , c);
            let src_offset = m.cell_to_offset(r, c);
            let dst_offset = m2.cell_to_offset(r - 1, c);
            //debug_assert!(src_offset < m.capacity() as int);
            //debug_assert!(dst_offset < m2.capacity() as int);
            unsafe {
                let v = *ps.offset(src_offset);
                //debug!("v = {}", v);
                *pd.offset(dst_offset) = v;
            }
        }
        let ai = m.get(0, c+1);
        let ai_minor_det =  det_naive(&m2);
        //debug!("sign: {}, ai: {}, Ai : {}", sign, ai, ai_minor_det);
        result = result + sign * ai * ai_minor_det;
        sign = -sign;
    }
    result
}


#[doc="Computes determinant using Gaussian 
elimination with partial pivoting
"]
pub fn det_ge<T:Number+Signed+Float>(a : &mut Matrix<T>)->T{
    assert!(a.is_square());
    let o = One::one();
    let mut result : T = o;
    let n = a.num_cols();
    // Iterate over rows
    for k in range(0, n){
        // We are working on k-th row.
        // Find the pivot position in the row
        let (_, cc) = a.max_abs_scalar_in_row(k, k, n);
        if cc > k {
            // We need to exchange columns of the submatrix.
            let mut l_tr = a.view(k, k, n - k, n - k);
            l_tr.eco_switch(0, cc - k);
            // The sign of determinant would change
            // depending on whether the permutation is
            // even or odd.
            let diff = cc - k;
            if (diff & 1) != 0 {
                // the gap in columns is odd.
                // we should change the sign of determinant.
                result = - result;
            }
        }
        // The top right part of L matrix
        let mut l_tr  = a.view(k, k, n - k, n -k);
        // Pick up the pivot
        let pivot = l_tr.get(0, 0);
        if pivot.is_zero() {
            // This is a singular matrix
            return pivot;
        }
        // Update determinant
        result = result * pivot;
        // bring 1 in the diagonal 
        l_tr.eco_scale(0, o/pivot);
        for c in range(1, l_tr.num_cols()){
            let first = l_tr.get(0, c);
            l_tr.eco_scale_add(c, 0, -first);
        }
    }
    result
}


/******************************************************
 *
 *   Unit tests 
 *
 ******************************************************/

#[cfg(test)]
mod test{
    use super::*;
    use matrix::*;
    use std::num;
    use testdata;

    #[test]
    fn test_det_0(){
        let m = matrix_rw_f64(2,2, [
            1., 2., 
            3., 4.]);
        let d = det_naive(&m);
        assert_eq!(d, -2.);
        let d = det_ge(&mut m.clone());
        assert_eq!(d, -2.);
    }


    #[test]
    fn test_det_1(){
        let m = matrix_rw_f64(3, 3, [1., 2., 3., 
            4., 5., 6., 
            7., 8., 9.]);
        let d = det_naive(&m);
        assert_eq!(d, 0.);
        let d = det_ge(&mut m.clone());
        assert_eq!(d, 0.);
    }

    #[test]
    fn test_det_hadamard(){
        let m = hadamard(4).unwrap();
        let d = det_naive(&m);
        assert_eq!(d, 16.0);
        let d = det_ge(&mut m.clone());
        assert_eq!(d, 16.0);

        let m = hadamard(8).unwrap();
        let d = det_naive(&m);
        assert_eq!(d, 4096.0);
        let d = det_ge(&mut m.clone());
        assert_eq!(d, 4096.0);
    }

    #[test]
    fn test_det_hilbert(){
        let sizes = vec![1u, 2, 4, 8];
        let determinants = vec![1.0, 0.083333333333333, 
            1.653439153439264e-07, 2.737050310006999e-33];
        let threshold =  1e-10;
        for (size, expected_value) in sizes.iter().zip(determinants.iter()){
            let m = hilbert(*size);
            let d = det_naive(&m);
            assert!((d - *expected_value).abs() < threshold);
            let d = det_ge(&mut m.clone());
            assert!((d - *expected_value).abs() < threshold);
        }
    }

    #[test]
    fn test_empty_mat_det(){
        let m : MatrixI64 = Matrix::new(0, 0);
        let d = m.det().unwrap();
        assert_eq!(d, 1);
    }

    #[test]
    fn test_examples_from_testdata(){
        assert_eq!(testdata::matrix::square_0().det().unwrap(), -13.);
        assert_eq!(testdata::matrix::square_1().det().unwrap(), 6.);
    }

}

/******************************************************
 *
 *   Benchmarks
 *
 ******************************************************/

#[cfg(test)]
mod bench{

    extern crate test;
    use self::test::Bencher;
    use super::*;
    use matrix::*;

    #[bench]
    fn bench_det_naive_hadamard_8 (b: &mut Bencher){
        let a = hadamard(8).unwrap();
        b.iter(|| {
            det_naive(&a);
        });
    }
    #[bench]
    fn bench_det_ge_hadamard_8 (b: &mut Bencher){
        let a = hadamard(8).unwrap();
        b.iter(|| {
            det_ge(&mut a.clone());
        });
    }
    #[bench]
    #[ignore]
    fn bench_det_naive_hadamard_16 (b: &mut Bencher){
        let a = hadamard(16).unwrap();
        b.iter(|| {
            det_naive(&a);
        });
    }
}
