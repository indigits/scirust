#![doc="Matrix inverse computation methods
"]



// std imports

// local imports
use matrix::matrix::{Matrix, MatrixF64};
use matrix::traits::{Shape, Search};
use matrix::eo::eo_traits::{ERO, ECO};
use error::SRError;


/// Computes the inverse of a matrix using elementary row operations
pub fn inverse_ero(a : &mut MatrixF64) ->  Result<MatrixF64, SRError>{
    if !a.is_square(){
        return Err(SRError::IsNotSquareMatrix);
    }
    let n = a.num_rows();
    let mut result  : MatrixF64 = Matrix::identity(n, n);
    // forward elimination
    for k in 0..n{
        let (_, rr) = a.max_abs_scalar_in_col(k, k, n);
        if rr > k {
            // TODO : we can switch only part of row
            a.ero_switch(k, rr);
            result.ero_switch(k, rr);
        }
        let mut v = a.view(k, k, n - k, n - k);
        // Pick the pivot
        let pivot  = unsafe { v.get_unchecked(0, 0) } ;
        if pivot == 0. {
            return Err(SRError::IsSingular);
        }
        // bring 1 in the diagonal 
        v.ero_scale(0, 1./pivot);
        result.ero_scale(k, 1./pivot);
        for r in (1..v.num_rows()){
            let first = unsafe { v.get_unchecked(r, 0) };
            v.ero_scale_add(r, 0, -first);
            // TODO: ignore 0 entries in k-th row of result
            result.ero_scale_add(k + r, k as isize, -first);
        }
        //println!("a: {}", a);
        //println!("b: {}", result);
    }
    //println!("a: {}", a);
    //println!("b: {}", result);
    // back substitution
    let mut k = n -1;
    loop {
        // We are using (k, k) entry in a and
        // updating k-th column in a.
        for r in 0..k{
            let factor = unsafe { a.get_unchecked(r, k) };
            result.ero_scale_add(r, k as isize, -factor);
        }
        if k == 0 {
            break;
        }
        k -= 1;
    }
    //println!("a: {}", a);
    //println!("result: {}", result);
    Ok(result)
}


/// Computes the inverse of a matrix using elementary column operations
pub fn inverse_eco(a : &mut MatrixF64) ->  Result<MatrixF64, SRError>{
    if !a.is_square(){
        return Err(SRError::IsNotSquareMatrix);
    }
    let n = a.num_rows();
    let mut result  : MatrixF64 = Matrix::identity(n, n);
    // forward elimination
    for k in (0..n){
        let (_, cc) = a.max_abs_scalar_in_row(k, k, n);
        if cc > k {
            // TODO : we can switch only part of column
            a.eco_switch(k, cc);
            result.eco_switch(k, cc);
        }
        let mut v = a.view(k, k, n - k, n - k);
        // Pick the pivot
        let pivot  = unsafe { v.get_unchecked(0, 0) };
        if pivot == 0. {
            return Err(SRError::IsSingular);
        }
        // bring 1 in the diagonal 
        v.eco_scale(0, 1./pivot);
        result.eco_scale(k, 1./pivot);
        for c in (1..v.num_cols()){
            let first = unsafe { v.get_unchecked(0, c) };
            v.eco_scale_add(c, 0, -first);
            // TODO: ignore 0 entries in k-th col of result
            result.eco_scale_add(k + c, k as isize, -first);
        }
        //println!("a: {}", a);
        //println!("b: {}", result);
    }
    //println!("a: {}", a);
    //println!("b: {}", result);
    // back substitution
    let mut k = n -1;
    loop {
        // We are using (k, k) entry in a and
        // updating k-th row in a.
        for c in (0..k){
            let factor = unsafe { a.get_unchecked(k, c) };
            result.eco_scale_add(c, k as isize, -factor);
        }
        if k == 0 {
            break;
        }
        k -= 1;
    }
    //println!("a: {}", a);
    //println!("result: {}", result);
    Ok(result)
}


/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/

#[cfg(test)]
mod test{
    use num;
    use super::*;
    use matrix::matrix::*;
    use matrix::constructors::*;
    use matrix::traits::*;

    #[test]
    fn test_inv_ero_0(){
        let a = matrix_rw_f64(2, 2, &[
            1., 0.,
            1., 1.]);
        let b = inverse_ero(&mut a.clone()).unwrap();
        let c = &a * &b;
        assert!(c.is_identity());

    }

    #[test]
    fn test_inv_eco_0(){
        let a = matrix_rw_f64(2, 2, &[
            1., 0.,
            1., 1.]);
        let b = inverse_eco(&mut a.clone()).unwrap();
        let c = &a * &b;
        assert!(c.is_identity());

    }

    #[test]
    fn test_inv_ero_hadamard(){
        for i in 2..6{
            let n = num::pow(2, i);
            let a = hadamard(n).unwrap();
            println!("n: {}", n);
            println!("a: {}", a);
            let b = inverse_ero(&mut a.clone()).unwrap();
            let c = &a * &b;
            println!("b: {}", b);
            println!("c: {}", c);

            assert!(c.is_identity());
        }

    }

    #[test]
    fn test_inv_eco_hadamard(){
        for i in 2..6{
            let n = num::pow(2, i);
            let a = hadamard(n).unwrap();
            println!("n: {}", n);
            println!("a: {}", a);
            let b = inverse_eco(&mut a.clone()).unwrap();
            let c = &a * &b;
            println!("b: {}", b);
            println!("c: {}", c);

            assert!(c.is_identity());
        }
    }


    #[test]
    fn test_inv_ero_hilbert(){
        for n in 2..10{
            let a = hilbert(n);
            println!("n: {}", n);
            println!("a: {}", a);
            let b = inverse_ero(&mut a.clone()).unwrap();
            let c = &a * &b;
            println!("b: {}", b);
            println!("c: {}", c);
            let i : MatrixF64  = Matrix::identity(n, n);
            /*
            Hilbert matrices are badly conditioned. Hence,
            the numerical accuracy fails. We don't really
            get a proper identity matrix. What we get
            is a matrix which is close to identity matrix.

            A suitable way to verify correctness is to 
            measure the maximum deviation of any entry
            isize the difference matrix below from identity. 
            */
            let diff = &i  - &c;
            let max = diff.max_abs_scalar_value();
            assert!(max < 1e-3);
        }


    }


    #[test]
    fn test_inv_eco_hilbert(){
        for n in 2..10{
            let a = hilbert(n);
            println!("n: {}", n);
            println!("a: {}", a);
            let b = inverse_eco(&mut a.clone()).unwrap();
            let c = &a * &b;
            println!("b: {}", b);
            println!("c: {}", c);
            let i : MatrixF64  = Matrix::identity(n, n);
            let diff = &i  - &c;
            let max = diff.max_abs_scalar_value();
            assert!(max < 1e-3);
        }

    }
}


/******************************************************
 *
 *   Bench marks follow.
 *
 *******************************************************/

#[cfg(test)]
mod bench{
    extern crate test;
    use self::test::Bencher;
    use super::*;
    use matrix::traits::*;
    use matrix::constructors::*;

    #[bench]
    fn bench_inverse_ero_hadamard_32 (b: &mut Bencher){
        let mut a = hadamard(32).unwrap();
        b.iter(|| {
            let result = inverse_ero(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }

    #[bench]
    fn bench_inverse_ero_hadamard_64 (b: &mut Bencher){
        let mut a = hadamard(64).unwrap();
        b.iter(|| {
            let result = inverse_ero(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }    

    #[bench]
    fn bench_inverse_ero_hadamard_128 (b: &mut Bencher){
        let mut a = hadamard(128).unwrap();
        b.iter(|| {
            let result = inverse_ero(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }    

    #[bench]
    #[ignore]
    fn bench_inverse_ero_hadamard_256 (b: &mut Bencher){
        let mut a = hadamard(256).unwrap();
        b.iter(|| {
            let result = inverse_ero(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }    

    #[bench]
    fn bench_inverse_eco_hadamard_32 (b: &mut Bencher){
        let mut a = hadamard(32).unwrap();
        b.iter(|| {
            let result = inverse_eco(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }

    #[bench]
    fn bench_inverse_eco_hadamard_64 (b: &mut Bencher){
        let mut a = hadamard(64).unwrap();
        b.iter(|| {
            let result = inverse_eco(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }    

    #[bench]
    fn bench_inverse_eco_hadamard_128 (b: &mut Bencher){
        let mut a = hadamard(128).unwrap();
        b.iter(|| {
            let result = inverse_eco(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }    

    #[bench]
    fn bench_inverse_eco_hadamard_256 (b: &mut Bencher){
        let mut a = hadamard(256).unwrap();
        b.iter(|| {
            let result = inverse_eco(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }    

    #[bench]
    fn bench_inverse_eco_hadamard_512 (b: &mut Bencher){
        let mut a = hadamard(512).unwrap();
        b.iter(|| {
            let result = inverse_eco(&mut a);
            match result {
                Ok(c) => assert!(c.is_square()),
                Err(e) => panic!(e.to_string()),
            }
            
        });
    }
}