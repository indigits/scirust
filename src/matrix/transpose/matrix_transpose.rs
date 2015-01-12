#![doc="Matrix transpose and multiplication implementations
"]

// std imports
use std::iter::range_step;
use std::cmp::min;
use std::ops;

// local imports
use number::{Number, Zero};
use matrix::traits::{Shape, MatrixBuffer, Strided};
use matrix::matrix::Matrix;
use matrix::transpose::traits::Transpose;
use error::{SRError, SRResult};


/// Checks if the matrices are transpose of each other
pub fn are_transpose<T:Number>(lhs: & Matrix<T>, rhs: & Matrix<T>) -> bool{
    if lhs.num_rows() != rhs.num_cols(){
        return false;
    }
    if lhs.num_cols() != rhs.num_rows(){
        return false;
    }
    for c in range(0, lhs.num_cols()){
        for r in range(0, lhs.num_rows()){
            if lhs.get(r, c) != rhs.get(c, r){
                return false;
            }
        }
    }
    // They are transpose of each other
    true
}

/// Simple implementation of transpose operation
pub fn transpose_simple<T:Number>(src: & Matrix<T>)->Matrix <T>{
    let rows = src.num_rows();
    let cols = src.num_cols();
    let mut result : Matrix<T> = Matrix::new(cols, rows);
    let mut psrc_col = src.as_ptr();
    let mut pdst_row = result.as_mut_ptr();
    let src_stride = src.stride() as int;
    let dst_stride = result.stride() as int;
    for _ in range(0, cols){
        let mut psrc = psrc_col;
        let mut pdst = pdst_row;
        for _ in range(0, rows){
            unsafe {
                *pdst = *psrc;
                psrc = psrc.offset(1i);
                pdst = pdst.offset(dst_stride);
            }
        }
        unsafe{
            // Move to next column in source
            psrc_col = psrc_col.offset(src_stride);
            // Move to next row in destination
            pdst_row = pdst_row.offset(1i);
        }
    }
    result
}


/// Block wise transpose
pub fn transpose_block<T:Number>(src: & Matrix<T>)->Matrix <T>{
    // Choose a block size
    let block_size = 32;
    // Construct the destination
    let rows = src.num_rows();
    let cols = src.num_cols();
    let mut result : Matrix<T> = Matrix::new(cols, rows);
    let psrc = src.as_ptr();
    let pdst = result.as_mut_ptr();
    let src_stride = src.stride() as int;
    let dst_stride = result.stride() as int;


    for cc in range_step(0, cols, block_size){
        for rr in range_step(0, rows, block_size){
            // We have to transpose a block of size blk x blk
            let blk_cols = min (block_size, cols - cc); 
            let blk_rows = min (block_size, rows - rr);
            unsafe{
                // Get our pointers to the beginning of the block
                let mut psrc_col = psrc.offset(src.cell_to_offset(rr, cc));
                let mut pdst_row = pdst.offset(result.cell_to_offset(cc, rr));

                // Now run through the block
                // iterate over columns of source block
                // iterate over rows of destination block
                for _ in range(0, blk_cols){
                    let mut psrc_blk = psrc_col;
                    let mut pdst_blk = pdst_row;
                    // iterate over the rows of a column in source block
                    // iterate over the columns of a row in destination block
                    for _ in range(0, blk_rows){
                        *pdst_blk = *psrc_blk;
                        // next row in source
                        psrc_blk = psrc_blk.offset(1i);
                        // next column in destination
                        pdst_blk = pdst_blk.offset(dst_stride);
                    }
                    // Move to next column in source block
                    psrc_col = psrc_col.offset(src_stride);
                    // Move to next row in destination block
                    pdst_row = pdst_row.offset(1i);
                }
            }
        }
    } 
    result
}



impl <T:Number> Transpose<T> for Matrix<T> {
    /// Computes the transpose of a matrix.
    /// This doesn't involve complex conjugation.
    /// Returns a new matrix
    #[inline]
    fn transpose(&self) -> Matrix <T>{
        transpose_block(self)
    }

    fn gram(&self) -> Matrix <T>{
        // simple implementation
        //let b = self.transpose();
        //b * *self
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(cols, cols);
        let ps = self.as_ptr();
        let stride = self.stride() as int;
        let z : T = Zero::zero();
        // We take advantage of the fact that the gram matrix
        // is symmetric. 
        // We only compute one half of it.
        for i in range(0, cols){
            for j in range(i, cols){
                unsafe{
                    let mut pi  = ps.offset((i as int)*stride);
                    let mut pj  = ps.offset((j as int)*stride);
                    let mut sum = z;
                    for _ in range(0, rows){
                        sum = sum + *pi * *pj;
                        // Move the pointer ahead
                        pi = pi.offset(1i);
                        pj = pj.offset(1i);
                    }
                    result.set(i, j, sum);
                    if i != j {
                        result.set(j, i, sum);
                    }
                }
            }
        }
        result

    }

}

/******************************************************
 *
 *   Matrix multiplication
 *
 *******************************************************/

 pub fn multiply_simple<T:Number>(lhs : &Matrix<T>, 
    rhs: &Matrix<T>)->SRResult<Matrix<T>>{
    // Validate dimensions match for multiplication
    if lhs.num_cols() != rhs.num_rows(){
        return Err(SRError::DimensionsMismatch);
    }
    let mut result : Matrix<T> = Matrix::new(lhs.num_rows(), rhs.num_cols());
    let pa = lhs.as_ptr();
    let pb = rhs.as_ptr();
    let pc = result.as_mut_ptr();
    let zero : T = Zero::zero();
    unsafe {
        for r in range(0, lhs.num_rows()){
            for c in range(0, rhs.num_cols()){
                let mut sum = zero;
                for j in range(0, lhs.num_cols()){
                    let lhs_offset = lhs.cell_to_offset(r, j);
                    let rhs_offset = rhs.cell_to_offset(j, c);
                    let term = *pa.offset(lhs_offset) * *pb.offset(rhs_offset);
                    sum = sum + term;
                }
                let dst_offset = result.cell_to_offset(r, c);
                *pc.offset(dst_offset)  = sum;
            }
        }
    }
    Ok(result)
}


#[doc="Computes A' * B. (transpose of lhs multiplied with rhs)

# Remarks

We don't have to compute the transpose directly. We note that this
essentially means taking the inner product of columns of A with
columns of B. We carry out the same. This is highly advantageous. 
Since the matrix is stored in column major order, hence multiplying
columns with each other is highly efficient from locality of data
perspective. This benefit shows up more as the matrix size increases.
"]
pub fn multiply_transpose_simple<T:Number>(lhs : &Matrix<T>, 
    rhs: &Matrix<T>)->SRResult<Matrix<T>>{
    // Validate dimensions match for multiplication
    if lhs.num_rows() != rhs.num_rows(){
        return Err(SRError::DimensionsMismatch);
    }
    let mut result : Matrix<T> = Matrix::new(lhs.num_cols(), rhs.num_cols());
    let pa = lhs.as_ptr();
    let pb = rhs.as_ptr();
    let pc = result.as_mut_ptr();
    let zero : T = Zero::zero();
    unsafe {
        for r in range(0, lhs.num_cols()){
            for c in range(0, rhs.num_cols()){
                let mut sum = zero;
                // Go to beginning of corresponding columns in lhs and rhs.
                let mut pl = pa.offset(lhs.cell_to_offset(0, r));
                let mut pr = pb.offset(rhs.cell_to_offset(0, c));
                for _ in range(0, lhs.num_rows()){
                    let term = *pl * *pr;
                    sum = sum + term;
                    // Move to next entry in column
                    pl = pl.offset(1i);
                    pr = pr.offset(1i);
                }
                // Write the inner product to destination
                let dst_offset = result.cell_to_offset(r, c);
                *pc.offset(dst_offset)  = sum;
            }
        }
    }
    Ok(result)
}


/// Matrix multiplication support
impl<T:Number> ops::Mul<Matrix<T>> for Matrix<T>{
    fn mul(&self, rhs: &Matrix<T>) -> Matrix<T> {
        let result = multiply_simple(self, rhs);
        match result {
            Ok(m) => m,
            Err(e) => panic!(e.to_string())
        }
    }
}


/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/


#[cfg(test)]
mod test{

    use std::iter::range_step;
    use matrix::matrix::*;
    use matrix::traits::*;
    use matrix::constructors::*;
    use super::*;

    #[test]
    fn test_transpose(){
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 2,  range(0, 4));
        assert_eq!(m.to_std_vec(), vec![0, 1, 2, 3]);
        assert_eq!(m.transpose().to_std_vec(), vec![0, 2, 1, 3]);
        assert_eq!(m.transpose().transpose().to_std_vec(), vec![0, 1, 2, 3]);
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 3,  range(0, 10));
        assert_eq!(m.transpose().to_std_vec(), vec![
            0, 2, 4, 1, 3, 5]);
        let m4 :  MatrixI64 = Matrix::from_iter_cw(2, 5, range(9, 100));
        let m5 = m4.transpose();
        println!("m4: {}", m4);
        println!("m5: {}", m5);
        assert_eq!(m5.transpose(), m4);
    }

    #[test]
    fn test_transpoe_simple_hilbert(){
        for n in range_step(100, 1024, 81){
            let h = hilbert(n);
            assert!(are_transpose(&h, &h.transpose()));
        }
    }

    #[test]
    fn test_transpoe_block_hilbert(){
        for n in range_step(100, 1024, 81){
            let h = hilbert(n);
            assert!(are_transpose(&h, & transpose_block(&h)));
        }
    }

    #[test]
    fn test_gram(){
        let m = matrix_rw_f32(4, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            4., 5., 6.,
            7., 8., 9.]);
        let g  = m.gram();
        assert_eq!(g, m.transpose() * m);
    }

    #[test]
    fn test_mult_1(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(2, 2, range(0, 4));
        let m2 : MatrixI64 = Matrix::from_iter_cw(2, 2, range(0, 4));
        let m3 = m1 * m2;
        let v = vec![2i64, 3, 6, 11];
        assert_eq!(m3.to_std_vec(), v);
        let m3 = multiply_simple(&m1, &m2).unwrap();
        let m4 = multiply_transpose_simple(&m1.transpose(), &m2).unwrap();
        assert_eq!(m3, m4);
    }


    #[test]
    fn test_mult_4(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(10, 20, range(0, 400));
        let m2 : MatrixI64 = Matrix::from_iter_cw(20, 5, range(0, 400));
        let m3 = multiply_simple(&m1, &m2).unwrap();
        let m4 = multiply_transpose_simple(&m1.transpose(), &m2).unwrap();
        assert_eq!(m3, m4);
    }

}



/******************************************************
 *
 *   Bench marks.
 *
 *******************************************************/


#[cfg(test)]
mod bench {
    extern crate test;
    use self::test::Bencher;
    use matrix::constructors::*;
    use matrix::traits::*;
    use super::*;

    #[bench]
    fn bench_transpose_simple(b: &mut Bencher){
        let a = hadamard(4096).unwrap();
        b.iter(|| {
                    transpose_simple(&a);
                });
    }

    #[bench]
    fn bench_transpose_block(b: &mut Bencher){
        let a = hadamard(4096).unwrap();
        b.iter(|| {
                    transpose_block(&a);
                });
    }


    #[bench]
    fn bench_gram_simple(b: &mut Bencher){
        let m = hadamard(128).unwrap();
        let m = m.repeat_matrix(1, 2);
        b.iter(|| {
                    // computation of gram matrix
                    m.transpose() * m;
                });

    }

    #[bench]
    fn bench_gram_optimized(b: &mut Bencher){
        let m = hadamard(128).unwrap();
        let m = m.repeat_matrix(1, 2);
        b.iter(|| {
                    // computation of gram matrix
                    m.gram();
                });

    }

    static MULTIPLY_MATRIX_SIZE : uint = 512;

    #[bench]
    fn bench_transpose_block_mult_mat_size(b: &mut Bencher){
        let a = hadamard(MULTIPLY_MATRIX_SIZE).unwrap();
        b.iter(|| {
                    transpose_block(&a);
                });
    }

    #[bench]
    fn bench_multiply_simple(b: &mut Bencher){
        let m = hadamard(MULTIPLY_MATRIX_SIZE).unwrap();
        b.iter(|| {
                    // Matrix multiplication
                    multiply_simple(&m, &m).is_ok();
                });
    }

    #[bench]
    fn bench_multiply_transpose_simple(b: &mut Bencher){
        let m = hadamard(MULTIPLY_MATRIX_SIZE).unwrap();
        b.iter(|| {
                    // take transpose
                    let m2 = m.transpose();
                    // Matrix multiplication with transpose
                    multiply_transpose_simple(&m2, &m).is_ok();
                });
    }
}
