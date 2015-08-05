#![doc="Matrix transpose and multiplication implementations
"]

// std imports
use std::cmp::min;
use std::ops;

// external imports
use num::traits::{Zero};

// local imports
use algebra::structure::{MagmaBase, CommutativeMonoidAddPartial, CommutativeMonoidMulPartial};
use matrix::traits::{Shape, MatrixBuffer, Strided};
use matrix::matrix::Matrix;
use matrix::transpose::traits::{Transpose, Frame};
use error::{SRError, SRResult};


/// Checks if the matrices are transpose of each other
pub fn are_transpose<T:MagmaBase>(lhs: & Matrix<T>, rhs: & Matrix<T>) -> bool{
    if lhs.num_rows() != rhs.num_cols(){
        return false;
    }
    if lhs.num_cols() != rhs.num_rows(){
        return false;
    }
    for c in 0..lhs.num_cols(){
        for r in 0..lhs.num_rows(){
            if lhs.get(r, c) != rhs.get(c, r){
                return false;
            }
        }
    }
    // They are transpose of each other
    true
}

/// Simple implementation of transpose operation
pub fn transpose_simple<T:MagmaBase>(src: & Matrix<T>)->Matrix <T>{
    let rows = src.num_rows();
    let cols = src.num_cols();
    let mut result : Matrix<T> = Matrix::new(cols, rows);
    let mut psrc_col = src.as_ptr();
    let mut pdst_row = result.as_mut_ptr();
    let src_stride = src.stride() as isize;
    let dst_stride = result.stride() as isize;
    for _ in 0..cols{
        let mut psrc = psrc_col;
        let mut pdst = pdst_row;
        for _ in 0..rows{
            unsafe {
                *pdst = *psrc;
                psrc = psrc.offset(1);
                pdst = pdst.offset(dst_stride);
            }
        }
        unsafe{
            // Move to next column in source
            psrc_col = psrc_col.offset(src_stride);
            // Move to next row in destination
            pdst_row = pdst_row.offset(1);
        }
    }
    result
}


/// Block wise transpose
pub fn transpose_block<T:MagmaBase>(src: & Matrix<T>)->Matrix <T>{
    // Choose a block size
    let block_size = 32;
    // Construct the destination
    let rows = src.num_rows();
    let cols = src.num_cols();
    let mut result : Matrix<T> = Matrix::new(cols, rows);
    let psrc = src.as_ptr();
    let pdst = result.as_mut_ptr();
    let src_stride = src.stride() as isize;
    let dst_stride = result.stride() as isize;


    for cc in (0..cols).step_by(block_size){
        for rr in (0..rows).step_by(block_size){
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
                for _ in 0..blk_cols{
                    let mut psrc_blk = psrc_col;
                    let mut pdst_blk = pdst_row;
                    // iterate over the rows of a column in source block
                    // iterate over the columns of a row in destination block
                    for _ in 0..blk_rows{
                        *pdst_blk = *psrc_blk;
                        // next row in source
                        psrc_blk = psrc_blk.offset(1);
                        // next column in destination
                        pdst_blk = pdst_blk.offset(dst_stride);
                    }
                    // Move to next column in source block
                    psrc_col = psrc_col.offset(src_stride);
                    // Move to next row in destination block
                    pdst_row = pdst_row.offset(1);
                }
            }
        }
    }
    result
}



impl <T:MagmaBase> Transpose<T> for Matrix<T> {
    /// Computes the transpose of a matrix.
    /// This doesn't involve complex conjugation.
    /// Returns a new matrix
    #[inline]
    fn transpose(&self) -> Matrix <T>{
        transpose_block(self)
    }
}

impl <T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial> Frame<T> for Matrix<T> {

    fn gram(&self) -> Matrix <T>{
        // simple implementation
        //let b = self.transpose();
        //b * *self
        let cols = self.num_cols();
        let rows = self.num_rows();
        let mut result = Matrix::new(cols, cols);
        let ps = self.as_ptr();
        let stride = self.stride() as isize;
        let z : T = Zero::zero();
        // We take advantage of the fact that the gram matrix
        // is symmetric.
        // We only compute one half of it.
        for i in 0..cols{
            for j in i..cols{
                unsafe{
                    let mut pi  = ps.offset((i as isize)*stride);
                    let mut pj  = ps.offset((j as isize)*stride);
                    let mut sum = z;
                    for _ in 0..rows{
                        sum = sum + *pi * *pj;
                        // Move the pointer ahead
                        pi = pi.offset(1);
                        pj = pj.offset(1);
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


 pub fn multiply_block<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial>(lhs : &Matrix<T>,
    rhs: &Matrix<T>)->SRResult<Matrix<T>>{
    use std::{cmp, slice};
    use std::ops::Add;
    const BLOCK_SIZE:usize = 16;
    // Validate dimensions match for multiplication
    if lhs.num_cols() != rhs.num_rows(){
        return Err(SRError::DimensionsMismatch);
    }
    let mut result : Matrix<T> = Matrix::new(lhs.num_rows(), rhs.num_cols());
    let zero : T = Zero::zero();
    let pa = lhs.as_ptr();
    let pb = rhs.as_ptr();
    let pc = result.as_mut_ptr();
    for row_start in (0..lhs.num_rows()).step_by(BLOCK_SIZE) {
        for col_start in (0..rhs.num_cols()).step_by(BLOCK_SIZE) {
            let row_end = cmp::min(row_start + BLOCK_SIZE, lhs.num_rows());
            let col_end = cmp::min(col_start + BLOCK_SIZE, rhs.num_cols());

            let mut res_block = [[zero; BLOCK_SIZE]; BLOCK_SIZE];

            for k_start in (0..lhs.num_cols()).step_by(BLOCK_SIZE) {
                let k_end = cmp::min(k_start + BLOCK_SIZE, lhs.num_cols());

                let mut lhs_block = [[zero; BLOCK_SIZE]; BLOCK_SIZE];
                let mut rhs_block = [[zero; BLOCK_SIZE]; BLOCK_SIZE];

                for k in k_start..k_end {
                    let k_offs = k - k_start;
                    let col_slice = unsafe {
                        slice::from_raw_parts(pa.offset(lhs.cell_to_offset(row_start, k)), row_end - row_start)
                    };
                    for (row_offs, value) in col_slice.iter().cloned().enumerate() {
                        lhs_block[row_offs][k_offs] = value;
                    }
                }
                for col in col_start..col_end {
                    let col_offs = col - col_start;
                    let col_slice = unsafe {
                        slice::from_raw_parts(pb.offset(rhs.cell_to_offset(k_start, col)), k_end - k_start)
                    };
                    for (k_offs, value) in  col_slice.iter().cloned().enumerate() {
                        // Note the inverted indices
                        rhs_block[col_offs][k_offs] = value;
                    }
                }
                for i in 0..(row_end - row_start) {
                    for j in 0..(col_end - col_start) {
                        let cell = lhs_block[i][0..(k_end - k_start)].iter().cloned()
                            .zip(rhs_block[j][0..(k_end - k_start)].iter().cloned())
                            .map(|(l, r)| l * r)
                            .fold(zero, Add::add);
                        res_block[i][j] = res_block[i][j] + cell;
                    }
                }
            }
            for col in col_start..col_end {
                let col_slice = unsafe {
                    slice::from_raw_parts_mut(pc.offset(result.cell_to_offset(row_start, col)), row_end - row_start)
                };
                for (row_offs, res_ptr) in col_slice.iter_mut().enumerate() {
                    *res_ptr = res_block[row_offs][col - col_start];
                }
            }
        }
    }
    return Ok(result);
}

 pub fn multiply_simple<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial>(lhs : &Matrix<T>,
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
        for r in 0..lhs.num_rows(){
            for c in 0..rhs.num_cols(){
                let mut sum = zero;
                for j in 0..lhs.num_cols(){
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
pub fn multiply_transpose_simple<T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial>(lhs : &Matrix<T>,
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
        for r in 0..lhs.num_cols(){
            for c in 0..rhs.num_cols(){
                let mut sum = zero;
                // Go to beginning of corresponding columns in lhs and rhs.
                let mut pl = pa.offset(lhs.cell_to_offset(0, r));
                let mut pr = pb.offset(rhs.cell_to_offset(0, c));
                for _ in 0..lhs.num_rows(){
                    let term = *pl * *pr;
                    sum = sum + term;
                    // Move to next entry in column
                    pl = pl.offset(1);
                    pr = pr.offset(1);
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
impl<'a, 'b, T:CommutativeMonoidAddPartial+CommutativeMonoidMulPartial> ops::Mul<&'b Matrix<T>> for &'a Matrix<T>{
    type Output = Matrix<T>;
    fn mul(self, rhs: &Matrix<T>) -> Matrix<T> {
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

    use matrix::matrix::*;
    use matrix::traits::*;
    use matrix::constructors::*;
    use super::*;

    #[test]
    fn test_transpose(){
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 2,  (0..4));
        assert_eq!(m.to_std_vec(), vec![0, 1, 2, 3]);
        assert_eq!(m.transpose().to_std_vec(), vec![0, 2, 1, 3]);
        assert_eq!(m.transpose().transpose().to_std_vec(), vec![0, 1, 2, 3]);
        let m  : MatrixI64 = Matrix::from_iter_cw(2, 3,  (0..10));
        assert_eq!(m.transpose().to_std_vec(), vec![
            0, 2, 4, 1, 3, 5]);
        let m4 :  MatrixI64 = Matrix::from_iter_cw(2, 5, (9..100));
        let m5 = m4.transpose();
        println!("m4: {}", m4);
        println!("m5: {}", m5);
        assert_eq!(m5.transpose(), m4);
    }

    #[test]
    fn test_transpoe_simple_hilbert(){
        for n in (100..1024).step_by(81){
            let h = hilbert(n);
            assert!(are_transpose(&h, &h.transpose()));
        }
    }

    #[test]
    fn test_transpoe_block_hilbert(){
        for n in (100..1024).step_by(81){
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
        assert_eq!(g, &(m.transpose()) * &m);
    }

    #[test]
    fn test_mult_1(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(2, 2, (0..4));
        let m2 : MatrixI64 = Matrix::from_iter_cw(2, 2, (0..4));
        let m3 = &m1 * &m2;
        let v = vec![2i64, 3, 6, 11];
        assert_eq!(m3.to_std_vec(), v);
        let m3 = multiply_simple(&m1, &m2).unwrap();
        let m4 = multiply_transpose_simple(&m1.transpose(), &m2).unwrap();
        assert_eq!(m3, m4);
    }




    #[test]
    fn test_mult_4(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(10, 20, (0..400));
        let m2 : MatrixI64 = Matrix::from_iter_cw(20, 5, (0..400));
        let m3 = multiply_simple(&m1, &m2).unwrap();
        let m4 = multiply_transpose_simple(&m1.transpose(), &m2).unwrap();
        assert_eq!(m3, m4);
    }

    #[test]
    fn test_mult_block_1(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(2, 2, (0..4));
        let m2 : MatrixI64 = Matrix::from_iter_cw(2, 2, (0..4));
        let m3 = &m1 * &m2;
        let v = vec![2i64, 3, 6, 11];
        assert_eq!(m3.to_std_vec(), v);
        let m3 = multiply_block(&m1, &m2).unwrap();
        let m4 = multiply_transpose_simple(&m1.transpose(), &m2).unwrap();
        assert_eq!(m3, m4);
    }

    #[test]
    fn test_mult_block_4(){
        let m1 : MatrixI64 = Matrix::from_iter_cw(10, 20, (0..400));
        let m2 : MatrixI64 = Matrix::from_iter_cw(20, 5, (0..400));
        let m3 = multiply_block(&m1, &m2).unwrap();
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
                    &(m.transpose()) * &m;
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

    static MULTIPLY_MATRIX_SIZE : usize = 512;

    #[bench]
    fn bench_transpose_block_mult_mat_size(b: &mut Bencher){
        let a = hadamard(MULTIPLY_MATRIX_SIZE).unwrap();
        b.iter(|| {
                    transpose_block(&a);
                });
    }

    macro_rules! make_multiply_bench {
        ($simple_name:ident, $block_name:ident, $matrix_size:expr) => (
            #[bench]
            fn $simple_name(b: &mut Bencher){
                let m = hadamard($matrix_size).unwrap();
                b.iter(|| {
                            // Matrix multiplication
                            multiply_simple(&m, &m).is_ok();
                        });
            }
            #[bench]
            fn $block_name(b: &mut Bencher){
                let m = hadamard($matrix_size).unwrap();
                b.iter(|| {
                            // Matrix multiplication
                            multiply_block(&m, &m).is_ok();
                        });
            }
        );
        // ($name_end:ident, $matrix_size:expr) => (
        //     make_multiply_bench(concat_idents!(bench_multiply_simple, $name_end),
        //         concat_idents!(bench_multiply_block, $name_end),
        //         $matrix_size);
        //     )
    }
    //the position of the numbers and the leading zeroes sort the benches intuitively
    make_multiply_bench!(bench_multiply_0001_simple, bench_multiply_0001_block, 1);
    make_multiply_bench!(bench_multiply_0002_simple, bench_multiply_0002_block, 2);
    make_multiply_bench!(bench_multiply_0004_simple, bench_multiply_0004_block, 4);
    make_multiply_bench!(bench_multiply_0008_simple, bench_multiply_0008_block, 8);
    make_multiply_bench!(bench_multiply_0016_simple, bench_multiply_0016_block, 16);
    make_multiply_bench!(bench_multiply_0032_simple, bench_multiply_0032_block, 32);
    make_multiply_bench!(bench_multiply_0064_simple, bench_multiply_0064_block, 64);
    make_multiply_bench!(bench_multiply_0128_simple, bench_multiply_0128_block, 128);
    make_multiply_bench!(bench_multiply_0256_simple, bench_multiply_0256_block, 256);
    make_multiply_bench!(bench_multiply_0512_simple, bench_multiply_0512_block, 512);
    make_multiply_bench!(bench_multiply_1024_simple, bench_multiply_1024_block, 1024);

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
