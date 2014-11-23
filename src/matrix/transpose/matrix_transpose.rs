


// local imports
use number::{Number, Zero};
use matrix::traits::{Shape, MatrixBuffer, Strided};
use matrix::matrix::Matrix;
use matrix::transpose::traits::Transpose;


impl <T:Number> Transpose<T> for Matrix<T> {
    /// Computes the transpose of a matrix.
    /// This doesn't involve complex conjugation.
    /// Returns a new matrix
    fn transpose(&self) -> Matrix <T>{
        let rows = self.num_rows();
        let cols = self.num_cols();
        let mut result : Matrix<T> = Matrix::new(cols, rows);
        let mut psrc_col = self.as_ptr();
        let mut pdst_row = result.as_mut_ptr();
        let src_stride = self.stride() as int;
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
 *   Unit tests
 *
 *******************************************************/


#[cfg(test)]
mod test{

    use matrix::matrix::*;
    use matrix::traits::*;
    use matrix::constructors::*;

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
    fn test_gram(){
        let m = matrix_rw_f32(4, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            4., 5., 6.,
            7., 8., 9.]);
        let g  = m.gram();
        assert_eq!(g, m.transpose() * m);
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

    #[bench]
    fn bench_transpose(b: &mut Bencher){
        let a = hadamard(1024).unwrap();
        b.iter(|| {
                    a.transpose();
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
}
