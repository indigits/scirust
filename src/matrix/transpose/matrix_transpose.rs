


// local imports
use number::Number;
use matrix::traits::{Shape, MatrixBuffer};
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
        let pa = self.as_ptr();
        let pc = result.as_mut_ptr();
        for r in range(0, rows){
            for c in range(0, cols){
                let src_offset = self.cell_to_offset(r, c);
                let dst_offset = result.cell_to_offset(c, r);
                unsafe {
                    *pc.offset(dst_offset) = *pa.offset(src_offset);
                }
            }
        }
        result
    }

    fn gram(&self) -> Matrix <T>{
        let b = self.transpose();
        b * *self
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



}



