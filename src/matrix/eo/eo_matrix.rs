

// std imports

// local imports
use algebra::Number;
use super::eo_traits::{ERO, ECO};
use matrix::matrix::Matrix;

/// Implementation of Elementary row operations.
impl<T:Number> ERO<T> for Matrix<T> {
}

/// Implementation of Elementary column operations.
impl<T:Number> ECO<T> for Matrix<T> {
}


/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/


#[cfg(test)]
mod test{
    use api::*;

    #[test]
    fn test_row_switch(){
        let mut m1 : MatrixI64 = Matrix::from_slice_cw(3, 3, vec![2, 3, 9, 2, 1, 7, 4, 2, 6].as_slice());
        println!("m1: {}", m1);
        let m2 = m1.clone();
        m1.ero_switch(1, 2);
        let m3 : MatrixI64 = Matrix::from_slice_cw(3, 3, vec![2, 9, 3, 2, 7, 1, 4, 6, 2].as_slice());
        assert_eq!(m1, m3);
        m1.ero_switch(2, 1);
        assert_eq!(m1, m2);
    }


    #[test]
    fn test_row_scale(){
        let mut m1 : MatrixI64 = Matrix::from_slice_cw(3, 3, vec![
            2, 3, 9, 
            2, 1, 7, 
            4, 2, 6].as_slice());
        println!("m1: {}", m1);
        m1.ero_scale(1, 2);
        let m3 : MatrixI64 = Matrix::from_slice_cw(3, 3, vec![
            2, 6, 9, 
            2, 2, 7, 
            4, 4, 6].as_slice());
        assert_eq!(m1, m3);
    }

    #[test]
    fn test_row_scale_slice_1(){
        let mut m1 = matrix_rw_i64(3,3,[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9].as_slice());
        m1.ero_scale_slice(1, 2, 1, 3);
        let m2 = matrix_rw_i64(3,3,[
            1, 2, 3,
            4, 10, 12,
            7, 8, 9].as_slice());
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_row_scale_add(){
        let mut m1 : MatrixI64 = Matrix::from_slice_cw(3, 3, vec![
            2, 3, 9, 
            2, 1, 7, 
            4, 2, 6].as_slice());
        println!("m1: {}", m1);
        m1.ero_scale_add(1, 2, 10);
        let m3 : MatrixI64 = Matrix::from_slice_cw(3, 3, vec![
            2, 93, 9, 
            2, 71, 7, 
            4, 62, 6].as_slice());
        assert_eq!(m1, m3);
    }

    #[test]
    fn test_col_switch(){
        let mut m1 = from_range_rw_i32(10, 6, 0, 500);
        println!("m1: {}", m1);
        let mut m2 = m1.transpose();
        m1.eco_switch(1, 2);
        m2.ero_switch(1, 2);
        assert_eq!(m1, m2.transpose());
    }


    #[test]
    fn test_col_scale(){
        let mut m1 = from_range_rw_i32(10, 6, 0, 500);
        println!("m1: {}", m1);
        let mut m2 = m1.transpose();
        m1.eco_scale(1, 2);
        m2.ero_scale(1, 2);
        assert_eq!(m1, m2.transpose());
    }

    #[test]
    fn test_col_scale_slice_1(){
        let mut m1 = matrix_rw_i64(3,3,[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9].as_slice());
        m1.eco_scale_slice(1, 2, 1, 3);
        let m2 = matrix_rw_i64(3,3,[
            1, 2, 3,
            4, 10, 6,
            7, 16, 9].as_slice());
        assert_eq!(m1, m2);
    }

    #[test]
    fn test_col_scale_add(){
        let mut m1 = from_range_rw_i32(10, 6, 0, 500);
        println!("m1: {}", m1);
        let mut m2 = m1.transpose();
        m1.eco_scale_add(1, 2, 3);
        m2.ero_scale_add(1, 2, 3);
        assert_eq!(m1, m2.transpose());
    }
}


/******************************************************
 *
 *   Bench marks
 *
 *******************************************************/


#[cfg(test)]
mod bench{
    extern crate test;
    use self::test::Bencher;
    use matrix::eo::eo_traits::*;
    use matrix::constructors::*;
    #[bench]
    fn bench_eo_col_switch(b: &mut Bencher){
        let mut m = from_range_rw_f64(1024, 1024, 0., 10000000.);
        b.iter(|| {
                    m.eco_switch(1, 100);
                });
    }

    #[bench]
    fn bench_eo_row_switch(b: &mut Bencher){
        let mut m = from_range_rw_f64(1024, 1024, 0., 10000000.);
        b.iter(|| {
                    m.ero_switch(1, 100);
                });
    }

    #[bench]
    fn bench_eo_col_scale(b: &mut Bencher){
        let mut m = from_range_rw_f64(1024, 1024, 0., 10000000.);
        b.iter(|| {
                    m.eco_scale(50, 10.);
                });
    }

    #[bench]
    fn bench_eo_row_scale(b: &mut Bencher){
        let mut m = from_range_rw_f64(1024, 1024, 0., 10000000.);
        b.iter(|| {
                    m.ero_scale(50, 10.);
                });
    }
}


