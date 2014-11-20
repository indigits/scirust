

/******************************************************
 *
 *   Unit tests
 *
 *******************************************************/
#[cfg(test)]
mod test{

    use super::*;
    use api::*;


    #[test]
    fn test_moment_sum_cw_2(){
        let m = matrix_rw_i32(3, 3, &[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9]);
        let v = m.view(1, 1, 2, 2);
        let s = sum_cw(&v);
        assert_eq!(s, matrix_cw_i32(1,2, &[13, 15]));
    }

    #[test]
    fn test_moment_sum_rw_2(){
        let m = matrix_rw_i32(3, 3, &[
            1, 2, 3,
            4, 5, 6,
            7, 8, 9]);
        let v = m.view(1, 1, 2, 2);
        let s = sum_rw(&v);
        assert_eq!(s, matrix_rw_i32(2,1, &[11, 17]));
    }

    #[test]
    fn test_moment_sum_sqr_cw_2(){
        let m = matrix_rw_i32(3, 3, &[
            1, 1, 2,
            2, 2, 1,
            3, 2, 2]);
        let v = m.view(1, 1, 2, 2);
        let s = sum_sqr_cw(&v);
        assert_eq!(s, matrix_cw_i32(1,2, &[8, 5]));
    }

    #[test]
    fn test_moment_sum_sqr_rw_2(){
        let m = matrix_rw_i32(3, 3, &[
            1, 1, 2,
            2, 2, 1,
            3, 2, 2]);
        let v = m.view(1, 1, 2, 2);
        let s = sum_sqr_rw(&v);
        assert_eq!(s, matrix_cw_i32(2,1, &[5, 8]));
    }

    #[test]
    fn test_moment_mean_cw_2(){
        let m = matrix_rw_f32(3, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            7., 8., 9.]);
        let v = m.view(1, 1, 2, 2);
        let s = mean_cw(&v);
        assert_eq!(s, matrix_cw_f32(1,2, &[13./2., 15./2.]));
    }

    #[test]
    fn test_moment_mean_rw_2(){
        let m = matrix_rw_f32(3, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            7., 8., 9.]);
        let v = m.view(1, 1, 2, 2);
        let s = mean_rw(&v);
        assert_eq!(s, matrix_rw_f32(2,1, &[11./2., 17./2.]));
    }

    #[test]
    fn test_moment_mean_sqr_cw_2(){
        let m = matrix_rw_f32(3, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            7., 8., 9.]);
        let v = m.view(1, 1, 2, 2);
        let s = mean_sqr_cw(&v);
        assert_eq!(s, matrix_cw_f32(1,2, &[44.5, 58.5]));
    }

    #[test]
    fn test_moment_mean_sqr_rw_2(){
        let m = matrix_rw_f32(4, 3, &[
            1., 2., 3.,
            4., 5., 6.,
            4., 5., 6.,
            7., 8., 9.]);
        let v = m.view(1, 1, 2, 2);
        let s = mean_sqr_rw(&v);
        assert_eq!(s, matrix_cw_f32(2,1, &[30.5, 30.5]));
    }


}
