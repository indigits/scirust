#![doc="Matrix rank computation methods
"]


// std imports

// local imports
use srmatrix::api::*;


/// Computes the rank of a matrix using elementary column operations
pub fn rank_eco(a : & MatrixF64) -> usize {
    let mut a = a.clone();
    let mut rank  = 0;
    let m = a.num_rows();
    let n = a.num_cols();
    // forward elimination row wise
    for k in 0..m{
        println!("a: {}", a);
        let (_, cc) = a.max_abs_scalar_in_row(k, k, n);
        if cc > k {
            // TODO : we can switch only part of column
            a.eco_switch(k, cc);
        }
        let mut v = a.view(k, k, m - k, n - k);
        // Pick the pivot
        let pivot  = v.get(0, 0).unwrap();
        if pivot == 0. {
            // Nothing to be done.
            continue;
        }
        // We have a non-zero pivot
        rank += 1;
        for c in 1..v.num_cols(){
            let first = v.get(0, c).unwrap();
            let factor = first/pivot;
            v.eco_scale_add(c, 0, -factor);
        }
    }
    rank
}


pub fn rank(a : & MatrixF64) -> usize{
    rank_eco(a)
}


/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/

#[cfg(test)]
mod test{
    use super::*;

    #[test]
    fn test_rank_eco_0(){
        let a = matrix_rw_f64(2, 2, &[
            1., 0.,
            1., 1.]);
        let r = rank(&a);
        assert_eq!(r, 2);
    }

    #[test]
    fn test_rank_eco_1(){
        let a = matrix_rw_f64(4, 4, &[
        4.0, 3.0, 4.0, 4.0,
        4.0, 1.0, 4.0, 2.0,
        1.0, 2.0, 1.0, 4.0,
        4.0, 3.0, 4.0, 1.0
        ]);
        let r = rank(&a);
        assert_eq!(r, 3);
    }

    #[test]
    fn test_rank_eco_2(){
        let a = matrix_rw_f64(3, 5, &[
        2.0, 1.0, 3.0, 2.0, 1.0,
        4.0, 2.0, 3.0, 3.0, 1.0,
        4.0, 2.0, 4.0, 3.0, 2.0
        ]);
        let r = rank(&a);
        assert_eq!(r, 3);
    }

    #[test]
    fn test_rank_eco_3(){
        let a = matrix_rw_f64(4, 5, &[
        2.0, 1.0, 3.0, 2.0, 1.0,
        4.0, 2.0, 3.0, 3.0, 1.0,
        4.0, 2.0, 4.0, 3.0, 2.0,
        10.0, 5.0, 9.0, 8.0, 3.0
        ]);
        let r = rank(&a);
        assert_eq!(r, 3);
    }

    #[test]
    fn test_rank_eco_4(){
        let a = matrix_rw_f64(2, 3, &[
        1.0000000000000, 1.0000000000000, 1.0000000000000,
        2.0000000000000, 2.0000000000000, 2.0000000000002
        ]);
        let r = rank(&a);
        assert_eq!(r, 2);
    }

    #[test]
    fn test_rank_eco_hilbert(){
        for i in 4..50{
            let m = hilbert(i);
            let r = rank(&m);
            assert_eq!(r, i);
        }
    }
}


