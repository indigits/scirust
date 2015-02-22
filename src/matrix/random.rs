// library imports
use rand;
use rand::Rng;
use rand::distributions::normal::{StandardNormal};
use super::matrix::*;
use matrix::traits::*;


#[allow(deprecated)]
/// Generate a random matrix of uniformly distributed numbers
pub fn rand_std_normal(rows: usize, cols : usize)-> MatrixF64 {
    let mut m : Matrix<f64> = Matrix::new(rows, cols);
    let mut rng = rand::thread_rng();
    for c in range(0, cols){
        for r in range (0, rows){
            let StandardNormal(n) = rng.gen::<StandardNormal>();
            m.set(r, c, n);
        }
    }
    m
}



/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/


#[cfg(test)]
#[allow(unused_imports)]
mod tests {

    use api::*;

    #[test]
    fn test_rand_std_normal(){
        let m1  = super::rand_std_normal(10, 10);
        assert_eq!(m1.num_cells(), 100);
    }

}