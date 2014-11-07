// std imports

// srmat imports

use matrix::MatrixF64;
use linalg::error::*;

/// A Gauss elimination problem specification
pub struct GaussElimination<'a, 'b>{
    pub a : &'a MatrixF64,
    pub b : &'b MatrixF64
}

impl<'a, 'b> GaussElimination<'a, 'b>{

    /// Setup of a new Gauss elimination problem.
    pub fn new(a : &'a MatrixF64, b : &'b MatrixF64) -> GaussElimination<'a, 'b>{
        assert!(a.is_square());
        assert_eq!(a.num_rows(), b.num_rows());
        GaussElimination{a : a , b : b}
    } 

    /// Carries out the procedure of Gauss elimination.
    pub fn solve(&self) -> Result<MatrixF64,LinearSystemError> {
        let mut m = self.a.clone();
        m.append_columns(self.b);
        let rows = m.num_rows();
        let cols = m.num_cols();
        // a vector to hold the positions
        //let mut v  = from_range_uint(rows, 1, 0, rows);
        // Forward elimination process.
        for k in range(0, rows){
            // We are working on k-th column.
            let rr = {
                // Create a view of the remaining elements in column
                let col_k_remaining = m.view(k, k, rows - k, 1);
                //println!("k={}, col_k_remaining: {}", k, col_k_remaining);
                // find the maximum value in this view
                let (_, rr, _) = col_k_remaining.max_abs_scalar();
                // translate rr to the overall row number
                rr + k
            };
            if rr > k {
                // We need to exchange rows of the submatrix.
                m.ero_switch(k, rr);
                // Lets keep this position change information in record
                //v.ero_switch(k, rr);
            }
            // Pick up the pivot
            let pivot = m.get(k, k);
            let mut lower_right  = m.view(k + 1, k, rows - k - 1, cols -k);
            //println!("Pivot: {}", pivot);
            //println!("lower_right: {}", lower_right);
            for r in range(0, lower_right.num_rows()){
                let first = lower_right.get(r, 0);
                let factor = first  / pivot;
                lower_right.ero_scale_add(r, -1, -factor);
            }
            //println!("m: {}", m);
        }
        // Backward substitution starts now.
        let mut b = m.view(0, self.a.num_cols(), 
            self.b.num_rows(), 
            self.b.num_cols());
        let mut r = m.num_rows() - 1;
        loop {
            let pivot = m.get(r, r);
            if pivot == 0. {
                // We have a problem here. We cannot find a solution.
                // TODO: make it more robust for under-determined systems.
                return Err(NoSolution);
            }
            b.ero_scale(r, 1.0/pivot);
            for j in range(r+1, m.num_rows()){
                let factor = m.get(r, j) / pivot;
                b.ero_scale_add(r, j as int, -factor);  
            }
            if r == 0 {
                break;
            }
            r -= 1;
        }
        //println!("m: {}", m);
        // We extract the result.
        Ok(b.to_matrix())
    }
    
}


/// Validates the equality Ax = b subject to maximum
/// absolute error being less than 1e-6.
pub fn validate_ax_eq_b(a : &MatrixF64, x : &MatrixF64, 
    b : &MatrixF64){
   validate_ax_eq_b_eps(a, x, b, 1e-6);
}

/// Validates the equality Ax = b subject to maximum
/// absolute error being less than a specified threshold.
pub fn validate_ax_eq_b_eps(a : &MatrixF64, x : &MatrixF64, 
    b : &MatrixF64, threshold: f64){
    let diff = a * *x  - *b;
    assert!(diff.max_abs_scalar_value() < threshold);
}



#[cfg(test)]
mod test{
    use super::*;
    use matrix::*;


    #[test]
    fn test_ge_0(){
        let a = matrix_f64(2,2, [1., 4., 2., 5.]);
        println!("{}", a);
        let b = vector_f64([3.0, 6.0]);
        let x = GaussElimination::new(&a, &b).solve().unwrap();
        println!("{}", x);
        assert_eq!(x, vector_f64([-1., 2.]));
    }

    #[test]
    fn test_ge_1(){
        let mut a = from_range(3, 3, 1.0, 100.0);
        let x = from_range(3, 1, 1.0, 100.0);
        // a above is rank-2.
        a.set(2,2, 11.0);
        let b  = a * x;
        println!("A: {}", a);
        println!("x: {}", x);
        println!("b: {}", b);
        let ge = GaussElimination::new(&a, &b);
        let z = ge.solve().unwrap();
        println!("z: {}", z);
        /*
        TODO: have better understanding of
        the roundoff error.
        In this case, it is greater than
        1e-15.
        */
        assert!((x - z).max_abs_scalar_value() < 1e-6);
    }  

    #[test]
    fn test_ge_2(){
        let a = matrix_f64(3,3, [1., 0., 0., 
            -1., 1., 0., 
            1., 1., 1.]);
        let b = vector_f64([1., 1., 1.]);
        let ge = GaussElimination::new(&a, &b);
        let x  = ge.solve().unwrap();
        println!("a: {}, x: {}, b: {}", a, x, b);
        // answer is [0, 0, 1]
        validate_ax_eq_b(&a, &x, &b); 
    }

    #[test]
    fn test_ge_3(){
        let a = matrix_f64(3,3, [2., 4., -2., 
            1., -6., 7., 
            1., 0., 2.]);
        let b = vector_f64([5., -2., 9.]);
        let ge = GaussElimination::new(&a, &b);
        let x  = ge.solve().unwrap();
        println!("a: {}, x: {}, b: {}", a, x, b);
        // answer is [1, 1, 2]
        //assert!(false);
        validate_ax_eq_b(&a, &x, &b); 
    }

    #[test]
    fn test_ge_no_solution(){
        let a = matrix_f64(2,2, [1., 4., 2., 8.]);
        println!("{}", a);
        let b = vector_f64([3.0, 6.0]);
        let result = GaussElimination::new(&a, &b).solve();
        match result {
            Ok(x) => {
                println!("{}", x);
                assert!(false);
            },
            Err(e) => println!("{}", e),
        }
        
    }  
}
