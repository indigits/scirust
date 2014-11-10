#![doc="LU factorization algorithms
"]



// std imports


// local imports
use matrix::*;


#[doc="LU factorization with partial
pivoting problem setup: PA  = LU
"]
pub struct LUPartialPivot<'a> {
    /// The matrix whose LU factorization is to be computed
    pub a : &'a MatrixF64,
    /// The corresponding permutation matrix
    pub p : MatrixF64,
    /// The corresponding lower triangular matrix
    pub l : MatrixF64,
    /// The corresponding upper triangular matrix
    pub u : MatrixF64,
}


impl<'a> LUPartialPivot<'a>{

    /// Setup of a new LU factorization with partial pivot problem
    pub fn new(a : &'a MatrixF64) -> LUPartialPivot<'a>{
        // We support only square matrices.
        assert!(a.is_square());
        let n = a.num_rows();
        LUPartialPivot{a : a, p : Matrix::identity(n, n),
            l :  Matrix::identity(n, n),
            u :  a.clone()
        }
    }

    /// Performs LU factorization with partial pivoting
    pub fn solve(&mut self){
        let n = self.a.num_rows();
        for k in range(0, n){
            // We are working on k-th column.
            let rr = {
                // Create a view of the remaining elements in column
                let col_k_remaining = self.u.view(k, k, n - k, 1);
                //println!("k={}, col_k_remaining: {}", k, col_k_remaining);
                // find the maximum value in this view
                let (_, rr, _) = col_k_remaining.max_abs_scalar();
                // translate rr to the overall row number
                rr + k
            };
            if rr > k {
                // We need to exchange rows of the submatrix.
                self.u.ero_switch(k, rr);
                // The corresponding change in permutation matrix also
                self.p.ero_switch(k, rr);
            }
            // Pick up the pivot
            let pivot = self.u.get(k, k);
            // The lower right part of U matrix
            let mut u_lr  = self.u.view(k + 1, k, n - k - 1, n -k);
            // The lower right part of L matrix
            let mut l_lr = self.l.view(k + 1, k, n - k - 1, n -k);
            for r in range(0, u_lr.num_rows()){
                let first = u_lr.get(r, 0);
                let factor = first  / pivot;
                u_lr.ero_scale_add(r, -1, -factor);
                l_lr.ero_scale_add(r, -1, factor);
            }
        }
    }

    /// Finds the maximum absolute entry in pa - lu
    pub fn max_abs_diff(&self) -> f64 {
        let d  = self.p * *(self.a) - self.l * self.u;
        d.max_abs_scalar_value()
    }


    pub fn print(&self){
        println!("a: {}", self.a);
        println!("p: {}", self.p);
        println!("pa: {}", self.p* *(self.a));
        println!("l: {}", self.l);
        println!("u: {}", self.u);
        println!("lu: {}", self.l * self.u);
    }
}


#[cfg(test)]
mod test{

    use super::*;
    use matrix::*;

    #[test]
    fn test_lu_pp_0(){
        let a = matrix_rw_f64(2, 2, [
            1., 2.,
            3., 8.]);
        let mut lu = LUPartialPivot::new(&a);
        lu.solve();
        lu.print();
        assert_eq!(lu.max_abs_diff(), 0.);
    }

    #[test]
    fn test_lu_pp_1(){
        let a = matrix_rw_f64(3, 3, [
            1., 1., 1.,
            1., 2., 2.,
            1., 2., 3.
            ]);
        let mut lu = LUPartialPivot::new(&a);
        lu.solve();
        lu.print();
        assert_eq!(lu.max_abs_diff(), 0.);
    }


    #[test]
    fn test_lu_pp_tridiag(){
        // tri-diagonal matrix
        let a = matrix_rw_f64(4, 4, [
             1., -1.,  0., 0.,
            -1.,  2., -1., 0.,
             0., -1.,  2., -1.,
             0.,  0., -1., 2.,
            ]);
        let mut lu = LUPartialPivot::new(&a);
        lu.solve();
        lu.print();
        assert_eq!(lu.max_abs_diff(), 0.);
    }
}