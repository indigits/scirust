#![doc="LU factorization algorithms
"]



// std imports


// local imports
use matrix::*;


#[doc="LU factorization with partial
pivoting problem setup: PA  = LU
"]
pub struct LUPartialPivot {
    /// The matrix whose LU factorization is to be computed
    ///  The factorization is done in place.
    a : MatrixF64,
    /// The corresponding permutation vector
    perm_vector : MatrixU16,
}


impl LUPartialPivot{

    /// Setup of a new LU factorization with partial pivot problem
    pub fn new(a : MatrixF64) -> LUPartialPivot {
        // We support only square matrices.
        assert!(a.is_square());
        let n = a.num_rows();
        LUPartialPivot{a : a,
            perm_vector : from_range_rw_u16(n, 1, 0, n as u16)
        }
    }

    /// Performs LU factorization with partial pivoting
    pub fn solve(&mut self){
        let a = &mut self.a;
        let p = &mut self.perm_vector;
        let n = a.num_rows();
        for k in range(0, n){
            // We are working on k-th column.
            let rr = {
                // Create a view of the remaining elements in column
                let col_k_remaining = a.view(k, k, n - k, 1);
                // find the maximum value in this view
                let (_, rr, _) = col_k_remaining.max_abs_scalar();
                // translate rr to the overall row number
                rr + k
            };
            if rr > k {
                // We need to exchange rows of the submatrix.
                let mut u_br = a.view(k, k, n - k, n - k);
                u_br.ero_switch(0, rr - k);
                // We will switch only those columns in l which have been filled up so far.
                let mut l_tl = a.view(0, 0, n, k);
                l_tl.ero_switch(k , rr);
                // The corresponding change in permutation matrix also
                p.ero_switch(k, rr);
            }
            // Pick up the pivot
            let pivot = a.get(k, k);
            // The lower right part of U matrix
            let mut u_br  = a.view(k + 1, k, n - k - 1, n -k);
            let mut l_bl = a.view(k+1, k, n -k -1, 1);
            // The lower right part of L matrix
            //let mut l_lr = self.l.view(k + 1, k, n - k - 1, n -k);
            for r in range(0, u_br.num_rows()){
                let first = u_br.get(r, 0);
                let factor = first  / pivot;
                u_br.ero_scale_add(r, -1, -factor);
                l_bl.set(r, 0, factor); 
            }
        }
    }

    /// Finds the maximum absolute entry in pa - lu
    pub fn max_abs_diff(&self, a : &MatrixF64) -> f64 {
        let d  = self.p() * *a - self.l() * self.u();
        d.max_abs_scalar_value()
    }


    pub fn l(&self) -> MatrixF64 {
        let mut l = self.a.lt();
        l.set_diagonal(1.);
        l
    }
    pub fn u(&self) -> MatrixF64 {
        self.a.ut()
    }
    pub fn p(&self) -> MatrixF64{
        let pv = &self.perm_vector;
        let n = pv.num_cells();
        let mut p : MatrixF64 = Matrix::zeros(n, n);
        for i in range(0, n){
            let index = pv.get(i, 0);
            p.set(i, index as uint, 1.);
        }
        p
    }


    pub fn print(&self){
        println!("p: {}", self.perm_vector);
        let l = self.l();
        let u = self.u();
        println!("l: {}", l);
        println!("u: {}", u);
        println!("lu: {}", l * u);
    }
}

/// Wrapper function to perform LU factorization
pub fn lu(a : &MatrixF64) -> (MatrixF64, MatrixF64){
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve();
        (lu.p().transpose() * lu.l(), lu.u())
}

pub fn lup(a : &MatrixF64) -> (MatrixF64, MatrixF64, MatrixF64){
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve();
        (lu.l(), lu.u(), lu.p())
}


#[cfg(test)]
mod test{
    extern crate test;
    use super::*;
    use matrix::*;
    use self::test::Bencher;

    #[test]
    fn test_lu_pp_0(){
        let a = matrix_rw_f64(2, 2, [
            1., 2.,
            3., 8.]);
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
    }

    #[test]
    fn test_lu_pp_1(){
        let a = matrix_rw_f64(3, 3, [
            1., 1., 1.,
            1., 2., 2.,
            1., 2., 3.
            ]);
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
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
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
        //assert!(false);
    }

    #[test]
    fn test_lu_pp_wrapper_method_0(){
        let a = from_range_rw_f64(4,4, 1., 100.);
        let mut lus = LUPartialPivot::new(a.clone());
        lus.solve();
        lus.print();
        assert_eq!(lus.max_abs_diff(&a), 0.);
        let (l, u) = lu(&a);
        assert_eq!(a, l*u);

        let (l, u, p) = lup(&a);
        assert_eq!(p*a, l*u);
    }

    #[test]
    //#[ignore]
    fn test_lu_hadamard(){
        let a = hadamard(8).unwrap();
        let mut lus = LUPartialPivot::new(a.clone());
        lus.solve();
        assert_eq!(lus.max_abs_diff(&a), 0.);
        let (l, u, p) = lup(&a);
        assert_eq!(p*a, l*u);
    }

    #[bench]
    fn bench_lu_hadamard_32 (b: &mut Bencher){
        let a = hadamard(32).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve();
        });
    }


    #[bench]
    fn bench_lu_hadamard_64 (b: &mut Bencher){
        let a = hadamard(64).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve();
        });
    }

    #[bench]
    fn bench_lu_hadamard_128 (b: &mut Bencher){
        let a = hadamard(128).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve();
        });
    }

    #[bench]
    fn bench_lu_hadamard_256 (b: &mut Bencher){
        let a = hadamard(256).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve();
        });
    }

    #[bench]
    #[ignore]
    fn bench_lu_hadamard_512 (b: &mut Bencher){
        let a = hadamard(512).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve();
        });
    }
}