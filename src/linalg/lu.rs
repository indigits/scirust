#![doc="LU factorization algorithms
"]



// std imports


// local imports
//use error::*;
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

    /// Performs LU factorization with partial pivoting using row operations
    pub fn solve_ero(&mut self){
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


    /// Performs LU factorization with partial pivoting using column operations
    pub fn solve_eco(&mut self){
        let a = &mut self.a;
        let p = &mut self.perm_vector;
        let n = a.num_cols();
        // Iterate over rows
        for k in range(0, n){
            //println!("lu_eco: {}", a);
            // We are working on k-th row.
            let cc = {
                // Create a view of the remaining elements in row
                let row_k_remaining = a.view(k, k, 1, n - k);
                // find the maximum value in this view
                let (_, _, cc) = row_k_remaining.max_abs_scalar();
                // translate cc to the overall row number
                cc + k
            };
            if cc > k {
                // We need to exchange columns of the submatrix.
                let mut l_tr = a.view(k, k, n - k, n - k);
                l_tr.eco_switch(0, cc - k);
                // We will switch only those rows in u which have been filled up so far.
                let mut u_bl = a.view(0, 0, k, n);
                u_bl.eco_switch(k , cc);
                // The corresponding change in permutation matrix also
                p.ero_switch(k, cc);
                //println!("lu_eco: {}", a);
            }
            // Pick up the pivot
            let pivot = a.get(k, k);
            if pivot == 0. {
                continue;
            }
            // The lower right part of L matrix
            let mut l_tr  = a.view(k, k + 1, n - k, n -k - 1);
            // The lower right part of U matrix
            let mut u_bl = a.view(k, k + 1, 1, n -k -1);
            //let mut l_lr = self.l.view(k + 1, k, n - k - 1, n -k);
            for c in range(0, l_tr.num_cols()){
                let first = l_tr.get(0, c);
                let factor = first  / pivot;
                l_tr.eco_scale_add(c, -1, -factor);
                u_bl.set(0, c, factor); 
            }
        }
    }

    /// Finds the maximum absolute entry in pa - lu
    pub fn max_abs_diff_ero(&self, a : &MatrixF64) -> f64 {
        let d  = self.p_ero() * *a - self.l_ero() * self.u_ero();
        d.max_abs_scalar_value()
    }

    /// Finds the maximum absolute entry in ap - lu
    pub fn max_abs_diff_eco(&self, a : &MatrixF64) -> f64 {
        let d  = *a * self.p_eco()  - self.l_eco() * self.u_eco();
        d.max_abs_scalar_value()
    }


    pub fn l_ero(&self) -> MatrixF64 {
        let mut l = self.a.lt();
        l.set_diagonal(1.);
        l
    }
    pub fn u_ero(&self) -> MatrixF64 {
        self.a.ut()
    }
    pub fn p_ero(&self) -> MatrixF64{
        let pv = &self.perm_vector;
        let n = pv.num_cells();
        let mut p : MatrixF64 = Matrix::zeros(n, n);
        for i in range(0, n){
            let index = pv.get(i, 0);
            p.set(i, index as uint, 1.);
        }
        p
    }

    pub fn l_eco(&self) -> MatrixF64 {
        self.a.lt()
    }
    pub fn u_eco(&self) -> MatrixF64 {
        let mut u = self.a.ut();
        u.set_diagonal(1.);
        u
    }
    pub fn p_eco(&self) -> MatrixF64{
        let pv = &self.perm_vector;
        let n = pv.num_cells();
        let mut p : MatrixF64 = Matrix::zeros(n, n);
        for i in range(0, n){
            let index = pv.get(i, 0);
            p.set(index as uint, i, 1.);
        }
        p
    }

    pub fn print_ero(&self){
        println!("p: {}", self.perm_vector);
        let l = self.l_ero();
        let u = self.u_ero();
        println!("l: {}", l);
        println!("u: {}", u);
        println!("lu: {}", l * u);
    }

    pub fn print_eco(&self){
        println!("p: {}", self.perm_vector);
        let l = self.l_eco();
        let u = self.u_eco();
        println!("u: {}", u);
        println!("l: {}", l);
        println!("lu: {}", l * u);
    }
}

///Performs LU factorization  A = LU 
pub fn lu_ero(a : &MatrixF64) -> (MatrixF64, MatrixF64){
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_ero();
        (lu.p_ero().transpose() * lu.l_ero(), lu.u_ero())
}

///Performs LU factorization  PA = LU 
pub fn lup_ero(a : &MatrixF64) -> (MatrixF64, MatrixF64, MatrixF64){
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_ero();
        (lu.l_ero(), lu.u_ero(), lu.p_ero())
}

///Performs LU factorization  A = LU 
pub fn lu_eco(a : &MatrixF64) -> (MatrixF64, MatrixF64){
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_eco();
        (lu.l_eco(), lu.u_eco() * lu.p_eco().transpose())
}

///Performs LU factorization  AP = LU 
pub fn lup_eco(a : &MatrixF64) -> (MatrixF64, MatrixF64, MatrixF64){
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_eco();
        (lu.l_eco(), lu.u_eco(), lu.p_eco())
}


#[cfg(test)]
mod test{
    extern crate test;
    use super::*;
    use matrix::*;
    use self::test::Bencher;

    #[test]
    fn test_lu_ero_0(){
        let a = matrix_rw_f64(2, 2, [
            1., 2.,
            3., 8.]);
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_ero();
        lu.print_ero();
        assert_eq!(lu.max_abs_diff_ero(&a), 0.);


    }

    #[test]
    fn test_lu_eco_0(){
        let a = matrix_rw_f64(2, 2, [
            1., 2.,
            3., 8.]);
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_eco();
        lu.print_eco();
        assert_eq!(lu.max_abs_diff_eco(&a), 0.);
    }


    #[test]
    fn test_lu_ero_1(){
        let a = matrix_rw_f64(3, 3, [
            1., 1., 1.,
            1., 2., 2.,
            1., 2., 3.
            ]);
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_ero();
        lu.print_ero();
        assert_eq!(lu.max_abs_diff_ero(&a), 0.);
    }

    #[test]
    fn test_lu_eco_1(){
        let a = matrix_rw_f64(3, 3, [
            1., 1., 1.,
            1., 2., 2.,
            1., 2., 3.
            ]);
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_eco();
        lu.print_eco();
        assert_eq!(lu.max_abs_diff_eco(&a), 0.);
    }


    #[test]
    fn test_lu_ero_tridiag(){
        // tri-diagonal matrix
        let a = matrix_rw_f64(4, 4, [
             1., -1.,  0., 0.,
            -1.,  2., -1., 0.,
             0., -1.,  2., -1.,
             0.,  0., -1., 2.,
            ]);
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_ero();
        lu.print_ero();
        assert_eq!(lu.max_abs_diff_ero(&a), 0.);
        //assert!(false);
    }

    #[test]
    fn test_lu_eco_tridiag(){
        // tri-diagonal matrix
        let a = matrix_rw_f64(4, 4, [
             1., -1.,  0., 0.,
            -1.,  2., -1., 0.,
             0., -1.,  2., -1.,
             0.,  0., -1., 2.,
            ]);
        let mut lu = LUPartialPivot::new(a.clone());
        lu.solve_eco();
        lu.print_eco();
        assert_eq!(lu.max_abs_diff_eco(&a), 0.);
        //assert!(false);
    }

    #[test]
    fn test_lu_ero_wrapper_method_0(){
        let a = from_range_rw_f64(4,4, 1., 100.);
        let mut lus = LUPartialPivot::new(a.clone());
        lus.solve_ero();
        lus.print_ero();
        assert_eq!(lus.max_abs_diff_ero(&a), 0.);
        let (l, u) = lu_ero(&a);
        assert_eq!(a, l*u);

        let (l, u, p) = lup_ero(&a);
        assert_eq!(p*a, l*u);
    }

    #[test]
    fn test_lu_eco_wrapper_method_0(){
        let a = from_range_rw_f64(4,4, 1., 100.);
        let mut lus = LUPartialPivot::new(a.clone());
        lus.solve_eco();
        lus.print_eco();
        assert_eq!(lus.max_abs_diff_eco(&a), 0.);
        let (l, u) = lu_eco(&a);
        assert_eq!(a, l*u);

        let (l, u, p) = lup_eco(&a);
        assert_eq!(a*p, l*u);
    }

    #[test]
    //#[ignore]
    fn test_lu_ero_hadamard(){
        let a = hadamard(8).unwrap();
        let mut lus = LUPartialPivot::new(a.clone());
        lus.solve_ero();
        assert_eq!(lus.max_abs_diff_ero(&a), 0.);
        let (l, u, p) = lup_ero(&a);
        assert_eq!(p*a, l*u);
    }

    #[test]
    fn test_lu_eco_hadamard(){
        let a = hadamard(16).unwrap();
        let mut lus = LUPartialPivot::new(a.clone());
        lus.solve_eco();
        assert_eq!(lus.max_abs_diff_eco(&a), 0.);
        let (l, u, p) = lup_eco(&a);
        assert_eq!(a*p, l*u);
    }

    #[bench]
    fn bench_lu_ero_hadamard_32 (b: &mut Bencher){
        let a = hadamard(32).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_ero();
        });
    }


    #[bench]
    fn bench_lu_ero_hadamard_64 (b: &mut Bencher){
        let a = hadamard(64).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_ero();
        });
    }

    #[bench]
    fn bench_lu_ero_hadamard_128 (b: &mut Bencher){
        let a = hadamard(128).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_ero();
        });
    }

    #[bench]
    fn bench_lu_ero_hadamard_256 (b: &mut Bencher){
        let a = hadamard(256).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_ero();
        });
    }

    #[bench]
    #[ignore]
    fn bench_lu_ero_hadamard_512 (b: &mut Bencher){
        let a = hadamard(512).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_ero();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_32 (b: &mut Bencher){
        let a = hadamard(32).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_eco();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_64 (b: &mut Bencher){
        let a = hadamard(64).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_eco();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_128 (b: &mut Bencher){
        let a = hadamard(128).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_eco();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_256 (b: &mut Bencher){
        let a = hadamard(256).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_eco();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_512 (b: &mut Bencher){
        let a = hadamard(512).unwrap();
        let mut lus = LUPartialPivot::new(a);
        b.iter(|| {
            lus.solve_eco();
        });
    }
}