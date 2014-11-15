#![doc="LU factorization algorithms
"]



// std imports


// local imports
//use error::*;
use matrix::*;
use error::*;


#[doc="LU factorization with partial
pivoting problem setup: PA  = LU
"]
pub struct LUDecomposition {
    /// The matrix whose LU factorization is to be computed
    ///  The factorization is done in place.
    a : MatrixF64,
    /// The corresponding permutation vector
    perm_vector : MatrixU16,
    /// The corresponding diagonal vector
    diag_vector : MatrixF64,
    /// Indicates if the permutation matrix is pre or post multiplied
    pre : bool
}


impl LUDecomposition{

    /// Setup of a new LU factorization with partial pivot problem
    pub fn new(a : MatrixF64) -> LUDecomposition {
        // We support only square matrices.
        assert!(a.is_square());
        let n = a.num_rows();
        LUDecomposition{a : a,
            perm_vector : from_range_rw_u16(n, 1, 0, n as u16),
            pre : false,
            diag_vector : Matrix::zeros(n, 1)
        }
    }

    /// Performs LU factorization with partial pivoting using row operations
    pub fn decompose_ero(&mut self){
        let a = &mut self.a;
        let p = &mut self.perm_vector;
        let d = &mut self.diag_vector;
        self.pre = true;
        let n = a.num_rows();
        for k in range(0, n){
            // We are working on k-th column.
            let (_, rr) = a.max_abs_scalar_in_col(k, k, n);
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
            // Put it in the diagonal vector
            d.set(k, 0, pivot);
            // The lower right part of U matrix
            let mut u_br  = a.view(k, k, n - k, n -k);
            u_br.ero_scale(0, 1./ pivot);
            // The lower left part of L matrix
            let mut l_bl = a.view(k+1, k, n -k -1, 1);
            for r in range(1, u_br.num_rows()){
                let first = u_br.get(r, 0);
                u_br.ero_scale_add(r, 0, -first);
                l_bl.set(r - 1, 0, first / pivot);
            }
        }
    }


    /// Performs LU factorization with partial pivoting using column operations
    pub fn decompose_eco(&mut self){
        let a = &mut self.a;
        let p = &mut self.perm_vector;
        let d = &mut self.diag_vector;
        self.pre = false;
        let n = a.num_cols();
        // Iterate over rows
        for k in range(0, n){
            // We are working on k-th row.
            // Find the pivot position in the row
            let (_, cc) = a.max_abs_scalar_in_row(k, k, n);
            if cc > k {
                // We need to exchange columns of the submatrix.
                let mut l_tr = a.view(k, k, n - k, n - k);
                l_tr.eco_switch(0, cc - k);
                // We will switch only those rows in u which have been filled up so far.
                let mut u_bl = a.view(0, 0, k, n);
                u_bl.eco_switch(k , cc);
                // The corresponding change in permutation matrix also
                p.ero_switch(k, cc);
            }
            // The top right part of L matrix
            let mut l_tr  = a.view(k, k, n - k, n -k);
            // Pick up the pivot
            let pivot = l_tr.get(0, 0);
            if pivot == 0. {
                continue;
            }
            // Put it in the diagonal vector
            d.set(k, 0, pivot);
            // bring 1 in the diagonal 
            l_tr.eco_scale(0, 1./pivot);
            // The lower right part of U matrix
            let mut u_bl = a.view(k, k + 1, 1, n -k -1);
            for c in range(1, l_tr.num_cols()){
                let first = l_tr.get(0, c);
                let factor = first  / pivot;
                l_tr.eco_scale_add(c, 0, -first);
                u_bl.set(0, c-1, factor); 
            }
        }
    }


    /// Performs LU factorization with partial pivoting using 
    /// Crout's algorithm
    pub fn decompose_crout(&mut self) -> Result<(), SRError>{
        let a = &mut self.a;
        let p = &mut self.perm_vector;
        let n = a.num_cols();
        // Iterate
        for p in range(0, n){
            for r in range(p, n){
                // We are computing l(r, p)
                let mut v = a.get(r, p);
                // subtract l(r, k) * u (k, p)
                for k in range(0, p){
                    v = v - a.get(r, k) * a.get(k, p);
                }
                a.set(r, p, v);
            }
            for c in range(p + 1, n){
                // u(p, p) = 1. Hence we don't compute it.
                // We are computing u(p, c)
                let mut v = a.get(p, c);
                for k in range(0, p){
                    // subtract l(p, k) * u (k, c)
                    v = v - a.get(p, k) * a.get(k, c);
                }
                // divide by l (p, p)
                v = v / a.get(p, p);
                a.set(p, c, v);
            }
        }
        /****
            Following example shows how the computations proceed
            for a 4x4 matrix.

            // lt 0-th column
            l(0, 0) = a(0, 0)
            l(1, 0) = a(1, 0)
            l(2, 0) = a(2, 0)
            l(3, 0) = a(3, 0)
            // ut 0-th row
            u(0, 0) = 1
            u(0, 1) = a(0, 1) / l(0, 0)
            u(0, 2) = a(0, 2) / l(0, 0)
            u(0, 3) = a(0, 3) / l(0, 0)
            // lt 1-st column
            l(1, 1) = a(1, 1) - l(1, 0) * u (0, 1)
            l(2, 1) = a(2, 1) - l(2, 0) * u (0, 1)
            l(3, 1) = a(3, 1) - l(3, 0) * u (0, 1)
            // ut 1-st row
            u(1, 1) = 1
            u(1, 2) = ( a(1, 2) - l(1, 0) * u(0, 2) ) / l (1, 1)  
            u(1, 3) = ( a(1, 2) - l(1, 0) * u(0, 3) ) / l (1, 1)


         */
        Ok(())
    }

    /// Finds the maximum absolute entry in a - ldu
    pub fn max_abs_diff(&self, a : &MatrixF64) -> f64 {
        let d = self.d();
        let p = self.p();
        let l = self.l();
        let u = self.u();
        let a = if self.pre {
            p * *a
        }
        else{
            *a * p
        };
        println!("p{}", a);
        println!("p a {}", a);
        let b = l * d * u;
        println!("l d u {}", b);
        let diff  = a  - b;
        diff.max_abs_scalar_value()
    }

    pub fn l(&self) -> MatrixF64 {
        self.a.lt()
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
            if self.pre {
                p.set(i, index as uint, 1.);
            }
            else{
              p.set(index as uint, i, 1.);  
            }
        }
        p
    }

    pub fn d(&self) -> MatrixF64{
        Matrix::diag_from_vec(&self.diag_vector)
    }

    pub fn print(&self){
        println!("p: {}", self.perm_vector);
        let d = self.d();
        let l = self.l();
        let u = self.u();
        println!("l: {}", l);
        println!("d: {}", d);
        println!("u: {}", u);
        println!("lu: {}", l * d * u);
    }
}

///Performs LU factorization  A = LU 
pub fn lu_ero(a : &MatrixF64) -> (MatrixF64, MatrixF64){
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_ero();
        (lu.p().transpose() * lu.l(), lu.d() * lu.u())
}

///Performs LU factorization  PA = LU 
pub fn lup_ero(a : &MatrixF64) -> (MatrixF64, MatrixF64, MatrixF64){
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_ero();
        (lu.l(), lu.d() * lu.u(), lu.p())
}

///Performs LU factorization  A = LU 
pub fn lu_eco(a : &MatrixF64) -> (MatrixF64, MatrixF64){
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_eco();
        (lu.l() * lu.d(), lu.u() * lu.p().transpose())
}

///Performs LU factorization  AP = LU 
pub fn lup_eco(a : &MatrixF64) -> (MatrixF64, MatrixF64, MatrixF64){
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_eco();
        (lu.l() * lu.d(), lu.u(), lu.p())
}


/******************************************************
 *
 *   Unit tests follow.
 *
 *******************************************************/
#[cfg(test)]
mod test{
    use super::*;
    use matrix::*;

    #[test]
    fn test_lu_ero_0(){
        let a = matrix_rw_f64(2, 2, [
            1., 2.,
            3., 8.]);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_ero();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);


    }

    #[test]
    fn test_lu_eco_0(){
        let a = matrix_rw_f64(2, 2, [
            1., 2.,
            3., 8.]);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_eco();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
    }


    #[test]
    fn test_lu_ero_1(){
        let a = matrix_rw_f64(3, 3, [
            1., 1., 1.,
            1., 2., 2.,
            1., 2., 3.
            ]);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_ero();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
    }

    #[test]
    fn test_lu_eco_1(){
        let a = matrix_rw_f64(3, 3, [
            1., 1., 1.,
            1., 2., 2.,
            1., 2., 3.
            ]);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_eco();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
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
        println!("a {}", a);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_ero();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
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
        println!("a {}", a);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_eco();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
        //assert!(false);
    }

    #[test]
    fn test_lu_ero_wrapper_method_0(){
        let a = from_range_rw_f64(4,4, 1., 100.);
        println!("a {}", a);
        let mut lus = LUDecomposition::new(a.clone());
        lus.decompose_ero();
        lus.print();
        assert!(lus.max_abs_diff(&a) <  1e-7);
        let (l, u) = lu_ero(&a);
        let d = a - l * u;
        assert!(d.max_abs_scalar_value() < 1e-7);

        let (l, u, p) = lup_ero(&a);
        let d = p*a -  l*u;
        assert!(d.max_abs_scalar_value() < 1e-7);
    }

    #[test]
    fn test_lu_eco_wrapper_method_0(){
        let a = from_range_rw_f64(4,4, 1., 100.);
        let mut lus = LUDecomposition::new(a.clone());
        lus.decompose_eco();
        lus.print();
        assert_eq!(lus.max_abs_diff(&a), 0.);
        let (l, u) = lu_eco(&a);
        assert_eq!(a, l*u);

        let (l, u, p) = lup_eco(&a);
        assert_eq!(a*p, l*u);
    }

    #[test]
    //#[ignore]
    fn test_lu_ero_hadamard(){
        let a = hadamard(8).unwrap();
        let mut lus = LUDecomposition::new(a.clone());
        lus.decompose_ero();
        assert_eq!(lus.max_abs_diff(&a), 0.);
        let (l, u, p) = lup_ero(&a);
        assert_eq!(p*a, l*u);
    }

    #[test]
    fn test_lu_eco_hadamard(){
        let a = hadamard(16).unwrap();
        let mut lus = LUDecomposition::new(a.clone());
        lus.decompose_eco();
        assert_eq!(lus.max_abs_diff(&a), 0.);
        let (l, u, p) = lup_eco(&a);
        assert_eq!(a*p, l*u);
    }

    #[test]
    #[ignore]
    fn test_lu_crout_1(){
        let a = matrix_rw_f64(3, 3, [
            1., 1., 1.,
            1., 2., 2.,
            1., 2., 3.
            ]);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_crout();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
    }

    #[test]
    #[ignore]
    fn test_lu_crout_2(){
        let a = matrix_rw_f64(3, 3, [
            25., 5., 1.,
            64., 8., 1.,
            144., 12., 1.
            ]);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_crout();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 20.);
    }


    #[test]
    #[ignore]
    fn test_lu_crout_tridiag(){
        // tri-diagonal matrix
        let a = matrix_rw_f64(4, 4, [
             1., -1.,  0., 0.,
            -1.,  2., -1., 0.,
             0., -1.,  2., -1.,
             0.,  0., -1., 2.,
            ]);
        let mut lu = LUDecomposition::new(a.clone());
        lu.decompose_crout();
        lu.print();
        assert_eq!(lu.max_abs_diff(&a), 0.);
        //assert!(false);
    }

    #[test]
    #[ignore]
    fn test_lu_crout_hadamard(){
        let a = hadamard(16).unwrap();
        let mut lus = LUDecomposition::new(a.clone());
        lus.decompose_crout();
        assert_eq!(lus.max_abs_diff(&a), 0.);
    }

}
/******************************************************
 *
 *   Bench marks follow.
 *
 *******************************************************/
#[cfg(test)]
mod bench {
    extern crate test;
    use super::*;
    use matrix::*;
    use self::test::Bencher;

    #[bench]
    fn bench_lu_ero_hadamard_32 (b: &mut Bencher){
        let a = hadamard(32).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_ero();
        });
    }


    #[bench]
    fn bench_lu_ero_hadamard_64 (b: &mut Bencher){
        let a = hadamard(64).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_ero();
        });
    }

    #[bench]
    fn bench_lu_ero_hadamard_128 (b: &mut Bencher){
        let a = hadamard(128).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_ero();
        });
    }

    #[bench]
    fn bench_lu_ero_hadamard_256 (b: &mut Bencher){
        let a = hadamard(256).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_ero();
        });
    }

    #[bench]
    #[ignore]
    fn bench_lu_ero_hadamard_512 (b: &mut Bencher){
        let a = hadamard(512).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_ero();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_32 (b: &mut Bencher){
        let a = hadamard(32).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_eco();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_64 (b: &mut Bencher){
        let a = hadamard(64).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_eco();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_128 (b: &mut Bencher){
        let a = hadamard(128).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_eco();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_256 (b: &mut Bencher){
        let a = hadamard(256).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_eco();
        });
    }

    #[bench]
    fn bench_lu_eco_hadamard_512 (b: &mut Bencher){
        let a = hadamard(512).unwrap();
        let mut lus = LUDecomposition::new(a);
        b.iter(|| {
            lus.decompose_eco();
        });
    }
}