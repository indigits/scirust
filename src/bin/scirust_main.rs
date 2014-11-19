#![feature(globs)]
extern crate scirust;
#[cfg(not(test))]
fn main(){
    use scirust::api::*;
    let a = matrix_cw_f64(2,2, [1., 4., 2., 8.].as_slice());
    println!("{}", a);
    let b = vector_f64([3.0, 6.0].as_slice());
    let result = GaussElimination::new(&a, &b).solve();
    match result {
        Ok(x) => println!("{}", x),
        Err(e) => println!("{}", e),
    }
    let m = matrix_cw_c64(2,2, [].as_slice());
    println!("{}", m);
}
