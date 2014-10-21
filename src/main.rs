extern crate srmat;


#[cfg(not(test))]
fn main() {
    use srmat::{Mat, MatI64};
    let m : MatI64  = Mat::new(4, 4);
    println!("{}", m);
    let m : MatI64 = Mat::from_iter(4, 4, range(1, 20));
    println!("{}", m);
    let m : MatI64 = Mat::from_iter(4, 4, range(1, 12));
    println!("{}", m);
    println!("{}", m.to_std_vec().as_slice());
    println!("{}", m.get(1,1));
    println!("{}", m.get(2,2));
}
