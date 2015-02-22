// std imports
use std::fmt;


pub fn print_slice<T : fmt::Display>(slice : &[T]){
    for v in slice.iter(){
        print!("{} ", v);
    }
    println!("\n");
}
